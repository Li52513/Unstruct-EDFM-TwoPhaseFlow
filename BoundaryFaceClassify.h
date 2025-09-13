#pragma once
#pragma once
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>
#include "Mesh.h"  // 需要 Mesh/Face/Node/Vector 接口

namespace BoundaryFaceClassify
{
    // ---------------------------- 结果类型 ----------------------------
// 统一 2D/3D 的六面分组（2D 时 z0/zL 为空）
    struct FaceGroups
    {
        std::vector<int> x0, xL;  // x=0 / x=Lx
        std::vector<int> y0, yL;  // y=0 / y=Ly
        std::vector<int> z0, zL;  // z=0 / z=Lz（Lz<=0 时为空）
    };

    // ---------------------------- 小工具 ----------------------------
    inline bool nearAbs(double v, double tgt, double tolAbs) noexcept
    {
        return std::abs(v - tgt) < tolAbs;
    }

    // 按几何尺度把相对容差转成绝对容差；L<=0 时给一个保底尺度
    inline double makeAbsTol(double L, double tolRel) noexcept
    {
        const double scale = (L > 0.0 ? L : 1.0);
        return std::max(1e-14, tolRel * scale);
    }

    // 收集“真边界面” id（Face::isBoundary()==true）
    inline std::vector<int> getBoundaryFaceIDs(const Mesh& mesh) {
        std::vector<int> ids;
        const auto& faces = mesh.getFaces();
        ids.reserve(faces.size() / 4u + 8u);
        for (const auto& f : faces) {
            if (f.isBoundary()) ids.push_back(f.id); // 约定 id 为 1-based
        }
        return ids;
    }

    // 从 Face.FaceNodeCoords 计算各轴向 min/max（不依赖节点表索引）
    inline void faceMinMaxXYZ_fromFace(const Face& f,
        double& xmin, double& xmax,
        double& ymin, double& ymax,
        double& zmin, double& zmax)
    {
        xmin = ymin = zmin = std::numeric_limits<double>::infinity();
        xmax = ymax = zmax = -std::numeric_limits<double>::infinity();

        // FaceNodeCoords: std::vector<Vector>（2D:2个点；3D:3或4个点）
        for (const Vector& p : f.FaceNodeCoords) {
            xmin = std::min(xmin, p.m_x); xmax = std::max(xmax, p.m_x);
            ymin = std::min(ymin, p.m_y); ymax = std::max(ymax, p.m_y);
            zmin = std::min(zmin, p.m_z); zmax = std::max(zmax, p.m_z);
        }
    }

    // 判定“整张面”是否位于 x/y/z = 常数 的坐标平面（min/max 同时贴近）
    inline bool faceOnPlaneX(const Face& f, double xplane, double tolAbs)
    {
        double xmin, xmax, ymin, ymax, zmin, zmax;
        faceMinMaxXYZ_fromFace(f, xmin, xmax, ymin, ymax, zmin, zmax);
        return nearAbs(xmin, xplane, tolAbs) && nearAbs(xmax, xplane, tolAbs);
    }
    inline bool faceOnPlaneY(const Face& f, double yplane, double tolAbs)
    {
        double xmin, xmax, ymin, ymax, zmin, zmax;
        faceMinMaxXYZ_fromFace(f, xmin, xmax, ymin, ymax, zmin, zmax);
        return nearAbs(ymin, yplane, tolAbs) && nearAbs(ymax, yplane, tolAbs);
    }
    inline bool faceOnPlaneZ(const Face& f, double zplane, double tolAbs)
    {
        double xmin, xmax, ymin, ymax, zmin, zmax;
        faceMinMaxXYZ_fromFace(f, xmin, xmax, ymin, ymax, zmin, zmax);
        return nearAbs(zmin, zplane, tolAbs) && nearAbs(zmax, zplane, tolAbs);
    }
    // ---------------------------- 主功能：六面分类（2D/3D统一） ----------------------------
/**
 * @brief 把所有“真边界面”按几何位置分到 x0/xL, y0/yL, z0/zL 六类。
 *        若某面同时落在多个平面（角/棱），按法向主导轴唯一分派。
 *
 * @param Lx,Ly,Lz  物理域尺寸；2D 情况传 Lz<=0 即可
 * @param tolRel    相对容差（建议 1e-9 ~ 1e-8）
 */
    inline FaceGroups ClassifyBoundaryFaces(
        const Mesh& mesh,
        double Lx, double Ly, double Lz,
        double tolRel = 1e-9)
    {
        FaceGroups groups;

        const double tolX = makeAbsTol(Lx, tolRel);
        const double tolY = makeAbsTol(Ly, tolRel);
        const double tolZ = makeAbsTol(Lz, tolRel);

        for (const auto& f : mesh.getFaces())
        {
            if (!f.isBoundary()) continue;

            // 是否整体落在某个坐标平面
            const bool onX0 = faceOnPlaneX(f, 0.0, tolX);
            const bool onXL = faceOnPlaneX(f, Lx, tolX);
            const bool onY0 = faceOnPlaneY(f, 0.0, tolY);
            const bool onYL = faceOnPlaneY(f, Ly, tolY);
            const bool onZ0 = (Lz > 0.0) ? faceOnPlaneZ(f, 0.0, tolZ) : false;
            const bool onZL = (Lz > 0.0) ? faceOnPlaneZ(f, Lz, tolZ) : false;

            if (!(onX0 || onXL || onY0 || onYL || onZ0 || onZL)) continue;

            // 角/棱的唯一分派：选择法向绝对值最大的轴
            const double ax = std::abs(f.normal.m_x);
            const double ay = std::abs(f.normal.m_y);
            const double az = std::abs(f.normal.m_z);
            const int fid = f.id;

            if (ax >= ay && ax >= az) {
                if (onX0) groups.x0.push_back(fid);
                else if (onXL) groups.xL.push_back(fid);
                else if (onY0) groups.y0.push_back(fid);
                else if (onYL) groups.yL.push_back(fid);
                else if (onZ0) groups.z0.push_back(fid);
                else if (onZL) groups.zL.push_back(fid);
            }
            else if (ay >= ax && ay >= az) {
                if (onY0) groups.y0.push_back(fid);
                else if (onYL) groups.yL.push_back(fid);
                else if (onX0) groups.x0.push_back(fid);
                else if (onXL) groups.xL.push_back(fid);
                else if (onZ0) groups.z0.push_back(fid);
                else if (onZL) groups.zL.push_back(fid);
            }
            else {
                if (onZ0) groups.z0.push_back(fid);
                else if (onZL) groups.zL.push_back(fid);
                else if (onX0) groups.x0.push_back(fid);
                else if (onXL) groups.xL.push_back(fid);
                else if (onY0) groups.y0.push_back(fid);
                else if (onYL) groups.yL.push_back(fid);
            }
        }
        return groups;
    }

    // 便捷重载：2D 情况（内部把 Lz 置 0）
    inline FaceGroups ClassifyBoundaryFaces(
        const Mesh& mesh,
        double Lx, double Ly,
        double tolRel = 1e-9)
    {
        return ClassifyBoundaryFaces(mesh, Lx, Ly, /*Lz=*/0.0, tolRel);
    }


}




///// 按几何位置把所有边界面分成 bottom/top/left/right 四组
   //struct FaceGroups 
   //{
   //    std::vector<int> bottom;
   //    std::vector<int> top;
   //    std::vector<int> left;
   //    std::vector<int> right;
   //};


   ///// 提取所有真边界面的 ID
   //inline std::vector<int> getBoundaryFaceIDs(const Mesh& mesh)
   //{
   //    std::vector<int> ids;
   //    for (auto const& f : mesh.getFaces())
   //    {
   //        if (f.isBoundary())
   //            ids.push_back(f.id);
   //    }
   //    return ids;
   //}
  ///**
   // * @param mesh      已构建好的 Mesh 对象
   // * @param lengthX   区域在 X 方向的长度（用于判定左/右）
   // * @param lengthY   区域在 Y 方向的长度（用于判定上/下）
   // * @param tol       坐标比较容差，默认 1e-8
   // */
   //inline FaceGroups ClassifySolidMatrixMeshBoundaryFaces(
   //    const Mesh& mesh,
   //    double lengthX,
   //    double lengthY,
   //    double tol = 1e-8)
   //{
   //    FaceGroups groups;
   //    auto bfaceIDs = getBoundaryFaceIDs(mesh);
   //    for (int fid : bfaceIDs) {
   //        const auto& f = mesh.getFaces()[fid - 1]; // face.id 从 1 开始
   //        double x = f.midpoint.m_x;
   //        double y = f.midpoint.m_y;
       //        if (std::abs(y - 0.0) < tol) groups.bottom.push_back(fid);
   //        else if (std::abs(y - lengthY) < tol) groups.top.push_back(fid);
   //        else if (std::abs(x - 0.0) < tol) groups.left.push_back(fid);
   //        else if (std::abs(x - lengthX) < tol) groups.right.push_back(fid);
   //    }
   //    return groups;
   //}