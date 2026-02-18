/**
 * @file Benchmark_BoundaryExport.cpp
 * @brief 边界条件系数导出测试 (Final Production Version)
 * @details
 * 功能描述：
 * 1. 构建 10x10x10 的 3D 基岩网格，自动进行边界分类。
 * 2. 利用 3D_BoundaryConditionManager 配置物理场边界：
 * - 压力：左边界施加线性静水压力梯度 (P = P0 + rho*g*z)。
 * - 温度/饱和度：施加定值或绝热条件。
 * 3. 遍历所有物理边界，安全获取 Face 对象几何中心 (midpoint)。
 * 4. 计算并导出每个边界面的 (Tag, Coord, a, b, c) 系数至 CSV，用于离散化验证。
 */

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <fstream>
#include <map>
#include <cmath>

 // =========================================================
 // 引入核心模块头文件
 // =========================================================
#include "3D_MeshManager.h"
#include "BoundaryConditionManager.h"
#include "MeshDefinitions.h"
#include "SolverContrlStrName_op.h"
// 引入基础数学类型 (Vector定义)
#include "UserDefineVarType.h" 

// 命名空间引用
using namespace BoundarySetting;
using namespace PhysicalProperties_string_op;

// =========================================================
// 辅助结构：用于 CSV 输出的单行数据缓存
// =========================================================
struct BoundaryDataRow {
    int faceID;         // 原始 ID (1-based)
    int internalIndex;  // 内存 Index (0-based, 调试用)
    int tagID;          // 物理边界 Tag
    double x, y, z;     // 面几何中心坐标

    // 压力方程系数 (a*P + b*Flux = c)
    double p_a = 0, p_b = 0, p_c = 0;
    bool p_set = false;

    // 温度方程系数
    double t_a = 0, t_b = 0, t_c = 0;
    bool t_set = false;

    // 饱和度方程系数
    double s_a = 0, s_b = 0, s_c = 0;
    bool s_set = false;
};

// =========================================================
// 核心测试函数实现
// =========================================================
void RunTest_BoundaryCondition_Export()
{
    std::cout << "\n========================================================================" << std::endl;
    std::cout << "   Test E: Boundary Coefficients Export (Production Ready)" << std::endl;
    std::cout << "========================================================================" << std::endl;

    // 1. 环境准备：构建网格
    // ---------------------------------------------------------
    std::cout << " -> [Mesh] Building 10x10x10 Grid for Validation..." << std::endl;
    // 参数说明: Lx=100, Ly=100, Lz=100, Nx=10, Ny=10, Nz=10
    // 使用 Prism=true, QuadBase=false (生成 Prism 单元) 或 根据您的默认设置
    MeshManager_3D meshMgr(100.0, 100.0, 100.0, 2, 2, 2, true, false);

    // 核心步骤：构建拓扑、计算几何信息、执行边界分类
    meshMgr.BuildSolidMatrixGrid_3D(NormalVectorCorrectionMethod::OrthogonalCorrection, "Bench_Export_Geom");

    //输出网格几何信息至txt
    meshMgr.exportMeshInfortoTxt("Test/BoundaryTest/Mesh_Infor");

    // 2. 实例化物理场 BC 管理器
    // ---------------------------------------------------------
    BoundaryConditionManager pBCMgr; // 压力场
    BoundaryConditionManager tBCMgr; // 温度场
    BoundaryConditionManager sBCMgr; // 饱和度场

    // 3. 配置测试场景 (Physical Setup)
    // ---------------------------------------------------------
    std::cout << " -> [Setup] Configuring Physics Boundary Conditions..." << std::endl;

    // --- 左边界 (Inlet, X=0, Tag=LEFT) ---
    // 场景：深部地层，压力随深度 Z 线性变化 (静水压力)
    // P(z) = RefVal + Grad * (z - RefCoord)
    // 设 Z=0 处压力 30MPa，梯度 -1e4 Pa/m (即每升高 1m 压力减小 1e4 Pa)
    pBCMgr.SetLinearDirichletBC(MeshTags::LEFT, 30.0e6, 0.0, -1.0e4, 2); // Axis 2 = Z方向 1 = Y方向 0 = z方向

    // 温度：定温注入 360K
    tBCMgr.SetDirichletBC(MeshTags::LEFT, 360.0);
    // 饱和度：纯水注入 Sw=1.0
    sBCMgr.SetDirichletBC(MeshTags::LEFT, 1.0);

    // --- 右边界 (Outlet, X=L, Tag=RIGHT) ---
    // 压力：定压生产 20MPa
    pBCMgr.SetDirichletBC(MeshTags::RIGHT, 20.0e6);
    // 温度/饱和度：绝热/无通量 (Neumann 0)
    tBCMgr.SetNeumannBC(MeshTags::RIGHT, 0.0);
    sBCMgr.SetNeumannBC(MeshTags::RIGHT, 0.0);

    // --- 其他封闭边界 (Bottom, Top, Front, Back) ---
    // 全部设为无通量 (Neumann 0)
    const std::vector<int> closedTags = {
        MeshTags::BOTTOM, MeshTags::TOP, MeshTags::TAG_FRONT, MeshTags::TAG_BACK
    };

    for (int tag : closedTags) {
        pBCMgr.SetNeumannBC(tag, 0.0);
        tBCMgr.SetNeumannBC(tag, 0.0);
        sBCMgr.SetNeumannBC(tag, 0.0); // 饱和度通常也封闭
    }

    // 4. 数据收集与导出 (Safe Access & Execution)
    // ---------------------------------------------------------
    std::string filename = "Test/BoundaryTest/Boundary_Coefficients_Export.csv";
    std::cout << " -> [Export] Generating " << filename << "..." << std::endl;

    std::ofstream csv(filename);
    if (!csv.is_open()) {
        std::cerr << "[Error] Cannot open file: " << filename << ". Please ensure directory exists." << std::endl;
        return;
    }

    // 4.1 写入 CSV 表头
    csv << "FaceID,InternalIdx,TagID,CenterX,CenterY,CenterZ,"
        << "P_Type,P_a,P_b,P_c,"
        << "T_Type,T_a,T_b,T_c,"
        << "S_Type,S_a,S_b,S_c\n";

    // 4.2 准备遍历数据
    // (A) 获取基于 Tag 的边界分组 (Face ID List)
    //     注意：调用 BoundaryFaceClassify_ByTag.h 中定义的结构
    const auto& groupsStruct = meshMgr.boundaryFaces_byTag();

    // (B) 获取 Mesh 原始数据引用 (用于 ID -> Object 转换)
    const Mesh& mesh = meshMgr.mesh();
    const auto& allFaces = mesh.getFaces();

    // (C) 构建统一遍历 Map (Tag -> ID List Pointer)
    //     这模拟了 Solver 中 "遍历所有已知边界类型" 的逻辑
    std::map<int, const std::vector<int>*> activeTags;
    if (!groupsStruct.x0.empty()) activeTags[MeshTags::LEFT] = &groupsStruct.x0;
    if (!groupsStruct.xL.empty()) activeTags[MeshTags::RIGHT] = &groupsStruct.xL;
    if (!groupsStruct.y0.empty()) activeTags[MeshTags::BOTTOM] = &groupsStruct.y0;
    if (!groupsStruct.yL.empty()) activeTags[MeshTags::TOP] = &groupsStruct.yL;
    if (!groupsStruct.z0.empty()) activeTags[MeshTags::TAG_FRONT] = &groupsStruct.z0;
    if (!groupsStruct.zL.empty()) activeTags[MeshTags::TAG_BACK] = &groupsStruct.zL;
    // 添加自定义边界 (如障碍物等)
    for (const auto& p : groupsStruct.custom) activeTags[p.first] = &p.second;

    int totalRows = 0;

    // 4.3 执行遍历
    for (const auto& pair : activeTags) {
        int tag = pair.first;
        const auto& faceIDs = *pair.second; // 这里存储的是 Face ID (1-based)

        // 遍历当前 Tag 下的所有 Face
        for (int faceID : faceIDs) {

            // [CRITICAL] 安全访问模式: ID -> Index -> Object
            // 严禁假设 index = id - 1，必须查询 Map
            int fIdx = mesh.getFaceIndex(faceID);

            // 越界与有效性检查
            if (fIdx == -1) {
                std::cerr << "[Warning] Face ID " << faceID << " (Tag " << tag << ") not found in Mesh Index Map!" << std::endl;
                continue;
            }
            if (fIdx < 0 || fIdx >= static_cast<int>(allFaces.size())) {
                std::cerr << "[Error] Face Index " << fIdx << " is out of vector bounds!" << std::endl;
                continue;
            }

            // 安全获取 Face 对象引用
            const Face& f = allFaces[fIdx];

            // [Fix] 获取面几何中心
            // 直接使用 face.h 中定义的 midpoint 成员 (在 ComputeGeometricInfor 中已计算)
            Vector center = f.midpoint;

            // 数据行填充
            BoundaryDataRow row;
            row.faceID = faceID;
            row.internalIndex = fIdx;
            row.tagID = tag;
            row.x = center.m_x;
            row.y = center.m_y;
            row.z = center.m_z;

            // 查询各物理场 BC (传入空间坐标以支持线性插值)
            // 压力 P
            if (pBCMgr.HasBC(tag)) {
                auto abc = pBCMgr.GetBCCoefficients(tag, center);
                row.p_set = true;
                row.p_a = abc.a; row.p_b = abc.b; row.p_c = abc.c;
            }
            // 温度 T
            if (tBCMgr.HasBC(tag)) {
                auto abc = tBCMgr.GetBCCoefficients(tag, center);
                row.t_set = true;
                row.t_a = abc.a; row.t_b = abc.b; row.t_c = abc.c;
            }
            // 饱和度 S
            if (sBCMgr.HasBC(tag)) {
                auto abc = sBCMgr.GetBCCoefficients(tag, center);
                row.s_set = true;
                row.s_a = abc.a; row.s_b = abc.b; row.s_c = abc.c;
            }

            // 4.4 写入 CSV 行
            csv << row.faceID << "," << row.internalIndex << "," << row.tagID << ","
                << row.x << "," << row.y << "," << row.z << ",";

            // P Output (Format: Type, a, b, c)
            if (row.p_set) csv << "BC," << row.p_a << "," << row.p_b << "," << row.p_c << ",";
            else          csv << "None,0,0,0,";

            // T Output
            if (row.t_set) csv << "BC," << row.t_a << "," << row.t_b << "," << row.t_c << ",";
            else          csv << "None,0,0,0,";

            // S Output
            if (row.s_set) csv << "BC," << row.s_a << "," << row.s_b << "," << row.s_c;
            else          csv << "None,0,0,0";

            csv << "\n";
            totalRows++;
        }
    }

    csv.close();
    std::cout << " -> [Result] Successfully exported " << totalRows << " boundary faces." << std::endl;
    std::cout << " -> [Verify] Please check '" << filename << "' for linear pressure gradients." << std::endl;
}