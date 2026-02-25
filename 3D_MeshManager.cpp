#include "3D_MeshManager.h"
#include "Mesh.h"
#include "FaceIndexedOctree.h" 

// 引入几何库
#include "EDFM_Geometry_3D.h"
#include "14_DOP.h"            // 14-DOP 几何剔除
#include "MatrixEdge.h"        // 基岩棱定义

#include <iostream>
#include <algorithm>
#include <chrono>
#include <unordered_set>
#include <limits>
#include <fstream>
#include <iomanip>
#include <numeric>

// 常量定义
static constexpr double MIN_INTERSECTION_AREA = 1e-12; // 交互面积阈值，小于此值忽略
static constexpr double POINT_MERGE_TOLERANCE = 1e-7;  // 点合并容差


//构造函数
MeshManager_3D::MeshManager_3D(double lx, double ly, double lz,
    int nx, int ny, int nz,
    bool usePrism, bool useQuadBase)
    : lx_(lx), ly_(ly), lz_(lz),
    nx_(nx), ny_(ny), nz_(nz),
    usePrism_(usePrism), useQuadBase_(useQuadBase)
{
    cout << "调用MeshManager构造函数并初始化" << "\n";
}

//【基岩】 构造3D基岩网格
void MeshManager_3D::BuildSolidMatrixGrid_3D(NormalVectorCorrectionMethod corr, std::string matrix_mesh_Name)
{
    // 1) 构建基础网格拓扑
    mesh_.GenerateMesh3D(lx_, ly_, lz_, nx_, ny_, nz_, usePrism_, useQuadBase_, matrix_mesh_Name);

    //2）为每个网格的网格面构建索引
    mesh_.buildFaceGlobalId2LocalIndexMap();

    // 3) 单元分类 (Identify Solid vs Boundary Cells)
    mesh_.ClassifySolidMatrixCells();

    // 4)  识别边界组基于 Tag 分类
    bcGroups_byTag_ = BoundaryFaceClassify_byTag::ClassifyBoundaryFacesByTag_3D(mesh_);

    // 5) 可选：统计所有边界面（用于自检）
    boundaryCount_ = 0;
    for (const auto& f : mesh_.getFaces()) if (f.isBoundary()) ++boundaryCount_;

    // 6) （可选）自检打印
    std::cout << "[BC] faces: x0="<<bcGroups_byTag_.x0.size()
        << " xL=" << bcGroups_byTag_.xL.size()
        << " y0=" << bcGroups_byTag_.y0.size()
        << " yL=" << bcGroups_byTag_.yL.size()
        << " z0=" << bcGroups_byTag_.z0.size()
        << " zL=" << bcGroups_byTag_.zL.size()
        << " | total boundary=" << boundaryCount_ << "\n";

    // 7) 计算几何信息 输入非正交向量修正方法
    mesh_.ComputeMatrixMeshFaceGeometricInfor(corr);

    // 8）计算所有基岩单元的 AABB，为后续求交做准备
    mesh_.computeAllCellAABBs();

    //  9) 构建光栅化背景网格索引
    // 这里使用 nx_, ny_, nz_ 作为划分分辨率，通常与网格密度一致效果最佳
    mesh_.buildCellBins(nx_, ny_, nz_);

}

// =========================================================
// [New] 全局索引初始化实现
// =========================================================
void MeshManager_3D::setupGlobalIndices()
{
    // 1. 获取基岩总网格数 (作为裂缝编号的起点)
    // 假设 mesh_.getCells() 返回的是所有有效基岩单元
    int nMatrix = static_cast<int>(mesh_.getCells().size());

    std::cout << "[MeshManager] Setting up Global Solver Indices..." << std::endl;
    std::cout << "              Matrix Cells (Offset): " << nMatrix << std::endl;

    // 2. 通知裂缝网络开始分配索引
    int totalDOF = frNet_2D_.distributeSolverIndices(nMatrix);

    std::cout << "              Total DOF (Matrix + Fracs): " << totalDOF << std::endl;
}

int MeshManager_3D::getTotalDOFCount() const
{
    int total = static_cast<int>(mesh_.getCells().size());
    for (const auto& frac : frNet_2D_.getFractures())
    {
        total += static_cast<int>(frac.fracCells.size());
    }
    return total;
}

// =========================================================
// DOF 编排器 (DOF Mapper) 实现
// =========================================================

void MeshManager_3D::setNumDOFs(int nDof)
{
    if (nDof < 1)
    {
        std::cerr << "[Error] MeshManager_3D::setNumDOFs - Number of DOFs must be >= 1. Force setting to 1." << std::endl;
        num_dofs_ = 1;
    }
    else
    {
        num_dofs_ = nDof;
        std::cout << "[System] 3D DOF Mapper configured: " << num_dofs_ << " DOFs per cell/element." << std::endl;
    }
}

int MeshManager_3D::getNumDOFs() const
{
    return num_dofs_;
}

int MeshManager_3D::getEquationIndex(int solverIndex, int dofOffset) const
{
    if (solverIndex < 0 || solverIndex >= getTotalDOFCount())
    {
        std::cerr << "[Error] MeshManager_3D::getEquationIndex - solverIndex ("
            << solverIndex << ") is out of valid bounds [0, "
            << getTotalDOFCount() - 1 << "]." << std::endl;
        return -1;
    }

    if (dofOffset < 0 || dofOffset >= num_dofs_)
    {
        std::cerr << "[Error] MeshManager_3D::getEquationIndex - dofOffset ("
            << dofOffset << ") is out of valid bounds [0, "
            << num_dofs_ - 1 << "]." << std::endl;
        return -1;
    }

    return solverIndex * num_dofs_ + dofOffset;
}

int MeshManager_3D::getTotalEquationDOFs() const
{
    return getTotalDOFCount() * num_dofs_;
}

// =========================================================
// [New] 交互数据清洗与去重 (Production-Ready)
// =========================================================
void MeshManager_3D::removeDuplicateInteractions()
{
    std::cout << "[MeshManager] Cleaning up Interaction Pairs (Ghost Removal & Smart Deduplication)..." << std::endl;

    if (interactionPairs_.empty()) return;

    // ---------------------------------------------------------
    // 1. 计算物理容差 (Physical Tolerance)
    // ---------------------------------------------------------
    // 计算单个基岩网格的尺寸
    double dx = (nx_ > 0) ? lx_ / nx_ : 0.0;
    double dy = (ny_ > 0) ? ly_ / ny_ : 0.0;
    double dz = (nz_ > 0) ? lz_ / nz_ : 0.0;

    // 计算网格体对角线的一半作为最大理论相交距离
    // 理论上，点到面的距离如果超过 (对角线/2)，则该点所在的网格不可能与面相交
    double gridHalfDiag = 0.5 * std::sqrt(dx * dx + dy * dy + dz * dz);

    // 引入安全系数 (Safety Factor 1.01) 避免浮点误差误删边界情况
    double distTolerance = gridHalfDiag * 1.01;

    std::cout << "              -> Grid Size: " << dx << " x " << dy << " x " << dz << std::endl;
    std::cout << "              -> Ghost Distance Tolerance: " << distTolerance << std::endl;

    // ---------------------------------------------------------
    // 2. 执行清洗循环
    // ---------------------------------------------------------
    std::vector<InteractionPair> cleanPairs;
    cleanPairs.reserve(interactionPairs_.size());

    // 用于去重的 Map: KeyString -> Index in cleanPairs
    std::unordered_map<std::string, size_t> uniqueMap;

    int ghostsCount = 0;
    int duplicatesCount = 0;

    for (const auto& pair : interactionPairs_)
    {
        // --- A. 健壮性检查 ---
        if (pair.matrixSolverIndex < 0 || pair.fracCellSolverIndex < 0) continue;

        // --- B. 几何过滤 (Ghost Removal) ---
        // 检查基岩中心到裂缝平面的距离
        if (pair.distMatrixToFracPlane > distTolerance) {
            ghostsCount++;
            continue; // 距离过远，判定为幽灵交互，直接丢弃
        }

        // --- C. 智能去重 (Smart Deduplication) ---
        // 构建唯一指纹: "MatID_FracID_Area"
        // 使用定点精度 (保留6位小数) 确保数值稳定性
        std::stringstream ss;
        ss << pair.matrixSolverIndex << "_"
            << pair.fracCellSolverIndex << "_"
            << std::fixed << std::setprecision(6) << pair.intersectionArea;

        std::string key = ss.str();

        if (uniqueMap.find(key) == uniqueMap.end())
        {
            // 新的唯一有效交互 -> 保留
            uniqueMap[key] = cleanPairs.size();
            cleanPairs.push_back(pair);
        }
        else
        {
            // ID相同 且 面积完全相同 -> 判定为算法重复计算 -> 丢弃
            duplicatesCount++;
        }
    }

    // ---------------------------------------------------------
    // 3. 覆盖原始数据
    // ---------------------------------------------------------
    interactionPairs_ = std::move(cleanPairs);

    std::cout << "              -> Removed " << ghostsCount << " ghost interactions (distance > tolerance)." << std::endl;
    std::cout << "              -> Removed " << duplicatesCount << " exact duplicates." << std::endl;
    std::cout << "              -> Final Pair Count: " << interactionPairs_.size() << std::endl;
}

// =========================================================
// [Fix] 解决共面/边界双重计算问题 (修复 .norm() 报错)
// =========================================================
void MeshManager_3D::resolveCoplanarInteractions()
{
    std::cout << "[MeshManager] Resolving Coplanar/Boundary Overlaps..." << std::endl;

    if (interactionPairs_.empty()) return;

    // 1. 构建临时索引：FracElemSolverIndex -> List of Pair Indices
    std::unordered_map<int, std::vector<size_t>> fracElemToPairsMap;
    for (size_t i = 0; i < interactionPairs_.size(); ++i) {
        if (interactionPairs_[i].fracCellSolverIndex >= 0) {
            fracElemToPairsMap[interactionPairs_[i].fracCellSolverIndex].push_back(i);
        }
    }

    std::vector<bool> toRemove(interactionPairs_.size(), false);
    int removedCount = 0;
    int scaledCount = 0;

    // 2. 遍历每一个裂缝单元
    for (auto& it : fracElemToPairsMap)
    {
        std::vector<size_t>& pIndices = it.second;
        if (pIndices.empty()) continue;

        // --- 核心逻辑：检测空间重叠 (Spatial Overlap Check) ---
        for (size_t i = 0; i < pIndices.size(); ++i) {
            if (toRemove[pIndices[i]]) continue;

            for (size_t j = i + 1; j < pIndices.size(); ++j) {
                if (toRemove[pIndices[j]]) continue;

                const InteractionPair& pA = interactionPairs_[pIndices[i]];
                const InteractionPair& pB = interactionPairs_[pIndices[j]];

                // 检查1: 面积是否极其相似 (差异 < 1%)
                if (std::abs(pA.intersectionArea - pB.intersectionArea) > 1e-5) continue;

                // 检查2: 几何中心距离是否极近 (重合)
                // [Fix] 手动计算距离，替代 .norm()
                double dx = pA.polygonCenter.m_x - pB.polygonCenter.m_x;
                double dy = pA.polygonCenter.m_y - pB.polygonCenter.m_y;
                double dz = pA.polygonCenter.m_z - pB.polygonCenter.m_z;
                double distCenters = std::sqrt(dx * dx + dy * dy + dz * dz);

                if (distCenters > 1e-4) continue;

                // [判定] 发现了双胞胎！执行仲裁
                // 计算 Matrix Center 到 Poly Center 的距离
                // 获取 Cell 对象 (确保 mesh_ 成员变量可访问)
                // 注意：这里需要从 GlobalID 转换到 Index，假设 GlobalID 1-based
                // 如果您有 getCellIndex 函数更好，这里直接遍历或用 map

                int idxA = -1, idxB = -1;
                // 尝试使用 mesh_ 的接口查找
                try {
                    // 假设 mesh_.getCellId2Index() 可用
                    // 如果不可用，请替换为 mesh_.getCellIndex(pA.matrixCellGlobalID)
                    // 这里使用更通用的方式：
                    const auto& cellMap = mesh_.getCellId2Index();
                    if (cellMap.find(pA.matrixCellGlobalID) != cellMap.end())
                        idxA = cellMap.at(pA.matrixCellGlobalID);

                    if (cellMap.find(pB.matrixCellGlobalID) != cellMap.end())
                        idxB = cellMap.at(pB.matrixCellGlobalID);
                }
                catch (...) { continue; }

                if (idxA == -1 || idxB == -1) continue;

                const auto& cellA = mesh_.getCells()[idxA];
                const auto& cellB = mesh_.getCells()[idxB];

                // [Fix] 手动计算 distA
                double dAx = cellA.center.m_x - pA.polygonCenter.m_x;
                double dAy = cellA.center.m_y - pA.polygonCenter.m_y;
                double dAz = cellA.center.m_z - pA.polygonCenter.m_z;
                double distA = std::sqrt(dAx * dAx + dAy * dAy + dAz * dAz);

                // [Fix] 手动计算 distB
                double dBx = cellB.center.m_x - pB.polygonCenter.m_x;
                double dBy = cellB.center.m_y - pB.polygonCenter.m_y;
                double dBz = cellB.center.m_z - pB.polygonCenter.m_z;
                double distB = std::sqrt(dBx * dBx + dBy * dBy + dBz * dBz);

                if (distA < distB) {
                    toRemove[pIndices[j]] = true; // B 较远，删 B
                }
                else {
                    toRemove[pIndices[i]] = true; // A 较远，删 A
                    break;
                }
                removedCount++;
            }
        }
    }

    // 3. 执行删除
    if (removedCount > 0) {
        std::vector<InteractionPair> cleanPairs;
        cleanPairs.reserve(interactionPairs_.size() - removedCount);
        for (size_t i = 0; i < interactionPairs_.size(); ++i) {
            if (!toRemove[i]) {
                cleanPairs.push_back(interactionPairs_[i]);
            }
        }
        interactionPairs_ = std::move(cleanPairs);
    }

    // 4. [Final Step] 面积归一化
    // 重建映射
    fracElemToPairsMap.clear();
    for (size_t i = 0; i < interactionPairs_.size(); ++i) {
        fracElemToPairsMap[interactionPairs_[i].fracCellSolverIndex].push_back(i);
    }

    for (auto& it : fracElemToPairsMap) {
        std::vector<size_t>& pIndices = it.second;
        if (pIndices.empty()) continue;

        int macroID = interactionPairs_[pIndices[0]].fracMacroID;
        int fracGID = interactionPairs_[pIndices[0]].fracElementGlobalID;

        // 查找真实面积
        const Fracture_2D* f = findFractureByID(fracture_network(), macroID);
        if (!f) continue;

        int fLocalIdx = f->getElemIndex(fracGID);
        if (fLocalIdx == -1) continue;

        double realElemArea = f->fracCells[fLocalIdx].area;
        double totalInterArea = 0.0;
        for (size_t idx : pIndices) totalInterArea += interactionPairs_[idx].intersectionArea;

        // 容差判定
        if (totalInterArea > realElemArea + 1e-4) {
            double scaleFactor = realElemArea / totalInterArea;
            for (size_t idx : pIndices) {
                interactionPairs_[idx].intersectionArea *= scaleFactor;
            }
            scaledCount++;
        }
    }

    std::cout << "              -> Resolved " << removedCount << " coplanar conflicts." << std::endl;
    std::cout << "              -> Normalized " << scaledCount << " elements." << std::endl;
}


// =========================================================
// [New] 拓扑映射构建实现
// =========================================================
void MeshManager_3D::buildTopologyMaps()
{
    std::cout << "[MeshManager] Building Topology Maps (Matrix <-> Pairs <-> Fracture)..." << std::endl;
    
    // 1. 重置并预分配 Matrix Map
    // Matrix Solver Index 范围是 [0, numCells - 1]
    int numMatrixCells = static_cast<int>(mesh_.getCells().size());
    mat2InteractionMap_.assign(numMatrixCells, std::vector<const InteractionPair*>());

    // 2. 重置 Frac Map
    frac2InteractionMap_.clear();

    // 3. 遍历扁平化列表，分发指针
    for (const auto& pair : interactionPairs_)
    {
        // 映射到 Matrix
        int mIdx = pair.matrixSolverIndex; 
        if (mIdx >= 0 && mIdx < numMatrixCells) {
            mat2InteractionMap_[mIdx].push_back(&pair);
        }
        else {
            // 理论上不应发生，除非 SolverIndex 越界
            std::cerr << "[Warning] InteractionPair has invalid Matrix SolverIndex: " << mIdx << std::endl;
        }

        // 映射到 Fracture
        int fIdx = pair.fracCellSolverIndex;
        if (fIdx != -1) {
            frac2InteractionMap_[fIdx].push_back(&pair);
        }
    }

    std::cout << "              -> Indexed " << interactionPairs_.size() << " pairs." << std::endl;
    std::cout << "              -> Matrix Cells with Interactions: " << std::count_if(mat2InteractionMap_.begin(), mat2InteractionMap_.end(), [](const auto& v) { return !v.empty(); }) << std::endl;
    std::cout << "              -> Fracture Elements with Interactions: " << frac2InteractionMap_.size() << std::endl;
}

const std::vector<const InteractionPair*>& MeshManager_3D::getInteractionsOfMatrix(int matrixSolverIdx) const
{
    if (matrixSolverIdx >= 0 && matrixSolverIdx < static_cast<int>(mat2InteractionMap_.size())) {
        return mat2InteractionMap_[matrixSolverIdx];
    }
    return emptyPairList_;
}

const std::vector<const InteractionPair*>& MeshManager_3D::getInteractionsOfFracture(int fracSolverIdx) const
{
    auto it = frac2InteractionMap_.find(fracSolverIdx);
    if (it != frac2InteractionMap_.end()) {
        return it->second;
    }
    return emptyPairList_;
}

const Fracture_2D* MeshManager_3D::findFractureByID(const FractureNetwork_2D& net, int fracID) {
    for (const auto& f : net.getFractures()) {
        if (f.id == fracID) return &f;
    }
    return nullptr;
}

//【基岩】 动态构建基岩单元的棱 (Edge)
void MeshManager_3D::_buildLocalMatrixEdges(const Cell& cell, std::vector<MatrixEdge>& outEdges) const
{
    outEdges.clear();
    // 使用 HashSet 去重
    std::unordered_set<MatrixEdge, MatrixEdgeHash> edgeSet;

    // 遍历单元的所有面
    for (int faceID : cell.CellFaceIDs)
    {
        int fIdx = mesh_.getFaceIndex(faceID);
        if (fIdx == -1) continue;
        const Face& face = mesh_.getFaces()[fIdx];

        const auto& nodes = face.FaceNodeIDs;
        size_t n = nodes.size();
        if (n < 2) continue;

        // 遍历面的边 (Node i -> Node i+1)
        for (size_t i = 0; i < n; ++i)
        {
            int n1 = nodes[i];
            int n2 = nodes[(i + 1) % n]; // 闭合
            edgeSet.insert(MatrixEdge(n1, n2)); // MatrixEdge 构造函数会自动排序 n1, n2
        }
    }
    // 转存到 vector
    outEdges.reserve(edgeSet.size());
    for (const auto& edge : edgeSet)
    {
        outEdges.push_back(edge);
    }
}

//【裂缝】 添加2D裂缝至裂缝网络
void MeshManager_3D::addFracturetoFractureNetwork(const Fracture_2D& frac)
{
    frNet_2D_.addFracture(frac);
}

//【裂缝】 对裂缝网络内部所有2D宏观裂缝进行网格独立划分并重新构建全局裂缝网格段索引，其中自行调用了 frNet_2D_.rebuildGlobalIndex()函数
void MeshManager_3D::meshAllFracturesinNetwork(int nU, int nV, NormalVectorCorrectionMethod method)
{
    frNet_2D_.meshAllFractures(nU, nV, method);
}

//【裂缝】 检测宏观裂缝网络内的裂缝-裂缝相交 (Element-wise),区分
void MeshManager_3D::DetectFractureFractureIntersectionsInNetwork(FFIntersectionStrategy strategy)
{
    frNet_2D_.DetectFractureFractureIntersections(strategy);
}

// 【后处理】
///【基岩】 输出基岩网格信息至Txt文件
void MeshManager_3D::exportMeshInfortoTxt(const std::string& prefix) const
{
    mesh_.exportToTxt(prefix);
}

///【裂缝】 输出裂缝网格信息至Txt文件
void MeshManager_3D::exportFracturesNetworkInfortoTxt(const std::string& prefix)const
{
    frNet_2D_.exportNetworkToTxt(prefix);
}

void MeshManager_3D::exportFracturesNetworkInfortoTxt_improved(const std::string& prefix)const
{
    frNet_2D_.exportNetworkToTxt_improved(prefix);
}

///【裂缝】导出特定宏观裂缝裂缝微观裂缝的网格面非正交信息至.csv
void MeshManager_3D::inspectFractureEdges_non_orthogonalInfor(int fracID, const std::string& prefix) const
{
    frNet_2D_.inspectFractureEdges(fracID, prefix);
}

///【裂缝】检查并导出 F-F 交线及微观网格对应关系至 CSV
void MeshManager_3D::inspectIntersections_FracturetoFracture(const std::string& prefix) const
{
    frNet_2D_.inspectIntersections(prefix);
}

// =========================================================
// 辅助函数: 多边形重构 (排序 & 面积计算) - [高精度曲面修正版]
// =========================================================
bool MeshManager_3D::_reconstructPolygon(const std::vector<Vector>& rawPoints, const Vector& referenceNormal, InteractionPair& outPair) const
{
    if (rawPoints.size() < 3) return false;

    // 1. 去重 (O(N^2) for small N)
    std::vector<Vector> uniquePoints;
    uniquePoints.reserve(rawPoints.size());
    for (const auto& p : rawPoints)
    {
        bool duplicate = false;
        for (const auto& exist : uniquePoints)
        {
            if ((p - exist).Mag() < POINT_MERGE_TOLERANCE)
            {
                duplicate = true;
                break;
            }
        }
        if (!duplicate) uniquePoints.push_back(p);
    }

    if (uniquePoints.size() < 3) return false;

    // 2. 计算几何中心 (Centroid)
    Vector center(0, 0, 0);
    for (const auto& p : uniquePoints) center = center + p;
    center = center / (double)uniquePoints.size();

    // 3. 排序 (Sort Counter-Clockwise)
    // 虽然是曲面，但局部来看投影到参考法向平面的拓扑序是正确的
    // 这一步仅用于确定点的连接顺序，不涉及面积计算，所以投影是安全的
    Vector axisZ = referenceNormal;
    Vector tempX = (std::abs(axisZ.m_x) < 0.9) ? Vector(1, 0, 0) : Vector(0, 1, 0);
    Vector axisY = (axisZ & tempX); axisY = axisY / axisY.Mag();
    Vector axisX = (axisY & axisZ); axisX = axisX / axisX.Mag();

    std::sort(uniquePoints.begin(), uniquePoints.end(), [&](const Vector& a, const Vector& b) {
        Vector va = a - center;
        Vector vb = b - center;
        double angA = std::atan2(va * axisY, va * axisX);
        double angB = std::atan2(vb * axisY, vb * axisX);
        return angA < angB;
        });

    // 4. [CRITICAL FIX] 计算 3D 曲面面积
    // 不再使用 (Cross * Normal) 的投影公式
    // 而是直接累加每个扇形三角形的模长 (Area = 0.5 * |Cross|)
    // 这能精确捕捉偏离平均平面的微小起伏面积
    double area3D = 0.0;
    Vector w_normal(0, 0, 0); // 加权法向

    size_t n = uniquePoints.size();
    for (size_t i = 0; i < n; ++i)
    {
        const Vector& p1 = uniquePoints[i];
        const Vector& p2 = uniquePoints[(i + 1) % n];

        // 构建微元三角形 (Center -> P1 -> P2)
        Vector v1 = p1 - center;
        Vector v2 = p2 - center;
        Vector cross = v1 & v2;

        double triArea = 0.5 * cross.Mag(); // 真实的 3D 面积
        area3D += triArea;

        // 累加法向以计算更准确的平均法向
        if (cross * referenceNormal < 0) cross = cross * -1.0; // 统一方向
        w_normal = w_normal + cross;
    }

    // 5. 阈值判断
    if (area3D < MIN_INTERSECTION_AREA) return false;

    // 6. 填充结果
    outPair.polygonPoints = uniquePoints;
    outPair.polygonCenter = center;

    // 更新为更精确的面积加权法向
    double w_mag = w_normal.Mag();
    outPair.polygonNormal = (w_mag > 1e-12) ? (w_normal / w_mag) : referenceNormal;

    outPair.intersectionArea = area3D; // [FIX] 使用 3D 面积

    return true;
}

// =========================================================
// Improved 3D-EDFM 核心求交实现 (Triangulation + Safe Access)
// =========================================================
void MeshManager_3D::SolveIntersection3D_improved(IntersectionStrategy strategy)
{
    std::cout << "=========================================================" << std::endl;
    std::cout << "[3D-EDFM] Starting SolveIntersection3D_improved..." << std::endl;
    std::cout << "          Strategy: " << (strategy == IntersectionStrategy::Octree_Optimized ? "Octree Optimized" : "Brute Force") << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    // 1. 清理旧数据
    interactionPairs_.clear();
    int totalInteractions = 0;

    // 2. 准备八叉树 (如果选用 Octree 策略)
    FaceIndexedOctree* matrixOctree = nullptr;
    if (strategy == IntersectionStrategy::Octree_Optimized)
    {
        std::cout << "  -> Strategy: Octree Optimized" << std::endl;
        std::cout << "  -> Calculating Global AABB for Matrix Mesh..." << std::endl;

        // 计算基岩网格的全局 AABB
        Vector minP(1e30, 1e30, 1e30);
        Vector maxP(-1e30, -1e30, -1e30);

        const auto& nodesMap = mesh_.getNodesMap();
        if (nodesMap.empty()) {
            std::cerr << "[Error] Matrix Mesh is empty! Aborting Intersection." << std::endl;
            return;
        }

        for (const auto& pair : nodesMap) {
            const Vector& p = pair.second.coord;
            if (p.m_x < minP.m_x) minP.m_x = p.m_x;
            if (p.m_y < minP.m_y) minP.m_y = p.m_y;
            if (p.m_z < minP.m_z) minP.m_z = p.m_z;

            if (p.m_x > maxP.m_x) maxP.m_x = p.m_x;
            if (p.m_y > maxP.m_y) maxP.m_y = p.m_y;
            if (p.m_z > maxP.m_z) maxP.m_z = p.m_z;
        }
        AABB globalBox(minP - Vector(0.1, 0.1, 0.1), maxP + Vector(0.1, 0.1, 0.1));

        std::cout << "  -> Building FaceIndexedOctree..." << std::endl;
        const auto& faces = mesh_.getFaces();
        // 构造八叉树 (存储的是 vector index)
        matrixOctree = new FaceIndexedOctree(globalBox, faces);
    }
    else {
        std::cout << "  -> Strategy: Brute Force (Slow)" << std::endl;
    }

    // 3. 遍历所有宏观裂缝
    // 假设 FractureNetwork_2D 提供了 getFractures()
    const auto& fractures = frNet_2D_.getFractures();

    std::cout << "  -> Processing " << fractures.size() << " Macro Fractures..." << std::endl;

    for (const auto& frac : fractures)
    {
        // 遍历该宏观裂缝下的所有微观单元 (Level 2)
        for (const auto& fracElem : frac.fracCells)
        {
            // ---------------------------------------------
            // Step A: Broad Phase (粗筛 - 确定候选基岩单元)
            // ---------------------------------------------
            std::unordered_set<int> candidateCellIDs; // 使用 Set 去重

            if (strategy == IntersectionStrategy::BruteForce)
            {
                const auto& cells = mesh_.getCells();
                for (const auto& cell : cells) {
                    candidateCellIDs.insert(cell.id);
                }
            }
            else if (matrixOctree)
            {
                // Octree 模式：查询裂缝微元 AABB
                std::vector<int> hitFaceIndices; // 注意：这里返回的是 Index
                matrixOctree->query(fracElem.boundingBox, hitFaceIndices);

                // 从 Face 反查 Cell (Owner & Neighbor)
                const auto& faces = mesh_.getFaces();
                for (int fid : hitFaceIndices) {
                    // 安全性检查
                    if (fid >= 0 && fid < (int)faces.size()) {
                        const Face& f = faces[fid]; // 直接使用 Index 访问
                        if (f.ownerCell != -1) candidateCellIDs.insert(f.ownerCell);
                        if (f.neighborCell != -1) candidateCellIDs.insert(f.neighborCell);
                    }
                }
            }

            // 准备 Frac 顶点的缓存
            std::vector<Vector> fracPolyCoords;
            for (int nodeIdx : fracElem.nodeIndices) {
                fracPolyCoords.push_back(frac.fracNodes[nodeIdx].coord);
            }

            // [核心改进] 三角化裂缝微元 (Triangulation Strategy)
            // 将四边形拆分为两个三角形，以适应非共面 Twisted 情况
            std::vector<std::vector<Vector>> fracTriangles;
            if (fracPolyCoords.size() == 4) {
                fracTriangles.push_back({ fracPolyCoords[0], fracPolyCoords[1], fracPolyCoords[2] });
                fracTriangles.push_back({ fracPolyCoords[0], fracPolyCoords[2], fracPolyCoords[3] });
            }
            else if (fracPolyCoords.size() == 3) {
                fracTriangles.push_back(fracPolyCoords);
            }

            // ---------------------------------------------
            // Step B: Narrow Phase (精细求交 - Type 1/2/3)
            // ---------------------------------------------
            for (int cellID : candidateCellIDs)
            {
                // [Safe Access] 获取 Cell Index
                int cellIdx = mesh_.getCellIndex(cellID);
                if (cellIdx == -1) continue;
                const Cell& matrixCell = mesh_.getCells()[cellIdx];

                // 二级 AABB 检测
                if (!matrixCell.boundingBox.overlaps(fracElem.boundingBox)) continue;

                std::vector<Vector> rawIntersectionPoints;

                // =========================================================
                // Type 0: Check for Coplanar Faces (处理完美共面)
                // =========================================================
                bool isCoplanarCase = false;
                for (int faceID : matrixCell.CellFaceIDs)
                {
                    // [Safe Access] Face ID -> Face Index
                    int fIdx = mesh_.getFaceIndex(faceID);
                    if (fIdx == -1) continue;
                    const Face& mFace = mesh_.getFaces()[fIdx];

                    // 1. 法向平行检查
                    double dot = mFace.normal * fracElem.normal;
                    if (std::abs(std::abs(dot) - 1.0) < 1e-3)
                    {
                        // 2. 共面检查 (点到平面距离 ~ 0)
                        Vector v = fracPolyCoords[0] - mFace.FaceNodeCoords[0];
                        if (std::abs(v * mFace.normal) < 1e-3)
                        {
                            // 确认共面！执行 2D 多边形重叠计算
                            if (EDFM_Geometry_3D::GetCoplanarPolygonIntersection(
                                mFace.FaceNodeCoords, fracPolyCoords, mFace.normal, rawIntersectionPoints))
                            {
                                isCoplanarCase = true;
                            }
                        }
                    }
                }

                // === 如果非共面，执行基于三角化的 Type 1/2/3 检测 ===
                if (!isCoplanarCase)
                {
                    // === Type 1: 裂缝节点在基岩单元内部 ===
                    for (const auto& pt : fracPolyCoords) {
                        if (EDFM_Geometry_3D::IsPointInConvexCell(pt, matrixCell, mesh_)) {
                            rawIntersectionPoints.push_back(pt);
                        }
                    }

                    // === Type 2: 裂缝边 (Segment) vs 基岩面 (Triangulated) ===
                    size_t numFN = fracPolyCoords.size();
                    for (size_t i = 0; i < numFN; ++i)
                    {
                        Vector pStart = fracPolyCoords[i];
                        Vector pEnd = fracPolyCoords[(i + 1) % numFN];

                        for (int faceID : matrixCell.CellFaceIDs)
                        {
                            // [Safe Access]
                            int fIdx = mesh_.getFaceIndex(faceID);
                            if (fIdx == -1) continue;
                            const Face& mFace = mesh_.getFaces()[fIdx];

                            // 为了鲁棒性，也将基岩面三角化 (虽然基岩面通常是平面的)
                            std::vector<Vector> mFacePts = mFace.FaceNodeCoords;
                            if (mFacePts.size() >= 3) {
                                // 简单的 Fan Triangulation (0-i-i+1)
                                for (size_t k = 1; k < mFacePts.size() - 1; ++k) {
                                    Vector triA = mFacePts[0];
                                    Vector triB = mFacePts[k];
                                    Vector triC = mFacePts[k + 1];

                                    Vector hit;
                                    // 使用线段-三角形求交
                                    if (EDFM_Geometry_3D::IntersectSegmentTriangle(
                                        pStart, pEnd, triA, triB, triC, mFace.normal, hit))
                                    {
                                        rawIntersectionPoints.push_back(hit);
                                    }
                                }
                            }
                        }
                    }

                    // === Type 3: 基岩棱 (Segment) vs 裂缝面 (Triangulated) ===
                    std::vector<MatrixEdge> matrixEdges;
                    // 注意：_buildLocalMatrixEdges 内部也必须使用 safe index
                    _buildLocalMatrixEdges(matrixCell, matrixEdges);

                    const auto& nodesMap = mesh_.getNodesMap();
                    for (const auto& mEdge : matrixEdges)
                    {
                        auto it1 = nodesMap.find(mEdge.n1);
                        auto it2 = nodesMap.find(mEdge.n2);
                        if (it1 == nodesMap.end() || it2 == nodesMap.end()) continue;

                        Vector mStart = it1->second.coord;
                        Vector mEnd = it2->second.coord;

                        // [核心] 只要 Matrix Edge 穿过任意一个 Frac Triangle，就算相交
                        for (const auto& tri : fracTriangles)
                        {
                            Vector hit;
                            // 动态计算该三角形的真实法向
                            Vector triNormal = (tri[1] - tri[0]) & (tri[2] - tri[0]);
                            double area = triNormal.Mag();
                            if (area > 1e-12) triNormal = triNormal / area;
                            else triNormal = fracElem.normal; // 退化回平均法向

                            if (EDFM_Geometry_3D::IntersectSegmentTriangle(
                                mStart, mEnd, tri[0], tri[1], tri[2], triNormal, hit))
                            {
                                rawIntersectionPoints.push_back(hit);
                            }
                        }
                    }
                } // End if(!isCoplanarCase)

                // ---------------------------------------------
                // Step C: Polygon Reconstruction
                // ---------------------------------------------
                if (rawIntersectionPoints.size() >= 3)
                {
                    InteractionPair pair(matrixCell.id, fracElem.id, frac.id);

                    // 2. [新增] 填充 Solver Index (Critical Step!)
                    // ---------------------------------------------------------
                    // Matrix Index: 
                    // 基岩的 SolverIndex 通常就是其在 vector 中的下标 (0 ~ Nm-1)
                    pair.matrixSolverIndex = cellIdx;

                    // Fracture Index:
                    // 直接从 fracElem 中读取 (前提是 distributeSolverIndices 已调用)
                    pair.fracCellSolverIndex = fracElem.solverIndex;

                    // [安全检查] 防止未初始化索引
                    if (pair.fracCellSolverIndex == -1) {
                        std::cerr << "[Error] Fracture Element " << fracElem.id
                            << " has invalid SolverIndex! Did you call setupGlobalIndices()?" << std::endl;
                        // 可以选择 continue 跳过或 throw 异常
                    }

                    // 注意：对于 Twisted 裂缝，投影法向使用 fracElem.normal 是可接受的近似
                    if (_reconstructPolygon(rawIntersectionPoints, fracElem.normal, pair))
                    {
                        // 计算 d_NNC
                        // 公式: d = | (Center_Matrix - Center_Poly) * Normal_Poly |
                        Vector vecDist = matrixCell.center - pair.polygonCenter;

                        // 点乘计算投影距离 (取绝对值)
                        pair.distMatrixToFracPlane = std::abs(vecDist * pair.polygonNormal);

                        // [安全性保障] 防止数值噪音导致的极小负值（虽然 abs 已保证非负，但可做极小值钳制）
                        // 这里的钳制主要为了后续做分母时的数值稳定，但在几何计算阶段如实记录即可
                        if (pair.distMatrixToFracPlane < 1e-15) pair.distMatrixToFracPlane = 1e-15;

                        interactionPairs_.push_back(pair);
                        totalInteractions++;
                    }
                }

            } // End Loop Candidate Cells
        } // End Loop Micro Elements
    } // End Loop Macro Fractures

    // 4. 清理资源
    if (matrixOctree) {
        delete matrixOctree;
        matrixOctree = nullptr;
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "[3D-EDFM] Intersection Analysis Completed." << std::endl;
    std::cout << "          Total Interaction Pairs Found: " << totalInteractions << std::endl;
    std::cout << "          Time Elapsed: " << elapsed.count() << " s." << std::endl;
    std::cout << "=========================================================" << std::endl;
}

// =========================================================
// [核心算法] 3D-EDFM 核心求交接口 (Twist Acceleration + Rasterization)
// =========================================================
void MeshManager_3D::SolveIntersection3D_improved_twist_accleration(IntersectionStrategy strategy)
{
    // 1. 初始化日志
    std::cout << "\n[EDFM-3D] Start Twist-Accelerated Intersection..." << std::endl;
    std::cout << "          Strategy: ";
    if (strategy == IntersectionStrategy::BruteForce) std::cout << "BruteForce";
    else if (strategy == IntersectionStrategy::Octree_Optimized) std::cout << "Octree (AABB)";
    else if (strategy == IntersectionStrategy::Octree_With_14DOP) std::cout << "Octree + 14-DOP (High Perf)";
    else if (strategy == IntersectionStrategy::Rasterization_14DOP) std::cout << "Rasterization + 14-DOP (Ultra Perf)";
    std::cout << std::endl;

    auto startTime = std::chrono::high_resolution_clock::now();

    interactionPairs_.clear();

    // 获取引用
    const auto& fractures = frNet_2D_.getFractures();
    const auto& matrixCells = mesh_.getCells();
    const auto& matrixFaces = mesh_.getFaces();
    const auto& matrixNodesMap = mesh_.getNodesMap();

    // =========================================================
    // 2. 空间索引准备 (Octree vs Rasterization)
    // =========================================================
    FaceIndexedOctree* localOctree = nullptr;

    // 判断是否使用 Octree
    bool useOctree = (strategy == IntersectionStrategy::Octree_Optimized ||
        strategy == IntersectionStrategy::Octree_With_14DOP);

    // 判断是否使用光栅化
    bool useRaster = (strategy == IntersectionStrategy::Rasterization_14DOP);

    // [Safety Check] 如果策略是 Rasterization 但背景网格未构建，回退到 BruteForce 或尝试构建
    if (useRaster && mesh_.getGridParams().nx == 0) {
        std::cout << "  [Warning] Background grid empty. Attempting to build..." << std::endl;
        // 尝试使用默认密度构建 (假设 nx_, ny_, nz_ 可用，或者使用硬编码)
        // mesh_.buildCellBins(50, 50, 20); // 示例
        // 如果无法构建，回退
        if (mesh_.getGridParams().nx == 0) {
            std::cout << "  [Fallback] Rasterization failed. Fallback to BruteForce." << std::endl;
            useRaster = false;
        }
    }

    if (useOctree) {
        Vector minP(1e30, 1e30, 1e30), maxP(-1e30, -1e30, -1e30);
        bool hasNode = false;

        for (const auto& pair : matrixNodesMap) {
            const Vector& p = pair.second.coord;
            minP.m_x = std::min(minP.m_x, p.m_x); minP.m_y = std::min(minP.m_y, p.m_y); minP.m_z = std::min(minP.m_z, p.m_z);
            maxP.m_x = std::max(maxP.m_x, p.m_x); maxP.m_y = std::max(maxP.m_y, p.m_y); maxP.m_z = std::max(maxP.m_z, p.m_z);
            hasNode = true;
        }

        if (hasNode) {
            AABB globalBox(minP - Vector(0.1, 0.1, 0.1), maxP + Vector(0.1, 0.1, 0.1));
            std::cout << "  -> Building Local Face-Indexed Octree..." << std::endl;
            localOctree = new FaceIndexedOctree(globalBox, matrixFaces);
        }
        else {
            useOctree = false;
        }
    }

    long long totalSubTriangles = 0;
    long long totalCandidateChecks = 0;
    long long totalExactChecks = 0;

    // 3. 遍历所有宏观裂缝
    for (const auto& frac : fractures)
    {
        const auto& elements = frac.fracCells;
        const auto& fracNodes = frac.fracNodes;

        // 4. 遍历裂缝微元
        for (const auto& elem : elements)
        {
            // 4.1 获取微元顶点
            std::vector<Vector> elemPoints;
            elemPoints.reserve(elem.nodeIndices.size());
            for (int nid : elem.nodeIndices) {
                if (nid >= 0 && nid < static_cast<int>(fracNodes.size())) {
                    elemPoints.push_back(fracNodes[nid].coord);
                }
            }
            if (elemPoints.size() < 3) continue;

            // 4.2 划分子三角形 (Per-Sub-Triangle)
            std::vector<std::vector<Vector>> subTriangles;
            if (elemPoints.size() == 3) subTriangles.push_back(elemPoints);
            else if (elemPoints.size() == 4) {
                subTriangles.push_back({ elemPoints[0], elemPoints[1], elemPoints[2] });
                subTriangles.push_back({ elemPoints[0], elemPoints[2], elemPoints[3] });
            }
            else {
                for (size_t i = 1; i < elemPoints.size() - 1; ++i) {
                    subTriangles.push_back({ elemPoints[0], elemPoints[i], elemPoints[i + 1] });
                }
            }

            // 5. 遍历子三角形
            for (const auto& triPts : subTriangles)
            {
                totalSubTriangles++;

                // 计算法向 (先叉乘再归一化)
                Vector triNormal = (triPts[1] - triPts[0]) & (triPts[2] - triPts[0]);
                triNormal.Normalize();

                // 5.1 构建 AABB (用于 Broad Phase)
                AABB triAABB;
                triAABB.min = Vector(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
                triAABB.max = Vector(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest());
                for (const auto& p : triPts) {
                    triAABB.min.m_x = std::min(triAABB.min.m_x, p.m_x); triAABB.min.m_y = std::min(triAABB.min.m_y, p.m_y); triAABB.min.m_z = std::min(triAABB.min.m_z, p.m_z);
                    triAABB.max.m_x = std::max(triAABB.max.m_x, p.m_x); triAABB.max.m_y = std::max(triAABB.max.m_y, p.m_y); triAABB.max.m_z = std::max(triAABB.max.m_z, p.m_z);
                }

                //// 解决完美平面裂缝 (Zero-Thickness) 导致的边界丢失或数值剔除问题
                //// 尤其是当裂缝恰好落在 Grid 边界 (x=50) 时，这一步至关重要！
                double eps = 1e-5;
                triAABB.min = triAABB.min - Vector(eps, eps, eps);
                triAABB.max = triAABB.max + Vector(eps, eps, eps);

                // 5.2 构建 14-DOP (仅当策略需要时)
                Box14DOP tri14DOP;
                if (strategy == IntersectionStrategy::Octree_With_14DOP || strategy == IntersectionStrategy::Rasterization_14DOP) {
                    tri14DOP.fromTriangle(triPts[0], triPts[1], triPts[2]);
                }

                // 5.3 获取候选基岩单元 (Broad Phase)
                std::vector<int> candidateCells;

                if (useRaster)
                {
                    // [Task 3] 光栅化查询
                    mesh_.getCandidateCells_Rasterization(triPts[0], triPts[1], triPts[2], candidateCells);
                }
                else if (useOctree && localOctree)
                {
                    // [Octree] 查询
                    std::vector<int> candidateFaces;
                    localOctree->query(triAABB, candidateFaces);

                    std::unordered_set<int> uniqueCellIDs;
                    for (int fid : candidateFaces) {
                        if (fid < 0 || fid >= static_cast<int>(matrixFaces.size())) continue;
                        const Face& f = matrixFaces[fid];
                        if (f.ownerCell >= 0) {
                            int idx = mesh_.getCellIndex(f.ownerCell);
                            if (idx != -1) uniqueCellIDs.insert(idx);
                        }
                        if (f.neighborCell >= 0) {
                            int idx = mesh_.getCellIndex(f.neighborCell);
                            if (idx != -1) uniqueCellIDs.insert(idx);
                        }
                    }
                    candidateCells.assign(uniqueCellIDs.begin(), uniqueCellIDs.end());
                }
                else
                {
                    // [Fallback] 暴力全量
                    candidateCells.resize(matrixCells.size());
                    std::iota(candidateCells.begin(), candidateCells.end(), 0);
                }

                // 6. 遍历候选基岩单元 (Narrow Phase)
                for (int cellIndex : candidateCells)
                {
                    totalCandidateChecks++;
                    const Cell& mCell = matrixCells[cellIndex];

                    // --- 14-DOP 剔除 ---
                    if (strategy == IntersectionStrategy::Octree_With_14DOP || strategy == IntersectionStrategy::Rasterization_14DOP) {
                        if (!tri14DOP.overlaps(mCell.boundingBox)) continue;
                    }

                    totalExactChecks++;
                    std::vector<Vector> rawInterPoints;

                    // ========================================================
                    // Type 0: Coplanar (共面处理)
                    // ========================================================
                    bool isCoplanar = false;
                    for (int faceID : mCell.CellFaceIDs) {
                        int fIdx = mesh_.getFaceIndex(faceID);
                        if (fIdx == -1) continue;
                        const Face& mFace = matrixFaces[fIdx];

                        // 检查法向平行 (|n1 . n2| ≈ 1)
                        if (std::abs(std::abs(mFace.normal * triNormal) - 1.0) < 1e-3) {
                            // 检查点面距离 (d ≈ 0)
                            if (std::abs((triPts[0] - mFace.FaceNodeCoords[0]) * mFace.normal) < 1e-3) {
                                // 执行共面裁剪
                                if (EDFM_Geometry_3D::GetCoplanarPolygonIntersection(
                                    mFace.FaceNodeCoords, triPts, mFace.normal, rawInterPoints)) {
                                    isCoplanar = true;
                                    // 注意：共面情况只需找到一个面重叠即可 break，
                                    // 或者是把所有共面交集都算上？通常只需一个有效重叠
                                    // 这里我们假设如果共面，则该 Interaction 仅由共面构成
                                    break;
                                }
                            }
                        }
                    }

                    // 如果不是共面，则进行常规的 3D 穿透检测
                    if (!isCoplanar)
                    {
                        // ========================================================
                        // Type 1: Triangle Vertex inside Matrix
                        // ========================================================
                        for (const auto& tp : triPts) {
                            if (EDFM_Geometry_3D::IsPointInConvexCell(tp, mCell, mesh_)) {
                                rawInterPoints.push_back(tp);
                            }
                        }

                        // ========================================================
                        // Type 2: Triangle Edge vs Matrix Face (Triangulated)
                        // ========================================================
                        std::vector<std::pair<Vector, Vector>> triEdges = {
                            {triPts[0], triPts[1]}, {triPts[1], triPts[2]}, {triPts[2], triPts[0]}
                        };

                        for (int faceID : mCell.CellFaceIDs) {
                            int fIdx = mesh_.getFaceIndex(faceID);
                            if (fIdx == -1) continue;
                            const Face& mFace = matrixFaces[fIdx];

                            const auto& mFacePts = mFace.FaceNodeCoords;
                            if (mFacePts.size() >= 3) {
                                for (size_t k = 1; k < mFacePts.size() - 1; ++k) {
                                    for (const auto& seg : triEdges) {
                                        Vector outPt;
                                        if (EDFM_Geometry_3D::IntersectSegmentTriangle(
                                            seg.first, seg.second,
                                            mFacePts[0], mFacePts[k], mFacePts[k + 1],
                                            mFace.normal, outPt))
                                        {
                                            rawInterPoints.push_back(outPt);
                                        }
                                    }
                                }
                            }
                        }

                        // ========================================================
                        // Type 3: Matrix Edge vs Triangle Face
                        // ========================================================
                        std::vector<MatrixEdge> cellEdges;
                        _buildLocalMatrixEdges(mCell, cellEdges);

                        for (const auto& edge : cellEdges) {
                            auto it1 = matrixNodesMap.find(edge.n1);
                            auto it2 = matrixNodesMap.find(edge.n2);
                            if (it1 == matrixNodesMap.end() || it2 == matrixNodesMap.end()) continue;

                            const Vector& p1 = it1->second.coord;
                            const Vector& p2 = it2->second.coord;
                            Vector outPt;
                            if (EDFM_Geometry_3D::IntersectSegmentTriangle(p1, p2,
                                triPts[0], triPts[1], triPts[2], triNormal, outPt)) {
                                rawInterPoints.push_back(outPt);
                            }
                        }
                    } // End if(!isCoplanar)

                    // 7. 重构交互多边形
                    if (rawInterPoints.size() >= 3) {
                        InteractionPair pair(mCell.id, elem.id, frac.id);

                        // [新增] 填充索引 ------------------------------
                        int matrixIdx = mesh_.getCellIndex(mCell.id); // 再次确认索引
                        pair.matrixSolverIndex = matrixIdx;
                        pair.fracCellSolverIndex = elem.solverIndex;

                        if (pair.fracCellSolverIndex == -1) {
                            std::cerr << "[Error] TwistAlgo: FracElem " << elem.id << " has no SolverIndex!" << std::endl;
                        }
                        // ---------------------------------------------
                        
                        // 使用子三角形法向 triNormal 辅助重构
                        if (_reconstructPolygon(rawInterPoints, triNormal, pair)) {

                            // [新增 Step 1] 计算 d_NNC
                            // 此时 pair.polygonNormal 已经是基于交互多边形精确计算的加权法向
                            // pair.polygonCenter 也是多边形的精确质心
                            Vector vecDist = mCell.center - pair.polygonCenter;
                            
                            pair.distMatrixToFracPlane = std::abs(vecDist * pair.polygonNormal);

                            // [安全性保障] 
                            // 对于基岩中心恰好落在裂缝面上的特殊情况（虽然概率极低），给予极小值防止除零
                            if (pair.distMatrixToFracPlane < 1e-15) pair.distMatrixToFracPlane = 1e-15;
                            interactionPairs_.push_back(pair);
                        }
                    }
                }
            }
        }
    }

    if (localOctree) {
        delete localOctree;
        localOctree = nullptr;
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> ms = endTime - startTime;

    std::cout << "  [Completed] Time: " << ms.count() << " ms." << std::endl;
    std::cout << "  - Sub-Triangles Processed : " << totalSubTriangles << std::endl;
    std::cout << "  - AABB Candidate Checks   : " << totalCandidateChecks << std::endl;
    std::cout << "  - Exact Geometric Checks  : " << totalExactChecks << std::endl;

    // 统计剔除率
    if (totalCandidateChecks > 0) {
        double reduction = 100.0 * (1.0 - (double)totalExactChecks / totalCandidateChecks);
        std::cout << "  - 14-DOP Culling Rate     : " << std::fixed << std::setprecision(2) << reduction << "%" << std::endl;
    }
    std::cout << "  - Interaction Pairs Found : " << interactionPairs_.size() << std::endl;
    std::cout << "=========================================================\n" << std::endl;
}


// =========================================================
// 导出 CSV 验证
// =========================================================
void MeshManager_3D::exportInteractionPairsToCSV(const std::string& filename) const
{
    std::ofstream file(filename);
    if (!file.is_open()) return;

    file << "MatrixID,FracElemID,FracMacroID,Area,Dist_NNC,CenterX,CenterY,CenterZ,NormalX,NormalY,NormalZ,NumPoints\n";
    file << std::fixed << std::setprecision(8);

    for (const auto& pair : interactionPairs_)
    {
        file << pair.matrixCellGlobalID << ","
            << pair.fracElementGlobalID << ","
            << pair.fracMacroID << ","
            << pair.intersectionArea << ","
            << pair.distMatrixToFracPlane << "," 
            << pair.polygonCenter.m_x << "," << pair.polygonCenter.m_y << "," << pair.polygonCenter.m_z << ","
            << pair.polygonNormal.m_x << "," << pair.polygonNormal.m_y << "," << pair.polygonNormal.m_z << ","
            << pair.polygonPoints.size() << "\n";
    }
    file.close();
    std::cout << "[Export] Interaction Pairs exported to " << filename << std::endl;

}

void MeshManager_3D::exportInteractionPolygonsToTxt(const std::string& filename) const
{
    std::ofstream ofs(filename);
    if (!ofs.is_open())
    {
        std::cerr << "[Error] Failed to open file for export: " << filename << std::endl;
        return;
    }

    // 设置高精度，防止多边形顶点重合导致的绘制闪烁
    ofs << std::fixed << std::setprecision(8);

    int count = 0;
    for (const auto& pair : interactionPairs_)
    {
        const auto& pts = pair.polygonPoints;
        if (pts.empty()) continue;

        // 格式: 顶点数 x1 y1 z1 x2 y2 z2 ...
        ofs << pts.size();
        for (const auto& p : pts)
        {
            ofs << " " << p.m_x << " " << p.m_y << " " << p.m_z;
        }
        ofs << "\n";
        count++;
    }

    ofs.close();
    std::cout << "[Export] " << count << " Interaction Polygons exported to " << filename << std::endl;
}

void MeshManager_3D::exportInteractionPolygonsToTxt_improved(const std::string& filename) const
{
    std::ofstream ofs(filename+"_InteractionPolygons.txt");
    if (!ofs.is_open())
    {
        std::cerr << "[Error] Failed to open file for export: " << filename << std::endl;
        return;
    }

    // 设置高精度
    ofs << std::fixed << std::setprecision(8);

    // [可选] 输出表头注释 (MATLAB 读取时需跳过或使用 comment style)
    // ofs << "% MatrixGID MatrixSID FracMacroID FracElemGID FracElemSID Area Dist NumPoints x1 y1 z1 ...\n";

    int count = 0;
    for (const auto& pair : interactionPairs_)
    {
        const auto& pts = pair.polygonPoints;
        if (pts.empty()) continue;

        // 1. 输出拓扑与物理属性
        ofs << pair.matrixCellGlobalID << " "
            << pair.matrixSolverIndex << " "
            << pair.fracMacroID << " "
            << pair.fracElementGlobalID << " "
            << pair.fracCellSolverIndex << " "
            << pair.intersectionArea << " "
            << pair.distMatrixToFracPlane << " ";

        // 2. 输出几何顶点
        ofs << pts.size();
        for (const auto& p : pts)
        {
            ofs << " " << p.m_x << " " << p.m_y << " " << p.m_z;
        }
        ofs << "\n";
        count++;
    }

    ofs.close();
    std::cout << "[Export] " << count << " Interaction Polygons exported to " << filename << std::endl;
}

// =========================================================
// [Debug Visualization] 导出搜索空间验证数据
// =========================================================

void MeshManager_3D::exportSearchSpaceToTxt(const std::string& filename, int targetFracID, IntersectionStrategy strategy)
{
    std::ofstream file(filename);
    if (!file.is_open()) return;

    // Status 1 (Yellow): Broad Phase 选中但被 14-DOP 剔除
    // Status 2 (Red)   : 最终进入 Exact Check
    file << std::fixed << std::setprecision(6);

    const auto& fractures = frNet_2D_.getFractures();
    const auto& matrixCells = mesh_.getCells();
    const auto& matrixFaces = mesh_.getFaces();

    // 1. 找到目标裂缝
    const Fracture_2D* pFrac = nullptr;
    for (const auto& f : fractures) {
        if (f.id == targetFracID) { pFrac = &f; break; }
    }
    if (!pFrac) return;

    std::cout << "[Debug] Exporting search space (" << filename << ")..." << std::endl;

    // 2. 准备 Octree (如果策略需要)
    FaceIndexedOctree* debugOctree = nullptr;
    bool useOctree = (strategy == IntersectionStrategy::Octree_Optimized ||
        strategy == IntersectionStrategy::Octree_With_14DOP);

    if (useOctree) {
        // ... (原有的 Octree 构建代码保持不变) ...
        Vector minP(1e30, 1e30, 1e30), maxP(-1e30, -1e30, -1e30);
        const auto& nodesMap = mesh_.getNodesMap();
        for (const auto& pair : nodesMap) {
            const Vector& p = pair.second.coord;
            minP.m_x = std::min(minP.m_x, p.m_x); minP.m_y = std::min(minP.m_y, p.m_y); minP.m_z = std::min(minP.m_z, p.m_z);
            maxP.m_x = std::max(maxP.m_x, p.m_x); maxP.m_y = std::max(maxP.m_y, p.m_y); maxP.m_z = std::max(maxP.m_z, p.m_z);
        }
        AABB globalBox(minP - Vector(0.1, 0.1, 0.1), maxP + Vector(0.1, 0.1, 0.1));
        debugOctree = new FaceIndexedOctree(globalBox, matrixFaces);
    }

    std::unordered_map<int, int> cellStatusMap;

    // 3. 遍历裂缝微元
    for (const auto& elem : pFrac->fracCells)
    {
        // ... (构建子三角形代码保持不变) ...
        std::vector<Vector> elemPoints;
        for (int nid : elem.nodeIndices) elemPoints.push_back(pFrac->fracNodes[nid].coord);
        std::vector<std::vector<Vector>> subTriangles;
        if (elemPoints.size() == 4) {
            subTriangles.push_back({ elemPoints[0], elemPoints[1], elemPoints[2] });
            subTriangles.push_back({ elemPoints[0], elemPoints[2], elemPoints[3] });
        }
        else { subTriangles.push_back(elemPoints); }

        for (const auto& triPts : subTriangles)
        {
            Box14DOP tri14DOP;
            tri14DOP.fromTriangle(triPts[0], triPts[1], triPts[2]);

            // ... (AABB 构建代码保持不变) ...
            AABB triAABB;
            triAABB.min = Vector(1e30, 1e30, 1e30); triAABB.max = Vector(-1e30, -1e30, -1e30);
            for (const auto& p : triPts) {
                triAABB.min.m_x = std::min(triAABB.min.m_x, p.m_x); triAABB.min.m_y = std::min(triAABB.min.m_y, p.m_y); triAABB.min.m_z = std::min(triAABB.min.m_z, p.m_z);
                triAABB.max.m_x = std::max(triAABB.max.m_x, p.m_x); triAABB.max.m_y = std::max(triAABB.max.m_y, p.m_y); triAABB.max.m_z = std::max(triAABB.max.m_z, p.m_z);
            }

            std::vector<int> candidates;

            // [新增] BruteForce 分支
            if (strategy == IntersectionStrategy::BruteForce) {
                // 全量添加
                candidates.resize(matrixCells.size());
                std::iota(candidates.begin(), candidates.end(), 0);
            }
            else if (strategy == IntersectionStrategy::Rasterization_14DOP) {
                mesh_.getCandidateCells_Rasterization(triPts[0], triPts[1], triPts[2], candidates);
            }
            else if (useOctree && debugOctree) {
                // ... (Octree 查询代码保持不变) ...
                std::vector<int> faces;
                debugOctree->query(triAABB, faces);
                std::unordered_set<int> cSet;
                for (int fid : faces) {
                    if (fid < 0 || fid >= matrixFaces.size()) continue;
                    const Face& f = matrixFaces[fid];
                    if (f.ownerCell >= 0) cSet.insert(mesh_.getCellIndex(f.ownerCell));
                    if (f.neighborCell >= 0) cSet.insert(mesh_.getCellIndex(f.neighborCell));
                }
                candidates.assign(cSet.begin(), cSet.end());
            }

            // 状态标记
            for (int cellID : candidates)
            {
                const Cell& cell = matrixCells[cellID];
                int currentStatus = 1;

                bool pass = false;

                // [Fix] BruteForce 和 OctreePure 都不进行 14-DOP 剔除
                if (strategy == IntersectionStrategy::BruteForce ||
                    strategy == IntersectionStrategy::Octree_Optimized) {
                    pass = true;
                }
                else {
                    if (tri14DOP.overlaps(cell.boundingBox)) {
                        pass = true;
                    }
                }

                if (pass) currentStatus = 2; // Red

                if (cellStatusMap.find(cellID) == cellStatusMap.end()) cellStatusMap[cellID] = currentStatus;
                else cellStatusMap[cellID] = std::max(cellStatusMap[cellID], currentStatus);
            }
        }
    }
    if (debugOctree) delete debugOctree;

    // 写入文件
    for (const auto& pair : cellStatusMap) {
        int cellID = pair.first;
        int status = pair.second;
        const AABB& box = matrixCells[cellID].boundingBox;
        file << cellID << " " << status << " "
            << box.min.m_x << " " << box.min.m_y << " " << box.min.m_z << " "
            << box.max.m_x << " " << box.max.m_y << " " << box.max.m_z << "\n";
    }
    file.close();
    std::cout << "  -> Exported to " << filename << std::endl;
}