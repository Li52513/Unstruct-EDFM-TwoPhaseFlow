#include "MeshManager.h"
#include "Mesh.h"
#include "FractureNetwork.h"
#include "FractureCommon.h"
#include "BoundaryFaceClassify.h"
#include "BoundaryFaceClassify_byTag.h"
#include <iostream>

//构造函数
MeshManager::MeshManager(double lx, double ly, double lz,
    int nx, int ny, int nz,
    bool usePrism, bool useQuadBase)
    : lx_(lx), ly_(ly), lz_(lz),
    nx_(nx), ny_(ny), nz_(nz),
    usePrism_(usePrism), useQuadBase_(useQuadBase)
{
    cout<<"调用MeshManager构造函数并初始化"<< "\n";
}

// 构造2D基岩网格
void MeshManager::BuildSolidMatrixGrid_2D(NormalVectorCorrectionMethod corr)
{
    // 1) 构建基础网格 (调用新的 2D 函数)
    mesh_.GenerateMesh2D(lx_, ly_, nx_, ny_, useQuadBase_,"2D_Mesh_test");

	mesh_.buildFaceGlobalId2LocalIndexMap(); // 构建 Face 索引映射

    // 2) 单元分类 (Identify Solid vs Boundary Cells)
    mesh_.ClassifySolidMatrixCells();

    // 3)  识别边界组：改为基于 Tag 分类
    bcGroups_byTag_ = BoundaryFaceClassify_byTag::ClassifyBoundaryFacesByTag_2D(mesh_);

    // 4) 统计与自检 (逻辑不变，但增加了 Custom 的输出)
    boundaryCount_ = 0;
    for (const auto& f : mesh_.getFaces()) {
        if (f.physicalGroupId != -1) ++boundaryCount_;
    }

    std::cout << "[BC] Groups based on Tags:\n"
        << "  Left  (Tag " << MeshTags::LEFT << "): " << bcGroups_byTag_.x0.size() << "\n"
        << "  Right (Tag " << MeshTags::RIGHT << "): " << bcGroups_byTag_.xL.size() << "\n"
        << "  Bottom(Tag " << MeshTags::BOTTOM << "): " << bcGroups_byTag_.y0.size() << "\n"
        << "  Top   (Tag " << MeshTags::TOP << "): " << bcGroups_byTag_.yL.size() << "\n";

    // 打印自定义边界（如障碍物）
    //for (auto const& [tag, list] : bcGroups_.custom) {
    //    std::cout << "  Custom(Tag " << tag << "): " << list.size() << "\n";
    //}
    std::cout << "  Total Boundary Faces: " << boundaryCount_ << "\n";

    // 6) 计算几何信息
    mesh_.ComputeMatrixMeshFaceGeometricInfor(corr);
}


void MeshManager::addFracture(const Vector& s, const Vector& e) 
{
    frNet_.addFracture(s, e);
}

void MeshManager::DetectAndSubdivideFractures(IntersectionSearchStrategy_2D searchStrategy)
{
    frNet_.invalidateFracElemIndex();

    // 1) 先找裂缝–裂缝的交点
    frNet_.DetectFracturetoFractureIntersections();

    int totalFaces = mesh_.getFaces().size();
    int totalCandidates = 0;            //用于记录精确检测的候选面
    int totalaabbCheckCount = 0;        //用于记录AABB碰撞检测的面数
    int totaltrueIntersectionCount = 0; //用于记录精确求交计算的面数
    IntersectionStatistics_2D stats;

    // 2) 然后找裂缝–网格面的交点  可以选取AABB加速算法
    for (auto& F : frNet_.fractures)
    {
        F.DetectFracturetoMeshFaceIntersections_improved(mesh_, mesh_.getCells(), mesh_.getNodesMap(), searchStrategy, stats);
        totalCandidates += F.candidateFaceCount_;
        totalaabbCheckCount += stats.aabbCheckCount;
        totaltrueIntersectionCount += stats.geometryCheckCount;

    }
    std::cout << "[加速效果] 总面数 = " << totalFaces
        << ", 总裂缝数 = " << frNet_.fractures.size()
        << ", 总候选面数 = " << totalCandidates
        << "进行AABB包围盒碰撞测试的面数 = " << totalaabbCheckCount
        << "进行精准相交测试的网格面数 = " << totaltrueIntersectionCount
        << ", 平均每条裂缝面数 = " << (double)totalCandidates / frNet_.fractures.size()
        << ", 剔除率 = " << (100.0 - 100.0 * totalCandidates / (totalFaces * frNet_.fractures.size()))
        << "%\n";

    // 3) 给全局裂缝–裂缝交点去重并编号
    frNet_.DeduplicateAndRenumberFractureToFractureIntersections();

    // 4) 把全局裂缝–裂缝交点也插入到每条裂缝的 intersections 中
    frNet_.DistributeFracture_FractureIntersectionsToGlobalInersections();

    // 5) 重新排序 & 编号，然后划分
    for (auto& F : frNet_.fractures)
    {
        F.sortAndRenumberIntersections();
        F.subdivide(mesh_.getCells(), mesh_.getNodesMap(), distanceMetric_);
    }
    frNet_.rebuildFracElemIndex(); // 重建裂缝段索引

    frNet_.printFractureInfo();
}


void MeshManager::ComputeFractureGeometryCouplingCoefficient()
{
    for (auto& F : frNet_.fractures)
        F.computeGeometryCouplingCoefficientgeomCIandgeomAlpha();
}

void MeshManager::setDFNRandomSeed(unsigned seed)
{
	frNet_.setRandomSeed(seed);
}

void MeshManager::generateDFN(int N,
    const Vector& minPoint,
    const Vector& maxPoint,
    double Lmin,
    double Lmax,
    double alpha,
    double kappa,
    bool avoidOverlap)
{
	frNet_.generateDFN(N, minPoint, maxPoint, Lmin, Lmax, alpha, kappa, avoidOverlap);
}

void MeshManager::exportMesh(const std::string& pref) const {
    mesh_.exportToTxt(pref);
}

void MeshManager::exportFractures(const std::string& pref) const {
    frNet_.exportToTxt(pref);
}

void MeshManager::printFractureInfo() const {
    frNet_.printFractureInfo();
}


// =========================================================
// 构建全局系统索引实现
// =========================================================
void MeshManager::BuildGlobalSystemIndexing()
{
    std::cout << "=========================================================" << std::endl;
    std::cout << "[System] Building Global Indexing for AX=B..." << std::endl;

    // 1. 锁定基岩自由度
    // 直接使用网格单元数量，默认 Row i 对应 cells_[i]
    n_matrix_dofs_ = mesh_.getGridCount();

    if (n_matrix_dofs_ == 0) {
        std::cerr << "[Warning] Matrix grid is empty! Indexing might be invalid." << std::endl;
    }

    std::cout << "       -> Matrix DOFs (0 ~ " << n_matrix_dofs_ - 1 << "): "
        << n_matrix_dofs_ << std::endl;

    // 2. 初始化全局查找表
    // 前 n_matrix_dofs_ 个位置对应基岩，填 nullptr
    globalSolverIndexToElementPtr_.clear();
    // 预估一个容量，避免频繁 realloc，虽然后面是 push_back
    // 假设裂缝单元数约为基岩数的 10% (仅预估)
    globalSolverIndexToElementPtr_.reserve(static_cast<size_t>(n_matrix_dofs_ * 1.1));

    // 填充基岩部分的占位符
    globalSolverIndexToElementPtr_.resize(n_matrix_dofs_, nullptr);

    // 3. 分配裂缝自由度
    // 裂缝编号从 n_matrix_dofs_ 开始连续递增
    int current_dof = n_matrix_dofs_;
    int frac_elem_count = 0;

    //遍历所有裂缝对象
    for (size_t fid = 0; fid < frNet_.fractures.size(); ++fid)
    {
        auto& frac = frNet_.fractures[fid];

        for (auto& elem : frac.elements)
        {
            // 赋予全局唯一解算器索引
            elem.solverIndex = current_dof++;

            // 赋予所属宏观裂缝 ID (0-based)
            elem.parentFractureID = static_cast<int>(fid);

            // [新增] 将当前单元的地址注册到全局查找表
            // 注意：必须取 elem 的地址 (&elem)
            globalSolverIndexToElementPtr_.push_back(&elem);

            frac_elem_count++;
        }
    }
    std::cout << "       -> Fracture DOFs (" << n_matrix_dofs_ << " ~ " << current_dof - 1 << "): "
        << frac_elem_count << std::endl;
    std::cout << "       -> Total System Size: " << current_dof << std::endl;
    std::cout << "       -> Lookup Table Size: " << globalSolverIndexToElementPtr_.size() << std::endl;
    
    // 4. 在索引分配完毕后，构建裂缝单元-基岩拓扑映射
    mesh_.rebuildFractureMap(frNet_.fractures);

    frNet_.buildSolverIndexCache(n_matrix_dofs_);

    std::cout << "=========================================================" << std::endl;
}

// =========================================================
// 查找接口实现
// =========================================================
const FractureElement* MeshManager::getFractureElementBySolverIndex(int solverIdx) const
{
    // 边界检查
    if (solverIdx < 0 || solverIdx >= static_cast<int>(globalSolverIndexToElementPtr_.size())) {
        return nullptr;
    }
    // 直接查表 (O(1))
    return globalSolverIndexToElementPtr_[solverIdx];
}


void MeshManager::BuildFracturetoFractureTopology()
{
    std::cout << "[Topology] Building Fracture-Fracture Connectivity..." << std::endl;
    
    // 1. 清空旧拓扑
    for (auto& frac : frNet_.fractures)
    {
        for (auto& elem : frac.elements)
        {
            elem.neighbors.clear();
        }
    }

    // 2. 构建 "交点ID -> 连接的裂缝单元列表" 的临时映射
    // Key: GlobalFFPoint ID
    // Value: List of elements connecting to this point
    struct ElementRef {
        // [寻址用] 用于在 vector 中找到对象
        int fracIdx;      // 宏观裂缝 vector 下标 (0-based)
        int elemVecIdx;   // 裂缝单元 vector 下标 (0-based)

        // [数据用] 用于写入 FractureNeighbor
        int solverIndex;  // 全局求解器索引
        int elemLocalID;  // 单元局部 ID (1-based, elem.id)
    };
    std::map<int, std::vector<ElementRef>> intersectionMap;

    // 遍历所有裂缝单元，注册它们连接的交点
    for (size_t fIdx = 0; fIdx < frNet_.fractures.size(); ++fIdx)
    {
        auto& frac = frNet_.fractures[fIdx];
        for (size_t eIdx = 0; eIdx < frac.elements.size(); ++eIdx)
        {
            const auto& elem = frac.elements[eIdx];

            // 准备引用数据
            ElementRef ref;
			ref.fracIdx = static_cast<int>(fIdx); //也等于elem.parentFractureID
            ref.elemVecIdx = static_cast<int>(eIdx);
            ref.solverIndex = elem.solverIndex;
            ref.elemLocalID = elem.id;

            // 检查起点
            if (elem.isFFatStart && elem.gIDstart != -1) {
                intersectionMap[elem.gIDstart].push_back(ref);
            }
            // 检查终点
            if (elem.isFFatEnd && elem.gIDend != -1) {
                intersectionMap[elem.gIDend].push_back(ref);
            }
        }
    }
    // 3. 建立连接：全连接逻辑 (Clique)
    int totalConnections = 0;
    for (auto const& kv : intersectionMap)
    {
        int ffPointID = kv.first;
        const auto& connectedElems = kv.second;

        if (connectedElems.size() < 2) continue; // 孤立点，无连接

        // 双重循环：让每个单元认识同一个交点上的其他所有单元
        for (size_t i = 0; i < connectedElems.size(); ++i)
        {
            for (size_t j = 0; j < connectedElems.size(); ++j)
            {
                if (i == j) continue; // 跳过自己

                const auto& selfRef = connectedElems[i];     // 我
                const auto& neighborRef = connectedElems[j]; // 邻居

                // [寻址] 找到“我”这个对象，以便写入 neighbors
                auto& selfElem = frNet_.fractures[selfRef.fracIdx].elements[selfRef.elemVecIdx];

                // [数据] 将邻居的信息写入我的 neighbors 列表
                // 注意：这里使用的是 neighborRef 的数据字段
                FractureElement::FractureNeighbor neighborInfo(
                    neighborRef.fracIdx,      // 邻居属于哪条大裂缝
                    neighborRef.elemLocalID,  // 邻居叫什么名字 (1-based ID)
                    neighborRef.solverIndex,  // 邻居在矩阵里住哪行 (Solver Index)
                    ffPointID                 // 我们在哪见的面 (Intersection ID)
                );

                selfElem.neighbors.push_back(neighborInfo);
                totalConnections++;
            }
        }
    }

    std::cout << "           -> Processed " << intersectionMap.size() << " intersection points." << std::endl;
    std::cout << "           -> Established " << totalConnections << " fracture-fracture links." << std::endl;
    std::cout << "=========================================================" << std::endl;

}

int MeshManager::getTotalDOFCount() const
{
    // Total = Matrix DOFs + Fracture DOFs
    // 由于我们在 BuildGlobalSystemIndexing 中是连续分配的，
    // 我们可以直接计算，或者如果要精确也可以遍历 frNet (但效率低)。
    // 这里我们动态计算一次作为校验，或者你可以增加一个成员变量 n_total_dofs_ 记录。
    // 为了不引入新变量，这里使用简单的计算：
    int total = n_matrix_dofs_;
    for (const auto& frac : frNet_.fractures) {
        total += static_cast<int>(frac.elements.size());
    }
    return total;
}

// =========================================================
// DOF 编排器 (DOF Mapper) 实现
// =========================================================
void MeshManager::setNumDOFs(int nDof)
{
    if (nDof < 1)
    {
        std::cerr << "[Error] MeshManager::setNumDOFs - Number of DOFs must be >= 1. Force setting to 1." << std::endl;
        num_dofs_ = 1;
    }
    else
    {
        num_dofs_ = nDof;
        std::cout << "[System] DOF Mapper configured: " << num_dofs_ << " DOFs per cell/element." << std::endl;
    }
}

int MeshManager::getNumDOFs() const
{
    return num_dofs_;
}

int MeshManager::getEquationIndex(int solverIndex, int dofOffset) const
{
    // 安全性检查：方程组装时高频调用，边界判定能有效拦截越界段错误
    if (solverIndex < 0 || solverIndex >= getTotalDOFCount())
    {
        std::cerr << "[Error] MeshManager::getEquationIndex - solverIndex ("
            << solverIndex << ") is out of valid bounds [0, "
            << getTotalDOFCount() - 1 << "]." << std::endl;
        return -1;
    }

    if (dofOffset < 0 || dofOffset >= num_dofs_)
    {
        std::cerr << "[Error] MeshManager::getEquationIndex - dofOffset ("
            << dofOffset << ") is out of valid bounds [0, "
            << num_dofs_ - 1 << "]." << std::endl;
        return -1;
    }

    // 核心映射机制：交织存储 (Interleaved Storage)
    return solverIndex * num_dofs_ + dofOffset;
}

int MeshManager::getTotalEquationDOFs() const
{
    return getTotalDOFCount() * num_dofs_;
}
