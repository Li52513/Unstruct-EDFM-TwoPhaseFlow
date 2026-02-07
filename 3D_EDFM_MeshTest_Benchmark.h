#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <map>
#include <unordered_map>
#include <algorithm>

#include "3D_MeshManager.h"

// =========================================================
// 内部辅助工具
// =========================================================
namespace {
    // 字符串格式化辅助
    template <typename T>
    std::string vec_to_csv_string1(const std::vector<T>& v) {
        if (v.empty()) return "[]";
        std::stringstream ss;
        ss << "[";
        for (size_t i = 0; i < v.size(); ++i) {
            ss << v[i] << (i == v.size() - 1 ? "" : ";");
        }
        ss << "]";
        return ss.str();
    }
}

// =========================================================
// 核心诊断导出函数实现
// =========================================================
void export3DTopologyDiagnosticCSV(MeshManager_3D& mgr, const std::string& caseName)
{
    std::string sep = ",";
    std::cout << "\n[Benchmark] >>> Starting 3D-EDFM Topology Diagnostic Export..." << std::endl;

    // ---------------------------------------------------------
        // 1. 导出 3D 基岩单元拓扑 (Matrix Cell Topology)
        // ---------------------------------------------------------
    std::string fn_matrix = caseName + "_Diagnostic_MatrixCells.csv";
    std::ofstream out_m(fn_matrix);
    if (!out_m.is_open()) {
        std::cerr << "[Error] Unable to create file: " << fn_matrix << std::endl;
        return;
    }

    // 表头：包含基础信息、邻居信息以及交互多边形信息
    out_m << "GlobalID" << sep
        << "LocalIndex" << sep
        << "SolverIndex" << sep
        << "Center_X" << sep << "Center_Y" << sep << "Center_Z" << sep
        << "Volume" << sep
        << "NumFaces" << sep
        << "FaceIDs_Locals" << sep
        << "FaceIDs_Globlas" << sep
        << "Neighbor_LocalIndices" << sep // 邻居 SolverIndex
        << "Neighbor_GlobalIDs" << sep    // 邻居 Gmsh Tag
        << "Interactions(PolyArea:Frac[GID:LID:SID])" // [New] 交互溯源列
        << "\n";

    const auto& cells = mgr.mesh().getCells();
    const auto& faces = mgr.mesh().getFaces();
    const auto& frNet = mgr.fracture_network(); // 用于辅助查询

    for (size_t i = 0; i < cells.size(); ++i)
    {
        const auto& cell = cells[i];
        int mySolverIndex = -1;
        try {
			mySolverIndex = mgr.mesh().getCellId2Index().at(cell.id);
        }
        catch (...) {
            std::cerr << "[Error] Cell " << cell.id << " not found in cellId2Index!" << std::endl;
        }

        // --- A. 收集网格面及邻居信息 ---
        std::vector<int> faceIDs_Locals;
        std::vector<int> faceIDs_Globals;
        std::vector<int> neigh_Locals;
        std::vector<int> neigh_Globals;

        for (int fid : cell.CellFaceIDs)
        {
            int fIdx = mgr.mesh().getFaceIndex(fid);
            if (fIdx == -1 || fIdx >= static_cast<int>(faces.size())) {
                std::cerr << "[Error] Face ID " << fid << " in Cell " << cell.id
                    << " maps to invalid index " << fIdx << std::endl;
                continue; // 跳过坏数据
            }
            if (fIdx != -1)
            {
                const auto& face = faces[fIdx];
                faceIDs_Locals.push_back(fIdx);
                faceIDs_Globals.push_back(face.id);

                int myIdx = mySolverIndex;
                int neighborIdx = (face.ownerCell_index == myIdx) ? face.neighborCell_index : face.ownerCell_index;

                neigh_Locals.push_back(neighborIdx); // SolverIndex / LocalIndex

                if (neighborIdx >= 0 && neighborIdx < static_cast<int>(cells.size())) {
                    neigh_Globals.push_back(cells[neighborIdx].id); // GlobalID
                }
                else {
                    neigh_Globals.push_back(-1);
                }
            }
        }

        // B. [Optimized] 收集交互多边形信息 (直接调用 Manager 接口)
         std::vector<std::string> interactionInfos;

         // 调用 mgr.getInteractionsOfMatrix 获取关联的 Pairs
         const auto& myPairs = mgr.getInteractionsOfMatrix(mySolverIndex);

         for (const auto* pair : myPairs) 
         {
                std::stringstream ss;
                // 格式: Area:FracID(GID:LID:SID)
                ss << std::fixed << std::setprecision(6) << pair->intersectionArea << ":Frac" << pair->fracMacroID << "(";
                
                // 获取裂缝单元的详细索引 (GID, LID, SID)
                // GID 和 SID 已在 pair 中，LID 需要通过 Fracture 对象查询
                int elemLID = -1;
                const Fracture_2D* f = mgr.findFractureByID(frNet, pair->fracMacroID);
                if (f) {
                    elemLID = f->getElemIndex(pair->fracElementGlobalID);
                }
                ss << pair->fracElementGlobalID << ":" << elemLID << ":" << pair->fracCellSolverIndex << ")";
                interactionInfos.push_back(ss.str());
         }
            out_m << cell.id << sep
                << mySolverIndex << sep
                << mySolverIndex << sep
                << cell.center.m_x << sep << cell.center.m_y << sep << cell.center.m_z << sep
                << cell.volume << sep
                << cell.CellFaceIDs.size() << sep
                << vec_to_csv_string1(faceIDs_Locals) << sep
                << vec_to_csv_string1(faceIDs_Globals) << sep
                << vec_to_csv_string1(neigh_Locals) << sep
                << vec_to_csv_string1(neigh_Globals) << sep
                << vec_to_csv_string1(interactionInfos)
                << "\n";
        
    }
    out_m.close();
    std::cout << "   -> Exported: " << fn_matrix << std::endl;

    // ---------------------------------------------------------
    // 2. 导出 2D 裂缝微元拓扑 (Fracture Micro Elements)
    // ---------------------------------------------------------
    std::string fn_frac = caseName + "_Diagnostic_FractureElements.csv";
    std::ofstream out_f(fn_frac);

    out_f << "MacroFracID" << sep
        << "Elem_GlobalID" << sep
        << "Elem_LocalIndex" << sep
        << "Elem_SolverIndex" << sep
        << "ParentID_Check" << sep
        << "Area" << sep
        << "Normal_X" << sep << "Normal_Y" << sep << "Normal_Z" << sep
        << "Interactions(PolyArea:FractureCell[GID:LID:SID])" // [New] 交互溯源列
        << "\n";

    const auto& fractures = mgr.fracture_network().getFractures();
    for (const auto& frac : fractures)
    {
        for (const auto& elem : frac.fracCells)
        {
            std::vector<std::string> interactionInfos;

            if (elem.solverIndex != -1)
            {
                // 调用 mgr.getInteractionsOfFracture
                const auto& myPairs = mgr.getInteractionsOfFracture(elem.solverIndex);

                for (const auto* pair : myPairs)
                {
                    std::stringstream ss;
                    // 格式: Area:Mat(GID:LID:SID)
                    ss << std::fixed << std::setprecision(6) << pair->intersectionArea << ":Mat(";

                    // Matrix 的 LID 通常等于 SID
                    int matLID = pair->matrixSolverIndex;

                    ss << pair->matrixCellGlobalID << ":" << matLID << ":" << pair->matrixSolverIndex << ")";
                    interactionInfos.push_back(ss.str());
                }
            }

            out_f << frac.id << sep
                << elem.id << sep
                << elem.localIndex << sep
                << elem.solverIndex << sep
                << elem.parentFractureID << sep
                << elem.area << sep
                << elem.normal.m_x << sep << elem.normal.m_y << sep << elem.normal.m_z << sep
                << vec_to_csv_string(interactionInfos)
                << "\n";
        }
    }
    out_f.close();
    std::cout << "   -> Exported: " << fn_frac << std::endl;

    // ---------------------------------------------------------
    // 3. 导出 2D 裂缝边拓扑 (Fracture Edges / Micro Faces)
    // ---------------------------------------------------------
    std::string fn_fedge = caseName + "_Diagnostic_FractureEdges.csv";
    std::ofstream out_fe(fn_fedge);

    out_fe << "Edge_GlobalID_0based" << sep
        << "MacroFracID" << sep
        << "Owner_GlobalID" << sep << "Owner_LocalIdx" << sep << "Owner_SolverIndex" << sep
        << "Neighbor_GlobalID" << sep << "Neighbor_LocalIdx" << sep << "Neighbor_SolverIndex" << sep
        << "Length" << "\n";

    // 优先使用 FVM 求解器实际使用的全局边列表
    const auto& globalEdges = mgr.fracture_network().getGlobalEdges();

    if (globalEdges.empty()) {
        std::cout << "   [Warning] Global Edges empty! Did you call rebuildEdgeProperties? Skipping Edge Export." << std::endl;
    }
    else {
        for (const auto& edge : globalEdges) {
            out_fe << edge.id << sep
                << edge.parentFractureID << sep
                << edge.ownerCellID << sep << edge.ownerCellLocalIndex << sep << edge.ownerCell_solverIndex << sep
                << edge.neighborCellID << sep << edge.neighborCellLocalIndex << sep << edge.neighborCell_solverIndex << sep
                << edge.length << "\n";
        }
    }
    out_fe.close();
    std::cout << "   -> Exported: " << fn_fedge << " (Total: " << globalEdges.size() << ")" << std::endl;

    // ---------------------------------------------------------
    // 4. 导出基岩-裂缝交互多边形列表 (M-F NNC Flat List)
    // ---------------------------------------------------------
    std::string fn_mf = caseName + "_Diagnostic_MF_Interactions.csv";
    std::ofstream out_mf_nnc(fn_mf);

    out_mf_nnc << "Matrix_GlobalID" << sep << "Matrix_SolverIndex" << sep
        << "Frac_MacroID" << sep
        << "Frac_Elem_GlobalID" << sep << "Frac_Elem_SolverIndex" << sep
        << "Interaction_Area" << sep
        << "Dist_MatrixToPlane" << sep
        << "Poly_Center_X" << sep << "Poly_Center_Y" << sep << "Poly_Center_Z" <<sep
        << "Normal_X" << sep << "Normal_Y" << sep << "Normal_Z" << sep
        << "\n";

    const auto& pairs = mgr.getInteractionPairs();
    for (const auto& pair : pairs)
    {
        out_mf_nnc << pair.matrixCellGlobalID << sep << pair.matrixSolverIndex << sep
            << pair.fracMacroID << sep
            << pair.fracElementGlobalID << sep << pair.fracCellSolverIndex << sep
            << pair.intersectionArea
            << pair.distMatrixToFracPlane << sep
            << pair.polygonCenter.m_x << sep << pair.polygonCenter.m_y << sep << pair.polygonCenter.m_z << sep
            << pair.polygonNormal.m_x << sep << pair.polygonNormal.m_y << sep << pair.polygonNormal.m_z << sep
            << "\n";
    }
    out_mf_nnc.close();
    std::cout << "   -> Exported: " << fn_mf << " (Total: " << pairs.size() << ")" << std::endl;

    // ---------------------------------------------------------
    // 5. 导出裂缝-裂缝交互微观线段 (F-F Interaction)
    // ---------------------------------------------------------
    std::string fn_ff = caseName + "_Diagnostic_FF_Interactions.csv";
    std::ofstream out_ff(fn_ff);

    out_ff << "FF_Obj_ID" << sep
        << "FracID_1" << sep << "FracID_2" << sep
        << "Seg_Length" << sep
        << "Frac_CellID_1" << sep << "Frac_CellSolverIndex_1" << sep
        << "Frac_CellID_2" << sep << "Frac_CellSolverIndex_2" << "\n";

    const auto& ffObjs = mgr.fracture_network().ffIntersections;

    for (const auto& obj : ffObjs)
    {
        for (const auto& seg : obj.segments)
        {
            out_ff << obj.id << sep
                << obj.fracID_1 << sep << obj.fracID_2 << sep
                << seg.length << sep
                << seg.cellID_1 << sep << seg.solverIndex_1 << sep // 直接读取！
                << seg.cellID_2 << sep << seg.solverIndex_2 << "\n"; // 直接读取！
        }
    }
    out_ff.close();
    std::cout << "   -> Exported: " << fn_ff << std::endl;

    std::cout << "[Benchmark] Diagnostic Export Completed Successfully.\n" << std::endl;
}

// =========================================================
// 自洽性检查逻辑实现
// =========================================================

bool perform3DTopologyConsistencyCheck(MeshManager_3D& mgr)
{
    std::cout << "\n>>> Performing 3D Topology Consistency Check..." << std::endl;
    bool all_passed = true;
    int error_count = 0;

    // 1. 检查所有裂缝单元是否都有有效的 SolverIndex 和 ParentID
    const auto& fractures = mgr.fracture_network().getFractures();
    for (const auto& frac : fractures) {
        for (const auto& elem : frac.fracCells) {
            // Check 1: SolverIndex
            if (elem.solverIndex == -1) {
                std::cerr << "[FAIL] Frac " << frac.id << " Elem " << elem.id << " has SolverIndex = -1" << std::endl;
                all_passed = false; error_count++;
            }
            // Check 2: ParentID
            if (elem.parentFractureID != frac.id) {
                std::cerr << "[FAIL] Frac " << frac.id << " Elem " << elem.id
                    << " has wrong ParentID: " << elem.parentFractureID << std::endl;
                all_passed = false; error_count++;
            }
            // Check 3: Map Consistency
            if (frac.getElemIndex(elem.id) != elem.localIndex) {
                std::cerr << "[FAIL] Frac " << frac.id << " Elem " << elem.id << " Map inconsistent." << std::endl;
                all_passed = false; error_count++;
            }
        }
    }

    // 2. 检查 InteractionPair 的索引一致性 (Matrix <-> Pair <-> Fracture)
    const auto& pairs = mgr.getInteractionPairs();
    for (const auto& pair : pairs) {
        // 检查 Matrix Index
        int mIdx = mgr.mesh().getCellIndex(pair.matrixCellGlobalID);
        if (mIdx != pair.matrixSolverIndex) {
            std::cerr << "[FAIL] MF Pair: Matrix Index mismatch! GID=" << pair.matrixCellGlobalID << std::endl;
            all_passed = false; error_count++;
        }

        // 检查 Fracture Index
        const Fracture_2D* f = mgr.findFractureByID(mgr.fracture_network(), pair.fracMacroID);
        if (f) {
            int fLid = f->getElemIndex(pair.fracElementGlobalID);
            if (fLid != -1) {
                int fSid = f->fracCells[fLid].solverIndex;
                if (fSid != pair.fracCellSolverIndex) {
                    std::cerr << "[FAIL] MF Pair: Frac Index mismatch! GID=" << pair.fracElementGlobalID
                        << " TrueSID=" << fSid << " Stored=" << pair.fracCellSolverIndex << std::endl;
                    all_passed = false; error_count++;
                }
            }
            else {
                std::cerr << "[FAIL] MF Pair: Frac Elem GID " << pair.fracElementGlobalID << " not found!" << std::endl;
                all_passed = false; error_count++;
            }
        }
        else {
            std::cerr << "[FAIL] MF Pair: Frac Macro ID " << pair.fracMacroID << " not found!" << std::endl;
            all_passed = false; error_count++;
        }
    }

    if (all_passed) std::cout << "   [PASS] All 3D Topology Consistency Checks Passed." << std::endl;
    else std::cout << "   [FAIL] Found " << error_count << " inconsistencies." << std::endl;

    return all_passed;
}

// =========================================================
// Benchmark 主程序实现
// =========================================================

void run_Benchmark_3D_EDFM(
    const std::string& caseName,
    const std::vector<Fracture_2D>& inputFractures,
    int nU, int nV,
    MeshManager_3D::IntersectionStrategy strategy
)
{
    std::cout << "\n\n#########################################################" << std::endl;
    std::cout << "Running 3D-EDFM BENCHMARK Test Case: " << caseName << std::endl;

    // 0. 参数设置 (模拟一个标准基岩块)
    double lengthX = 100.0, lengthY = 100.0, lengthZ = 100.0;
    int nx = 1, ny = 1, nz = 1;

    // 1. 生成基岩 (Matrix Mesh)
    std::cout << " 1. Generating 3D Matrix Mesh (" << nx << "x" << ny << "x" << nz << ")..." << std::endl;
    MeshManager_3D mgr(lengthX, lengthY, lengthZ, nx, ny, nz, true, false);
    mgr.BuildSolidMatrixGrid_3D();

    // 2. 添加并划分裂缝网格 (Add & Mesh Fractures)
    std::cout << " 2. Adding " << inputFractures.size() << " Fractures and Meshing (" << nU << "x" << nV << ")..." << std::endl;
    for (const auto& f : inputFractures) {
        // 通过 Manager 添加裂缝到网络
        mgr.addFracturetoFractureNetwork(f);
    }

    // 对所有裂缝进行独立网格划分
    // 此时会自动填充 localIndex 和 parentID
    mgr.meshAllFracturesinNetwork(nU, nV, NormalVectorCorrectionMethod::OrthogonalCorrection);

    // 3. 构建全局索引 (Setup Global Indices) - [关键步骤!]
    // 这一步必须在求交之前调用，用于给基岩网格和裂缝element分配 SolverIndex
    std::cout << " 3. Setting up Global Solver Indices..." << std::endl;
	mgr.setupGlobalIndices();
    // 4. M-F 几何求交 (Intersection)
    std::cout << " 4. Running M-F Intersection (Strategy: " << (int)strategy << ")..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    mgr.SolveIntersection3D_improved_twist_accleration(strategy);

    // 5. F-F 几何求交
    mgr.DetectFractureFractureIntersectionsInNetwork(FFIntersectionStrategy::Octree_Optimized);

    // 6. 重构裂缝边属性 (Rebuild Edge Properties) - [关键步骤!]
    // 这一步生成 FVM 所需的 globalEdges，并填充 SolverIndex
    std::cout << " 6. Rebuilding Edge Properties for FVM..." << std::endl;  
    mgr.fracture_network().rebuildEdgeProperties();

    // =========================================================
    //  7. 构建拓扑映射 (关键步骤!)
    // =========================================================
    // 在导出诊断或运行 Solver 之前，构建加速索引
    // 这将填充 mat2InteractionMap_ 和 frac2InteractionMap_
    std::cout << " 7. Building Topology Maps..." << std::endl;
    mgr.buildTopologyMaps();
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "    -> Geometry & Topology Processing Time: " << elapsed.count() << " s" << std::endl;
    
    // 8. 导出诊断与检查
    std::cout << " 8. Exporting Diagnostics..." << std::endl;
    export3DTopologyDiagnosticCSV(mgr, caseName);

    perform3DTopologyConsistencyCheck(mgr);

    // 9. 导出txt结果用于matlab可视化
    mgr.exportMeshInfortoTxt(caseName);

    mgr.exportFracturesNetworkInfortoTxt_improved(caseName);

    mgr.exportInteractionPolygonsToTxt_improved(caseName);

    std::cout << "Benchmark Completed." << std::endl;
    std::cout << "#########################################################\n" << std::endl;
}

/**
 * @brief 快速创建裂缝对象的辅助函数
 * @param id 裂缝 ID
 * @param p1, p2, p3, p4 裂缝的四个角点 (逆时针顺序需满足围成环)
 * @return Fracture_2D 对象
 */
Fracture_2D CreateFracture1(int id, const Vector& p1, const Vector& p2, const Vector& p3, const Vector& p4)
{
    std::vector<Vector> boundaryVertices = { p1, p2, p3, p4 };
    // [Assumption]: Fracture_2D 构造函数接受 ID 和顶点列表
    return Fracture_2D(id, boundaryVertices);
}

int Improved_3D_EDFM_MeshTest()
{
    // -----------------------------------------------------
    // Case 1: 单一非扭曲裂缝 (Single Planar)
    // -----------------------------------------------------
    // 垂直平面 x=50, y=[20,80], z=[10,90]
    std::vector<Fracture_2D> case1_fracs;
    case1_fracs.push_back(CreateFracture1(1,
        Vector(50, 10, 10), Vector(50, 90, 10),
        Vector(50, 90, 90), Vector(50, 10, 90)));
	run_Benchmark_3D_EDFM("Test/MeshTest/GeomIndexTest/3D_EDFM/Benchmark_Case1_SinglePlanar", case1_fracs, 2, 2, MeshManager_3D::IntersectionStrategy::Rasterization_14DOP);


    // -----------------------------------------------------
    // Case 2: 单一扭曲裂缝 (Single Twisted)
    // -----------------------------------------------------
    // 基础为平面 z=50，但将一个角点拉高到 z=90，形成双曲抛物面
    // P1(20,20,50), P2(80,20,50), P3(80,80,90) <--- Twisted, P4(20,80,50)
    std::vector<Fracture_2D> case2_fracs;
    case2_fracs.push_back(CreateFracture1(1,
        Vector(20, 20, 22), Vector(80, 20, 22),
        Vector(80, 80, 90), Vector(20, 80, 22)));
    run_Benchmark_3D_EDFM("Test/MeshTest/GeomIndexTest/3D_EDFM/Benchmark_Case2_SingleTwisted", case2_fracs, 10, 10, MeshManager_3D::IntersectionStrategy::Rasterization_14DOP);

    // -----------------------------------------------------
    // Case 3: 两根扭转相交裂缝 (Multi Intersecting)
    // -----------------------------------------------------
    std::vector<Fracture_2D> case3_fracs;
    case3_fracs.push_back(CreateFracture1(1,
        Vector(20, 10, 10), Vector(80, 10, 75),
        Vector(80, 70, 10), Vector(20, 80, 45)));
    case3_fracs.push_back(CreateFracture1(2,
        Vector(55, 10, 10), Vector(55, 70, 10),
        Vector(30, 70, 88), Vector(70, 10, 68)));
    run_Benchmark_3D_EDFM("Test/MeshTest/GeomIndexTest/3D_EDFM/Benchmark_Case3_TwoVetrix", case3_fracs, 2, 2, MeshManager_3D::IntersectionStrategy::Rasterization_14DOP);

    // -----------------------------------------------------
    // Case 4: 多根相交裂缝 (Multi Intersecting)
    // -----------------------------------------------------
    // 十字交叉: Fracture A (x=50), Fracture B (y=50)
    std::vector<Fracture_2D> case4_fracs;
    case4_fracs.push_back(CreateFracture(1,
        Vector(50, 10, 10), Vector(50, 90, 10),
        Vector(50, 90, 90), Vector(50, 10, 90))); // Plane X=50
    case4_fracs.push_back(CreateFracture(2,
        Vector(10, 50, 10), Vector(90, 50, 10),
        Vector(90, 50, 90), Vector(10, 50, 90))); // Plane Y=50
    run_Benchmark_3D_EDFM("Test/MeshTest/GeomIndexTest/3D_EDFM/Benchmark_Case4_MultiIntersecting", case4_fracs, 2, 2, MeshManager_3D::IntersectionStrategy::Rasterization_14DOP);


    return 0;

}