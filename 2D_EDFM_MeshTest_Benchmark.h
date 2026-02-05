/**
 * @file EDFM_Benchmark.h
 * @brief EDFM 拓扑关系基准测试与诊断工具 (Enhanced Version)
 * @details 包含基岩、裂缝及网格面的全拓扑关系导出与自洽性检查
 * @author Computational Physicist
 */

#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <map>
#include <set>
#include <algorithm>

#include "MeshManager.h"

 // ————————————————————————————————————————————————————————————————————————
 // 辅助工具
 // ————————————————————————————————————————————————————————————————————————
/**
 * @brief 将 vector 转换为 CSV 友好的字符串格式 [1;2;3]
 */
template <typename T>
std::string vec_to_csv_string(const std::vector<T>& v) {
    if (v.empty()) return "[]";
    std::stringstream ss;
    ss << "[";
    for (size_t i = 0; i < v.size(); ++i) {
        ss << v[i] << (i == v.size() - 1 ? "" : ";");
    }
    ss << "]";
    return ss.str();
}

/**
 * @brief 将 set 转换为 CSV 友好的字符串格式
 */
template <typename T>
std::string set_to_csv_string(const std::set<T>& s) {
    if (s.empty()) return "[]";
    std::vector<T> v(s.begin(), s.end());
    return vec_to_csv_string(v);
}

// ————————————————————————————————————————————————————————————————————————
// 核心诊断导出函数
// ————————————————————————————————————————————————————————————————————————

void exportTopologyDiagnosticCSV(MeshManager& mgr, const std::string& caseName)
{
    std::string sep = ",";//CSV 分隔符
    // =========================================================
    // 0. 预计算与数据引用 (Pre-computation)
    // =========================================================

    // =========================================================
    // 1. 导出基岩单元拓扑报告 (Matrix Topology Report)
    // =========================================================
    std::string fn_matrix = caseName + "_Diagnostic_Matrix.csv";
    std::ofstream out_m(fn_matrix);

    if (!out_m.is_open()) {
        std::cerr << "[Error] Unable to create file: " << fn_matrix << std::endl;
        return;
    }


    out_m << "LocalIndex" << sep            // [查表值] 用于验证 Map 一致性
        << "SolverIndex" << sep           // [查表值] 
        << "GlobalID" << sep              // Gmsh Tag
        << "Center_X" << sep << "Center_Y" << sep << "Center_Z" << sep
        << "FaceIDs_Locals" << sep
        << "FaceIDs_Globlas" << sep
        << "Neighbor_LocalIndices" << sep // 邻居 SolverIndex
        << "Neighbor_GlobalIDs" << sep    // 邻居 Gmsh Tag
        << "Embedded_FracElem_SolverIndices" << sep
        << "Embedded_FracElem_LocalIndices" << sep
        << "Contained_ParentFracIDs"
        << "\n";

    const auto& cells = mgr.mesh().getCells();
    const auto& faces = mgr.mesh().getFaces();

    for (size_t i = 0; i < cells.size(); ++i)
    {
		const auto& cell = cells[i];

        int mappedLocalIndex = -1;
        try {
			mappedLocalIndex = mgr.mesh().getCellId2Index().at(cell.id); //获取 LocalIndex
        }
        catch (...) {
            std::cerr << "[Error] Cell " << cell.id << " not found in cellId2Index!" << std::endl;
        }

        // --- A. 收集邻居信息 ---
        std::vector<int> faceIDs_Locals;
        std::vector<int> faceIDs_Globals;
        std::vector<int> neigh_Locals;
		std::vector<int> neigh_Globals;

        for (int fid : cell.CellFaceIDs)
        {
			int faceVecIdx = mgr.mesh().getFaceIndex(fid); // 从 GlobalID 映射到 faces 数组索引
            if (faceVecIdx == -1 || faceVecIdx >= static_cast<int>(faces.size())) {
                std::cerr << "[Error] Face ID " << fid << " in Cell " << cell.id
                    << " maps to invalid index " << faceVecIdx << std::endl;
                continue; // 跳过坏数据
            }

            const auto& face = faces[faceVecIdx]; // 安全访问
			faceIDs_Globals.push_back(face.id);                             // GlobalID
			faceIDs_Locals.push_back(faceVecIdx);    // LocalIndex

            int myIdx = mappedLocalIndex;
			int neighborIdx = (face.ownerCell_index == myIdx) ? face.neighborCell_index : face.ownerCell_index;

			neigh_Locals.push_back(neighborIdx); // SolverIndex / LocalIndex

            if (neighborIdx >= 0 && neighborIdx < static_cast<int>(cells.size())) {
				neigh_Globals.push_back(cells[neighborIdx].id); // GlobalID
            }
            else {
                neigh_Globals.push_back(-1);
            }

        }

        // --- B. 收集每个cell内的嵌入裂缝信息 ---
		std::vector<int> embeddedSegs = 
            mgr.mesh().getFracElemSolverIndexfromCellGlobalId(cell.id); //输入 GlobalID，输出 FracElem SolverIndex
       
        std::set<int> parentFracIDs;
		std::set<int> embeddedLocalIndices;
        for (int segIdx : embeddedSegs) {
            // [使用新接口] 通过 SolverIndex 查找 Element 指针
            const FractureElement* elemPtr = 
                mgr.getFractureElementBySolverIndex(segIdx);
            if (elemPtr != nullptr) {
                parentFracIDs.insert(elemPtr->parentFractureID);
				embeddedLocalIndices.insert(elemPtr->id); // LocalIndex
            }
        }

        out_m<<mappedLocalIndex<<sep
            <<mappedLocalIndex<<sep
            <<cell.id<<sep
			<< cell.center.m_x << sep << cell.center.m_y << sep<<cell.center.m_z<<sep
			<< vec_to_csv_string(faceIDs_Locals) << sep
			<< vec_to_csv_string(faceIDs_Globals) << sep
			<< vec_to_csv_string(neigh_Locals) << sep
			<< vec_to_csv_string(neigh_Globals) << sep
			<< vec_to_csv_string(embeddedSegs) << sep
			<< set_to_csv_string(embeddedLocalIndices) << sep
			<< set_to_csv_string(parentFracIDs)
			<< "\n";
    }
	out_m.close();
    std::cout << "   -> [Output] Matrix Diagnostic CSV: " << fn_matrix << std::endl;

    // =========================================================
    // 2. 导出网格面拓扑报告 (Face Topology Report)
    // =========================================================
    std::string fn_face = caseName + "_Diagnostic_Face.csv";
    std::ofstream out_face(fn_face);

    out_face << "LocalIndex" << sep
        << "GlobalID" << sep
        << "Center_X" << sep << "Center_Y" << sep << "Center_Z" << sep
        << "Owner_GlobalID" << sep
        << "Neighbor_GlobalID" << sep
        << "Owner_LocalIndex" << sep
        << "Neighbor_LocalIndex"
        << "\n";

    for (size_t i = 0; i < faces.size(); ++i)
    {
        const auto& f = faces[i];

        int owner_GID = (f.ownerCell_index >= 0) ? cells[f.ownerCell_index].id : -1;
		int neighbor_GID = (f.neighborCell_index >= 0) ? cells[f.neighborCell_index].id : -1;

		int face_LocalIndex = mgr.mesh().getFaceIndex(f.id);

        out_face << face_LocalIndex  << sep
            << f.id << sep
            << f.midpoint.m_x<< sep << f.midpoint.m_y << sep << f.midpoint.m_z << sep
            << owner_GID << sep
            << neighbor_GID << sep
            << f.ownerCell_index << sep
            << f.neighborCell_index
            << "\n";
    }
    out_face.close();
    std::cout << "   -> [Output] Face Diagnostic CSV:   " << fn_face << std::endl;

    // =========================================================
    // 3. 导出裂缝单元拓扑报告 (Fracture Topology Report)
    // =========================================================
    std::string fn_frac = caseName + "_Diagnostic_Fracture.csv";
    std::ofstream out_f(fn_frac);

    out_f << "SolverIndex" << sep
        << "ParentFracID" << sep
        << "ElemLocalID" << sep
        << "HostCell_GlobalID" << sep
        << "HostCell_SolverIndex" << sep
        << "Length" << sep << "Aperture" << sep<<"cellToFracture_Distance"<< sep
        << "Start_X" << sep << "Start_Y" << sep
        << "End_X" << sep << "End_Y" << sep
        << "Intersections" // [新增] 裂缝相交信息列
        << "\n";

    const auto& fractures = mgr.fracture_network().fractures;

    for (size_t fid = 0; fid < fractures.size(); ++fid)
    {
        const auto& frac = fractures[fid];
        for (const auto& elem : frac.elements)
        {
            Vector fracDir = frac.end - frac.start;
            Vector p_start = frac.start + fracDir * elem.param0;
            Vector p_end = frac.start + fracDir * elem.param1;

            int hostSolverIdx = -1;
            auto it = mgr.mesh().getCellId2Index().find(elem.cellID);
            if (it != mgr.mesh().getCellId2Index().end()) {
				hostSolverIdx = it->second;
            }

            // [新增] 提取相交信息 (Fracture-Fracture Intersections)
            std::vector<std::string> intersectInfos;
            for (const auto& nb : elem.neighbors)
            {
                // 格式: PeerSolverIndex(PeerFracID:PeerLocalID @IntersectionID)
                std::stringstream ss;
                ss << nb.solverIndex << "("
                    << nb.fracID << ":"
                    << nb.elemID << " @"
                    << nb.intersectionID << ")";
                intersectInfos.push_back(ss.str());
            }

			out_f << elem.solverIndex << sep
				<< elem.parentFractureID << sep
				<< elem.id << sep
				<< elem.cellID << sep
				<< hostSolverIdx << sep
				<< elem.length << sep
				<< elem.aperture << sep
                << elem.avgDistance<<sep
				<< p_start.m_x << sep << p_start.m_y << sep
                << p_end.m_x << sep << p_end.m_y << sep
                << vec_to_csv_string(intersectInfos) 
				<< "\n";
        }
    }
    out_f.close();
    std::cout << "   -> [Output] Fracture Diagnostic CSV: " << fn_frac << std::endl;
}

// ————————————————————————————————————————————————————————————————————————
// 自洽性检查逻辑
// ————————————————————————————————————————————————————————————————————————

bool performTopologyConsistencyCheck(MeshManager& mgr)
{
    std::cout << "\n   >>> Performing Internal Topology Self-Check..." << std::endl;
    bool all_passed = true;
    int error_count = 0;

    struct ElemRef { int fid; int localID; int cellID; };
    std::map<int, ElemRef> solverIdxToElem;

    const auto& fractures = mgr.fracture_network().fractures;
    for (size_t fid = 0; fid < fractures.size(); ++fid) {
        for (const auto& elem : fractures[fid].elements) {
            if (elem.solverIndex != -1) {
                solverIdxToElem[elem.solverIndex] = { (int)fid, elem.id, elem.cellID };
            }
        }
    }

    const auto& cells = mgr.mesh().getCells();
    for (size_t i = 0; i < cells.size(); ++i)
    {
        const auto& cell = cells[i];
        std::vector<int> embeddedIndices = mgr.mesh().getFracElemSolverIndexfromCellGlobalId(cell.id);

        for (int solverIdx : embeddedIndices)
        {
            if (solverIdxToElem.find(solverIdx) == solverIdxToElem.end()) {
                std::cerr << "      [FAIL] Matrix Cell " << cell.id << " claims to contain SolverIndex "
                    << solverIdx << ", but it does not exist in FractureNetwork!" << std::endl;
                all_passed = false; error_count++;
                continue;
            }

            ElemRef ref = solverIdxToElem[solverIdx];
            if (ref.cellID != cell.id) {
                std::cerr << "      [FAIL] Topology Mismatch! Matrix Cell " << cell.id
                    << " claims FractureElem (SolverIdx=" << solverIdx << "), "
                    << "but that Element thinks its Host is Cell " << ref.cellID << std::endl;
                all_passed = false; error_count++;
            }
        }
    }

    if (all_passed) {
        std::cout << "      [PASS] Topology Consistency Check Passed." << std::endl;
    }
    else {
        std::cout << "      [FAIL] Found " << error_count << " topology inconsistencies!" << std::endl;
    }
    return all_passed;
}
// ————————————————————————————————————————————————————————————————————————
// Benchmark 主程序
// ————————————————————————————————————————————————————————————————————————

void run_Benchmark_Topology_Test(
    const std::string& caseName,
    const std::vector<Fracture>& fractures,
    const NormalVectorCorrectionMethod& faceVector_method,
    const DistanceMetric& dis_method,
    const IntersectionSearchStrategy_2D& search_strategy
)
{
    std::cout << "\n\n#########################################################" << std::endl;
    std::cout << "Running BENCHMARK Test Case: " << caseName << std::endl;

    double lengthX = 1, lengthY = 1, lengthZ = 0;
    int sectionNumX = 2, sectionNumY = 2, sectionNumZ = 0;

    // 1. 生成基岩
    std::cout << " 1. Generating Matrix Mesh..." << std::endl;
    MeshManager mgr(lengthX, lengthY, lengthZ, sectionNumX, sectionNumY, sectionNumZ, true, false);
    mgr.BuildSolidMatrixGrid_2D(faceVector_method);
    mgr.exportMesh(caseName + "_MatrixMesh");

    // 2. 添加裂缝
    std::cout << " 2. Adding " << fractures.size() << " Fractures..." << std::endl;
    for (const auto& f : fractures) {
        mgr.addFracture(f.start, f.end);
    }

    // 3. 几何切割
    std::cout << " 3. Subdividing Fractures..." << std::endl;
    mgr.setDistanceMetric(dis_method);
    auto start = std::chrono::high_resolution_clock::now();
    mgr.DetectAndSubdivideFractures(search_strategy);

    // 4. 构建索引 (System Indexing)
    // 注意：BuildGlobalSystemIndexing 内部已调用 rebuildFractureMap
    std::cout << " 4. Building Global System Indexing..." << std::endl;
    mgr.BuildGlobalSystemIndexing();
	mgr.BuildFracturetoFractureTopology();
    auto end = std::chrono::high_resolution_clock::now();

    // 5. 导出
    mgr.ComputeFractureGeometryCouplingCoefficient();
    mgr.exportFractures(caseName + "_FractureMesh");

    std::cout << " 5. Exporting Benchmark Diagnostics..." << std::endl;
    exportTopologyDiagnosticCSV(mgr, caseName);

    performTopologyConsistencyCheck(mgr);

    std::cout << "Benchmark Completed." << std::endl;
    std::cout << "#########################################################\n" << std::endl;
}

int EDFM_DFN_Geomtest_2D()
{
    // =========================================================
    // Step 1: 预生成统一的随机裂缝网络 (DFN)
    // =========================================================
    std::cout << ">>> Generating Random DFN Network..." << std::endl;

    // 创建一个临时的 Manager 仅用于生成裂缝数据
    // 尺寸必须与 run2D_EDFM_test 中的 lengthX, lengthY 一致 (1.0 x 1.0)
    double Lx = 1.0, Ly = 1.0;
    MeshManager genMgr(Lx, Ly, 0.0, 1, 1, 0, false, false);

    // DFN 参数设置
    int N_fractures = 5;           // 裂缝数量 (建议 50-200 以测试性能)
    unsigned int seed = 2024;       // 固定种子，保证每次运行结果一致
    double Lmin = 0.7, Lmax = 0.8;  // 裂缝长度范围
    double alpha = 2.5;             // 幂律指数
    double kappa = 0.0;             // 0.0 表示完全随机方向
    bool avoidOverlap = true;       // 【关键】开启避免重合检测

    genMgr.setDFNRandomSeed(seed);
    genMgr.generateDFN(
        N_fractures,
        Vector(0, 0, 0),   // 最小坐标 (基岩左下角)
        Vector(Lx, Ly, 0), // 最大坐标 (基岩右上角)
        Lmin, Lmax, alpha, kappa, avoidOverlap
    );

    // 提取生成的裂缝列表
    // 注意：这里我们拷贝一份 fracture_network 中的裂缝数据
    std::vector<Fracture> dfnFractures = genMgr.fracture_network().fractures;

    if (dfnFractures.empty()) {
        std::cerr << "[Error] DFN generation failed or generated 0 fractures." << std::endl;
        return -1;
    }
    std::cout << ">>> DFN Generated. Total fractures: " << dfnFractures.size() << "\n" << std::endl;


    // =========================================================
    // Step 2: 使用同一套 DFN 数据运行对比测试
    // =========================================================

	run_Benchmark_Topology_Test(
		"Test/MeshTest/GeomIndexTest/2D_EDFM/Benchmark_2D_EDFM_Topology_Test",
		dfnFractures,
		NormalVectorCorrectionMethod::OrthogonalCorrection,
		DistanceMetric::CrossAwareGauss,
		IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA
	);
	return 0;

}