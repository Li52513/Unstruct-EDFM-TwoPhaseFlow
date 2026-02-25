/**
 * @file Test_DOFMapper.cpp
 * @brief DOF编排器及交织存储策略的生产级独立验证程序
 * @details 包含完整的 2D/3D 网格生成、裂缝拓扑构建及全局稀疏矩阵维度/索引的严密断言测试。
 */

#include <iostream>
#include <cassert>
#include <vector>

 // 引入核心管理类
#include "MeshManager.h"
#include "3D_MeshManager.h"
#include "FractureCommon.h"

// =====================================================================
// 测试用例 1：2D EDFM 单相热流耦合场景验证 (P, T -> Ndof = 2)
// =====================================================================
/**
 * @brief 验证 2D 网格生成及单相场景下的自由度交织存储
 */
void Verify_2D_SinglePhase_DOFMapping()
{
	std::cout << "\n========== [Test 1] 2D MeshManager: Single-Phase (Ndof = 2) ==========" << std::endl;

	// 1. 实例化 2D MeshManager (域大小 10x10，网格 5x5)
	MeshManager mgr2D(10.0, 10.0, 0.0, 5, 5, 0, false, false);

    // 2. 真实生成基岩网格
    std::cout << "-> Generating 2D Matrix Grid..." << std::endl;
    mgr2D.BuildSolidMatrixGrid_2D();

    // 3. 真实添加裂缝并离散化
    std::cout << "-> Adding Fracture and subdividing..." << std::endl;
    mgr2D.addFracture(Vector(2.0, 2.0, 0.0), Vector(8.0, 8.0, 0.0));
    mgr2D.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::BruteForce);

    // 4. 构建系统全局索引
    mgr2D.BuildGlobalSystemIndexing();

    // 5. [核心] 设置单相热流耦合自由度 (P, T)
    mgr2D.setNumDOFs(2);

    // 6. 获取拓扑信息
    int matrixDOFCount = mgr2D.getMatrixDOFCount(); // 基岩块数
    int totalBlockCount = mgr2D.getTotalDOFCount(); // 总块数 (基岩+裂缝)
    int expectedTotalEquations = totalBlockCount * 2;

    // ------------------- [严密断言测试] -------------------

    // 断言 A: 基本属性正确
    assert(mgr2D.getNumDOFs() == 2 && "DOF per block should be 2 for single-phase.");
    assert(mgr2D.getTotalEquationDOFs() == expectedTotalEquations && "Total equation size mismatch!");

    // 断言 B: 越界安全拦截测试
    assert(mgr2D.getEquationIndex(-1, 0) == -1 && "Should reject negative block index.");
    assert(mgr2D.getEquationIndex(totalBlockCount, 0) == -1 && "Should reject out-of-bound block index.");
    assert(mgr2D.getEquationIndex(0, 2) == -1 && "Should reject out-of-bound DOF offset (>= 2).");

    // 断言 C: 基岩网格的交织排布验证
    if (matrixDOFCount > 0)
    {
        // 测试第一个基岩网格块 (Block 0)
        assert(mgr2D.getEquationIndex(0, 0) == 0 && "Matrix Cell 0, P -> Row 0");
        assert(mgr2D.getEquationIndex(0, 1) == 1 && "Matrix Cell 0, T -> Row 1");

        // 测试任意基岩块 (如 Block 5)
        if (matrixDOFCount > 5) {
            assert(mgr2D.getEquationIndex(5, 0) == 10 && "Matrix Cell 5, P -> Row 10");
            assert(mgr2D.getEquationIndex(5, 1) == 11 && "Matrix Cell 5, T -> Row 11");
        }
    }

    // 断言 D: 裂缝网格的偏移交织排布验证 (裂缝编号紧随基岩)
    if (totalBlockCount > matrixDOFCount)
    {
        int firstFracBlock = matrixDOFCount;
        int expectedRowP = firstFracBlock * 2;
        int expectedRowT = firstFracBlock * 2 + 1;

        assert(mgr2D.getEquationIndex(firstFracBlock, 0) == expectedRowP && "First Frac P row mismatch.");
        assert(mgr2D.getEquationIndex(firstFracBlock, 1) == expectedRowT && "First Frac T row mismatch.");
    }

    std::cout << "[PASS] 2D Single-Phase DOF Mapping verified successfully." << std::endl;
}

// =====================================================================
// 测试用例 2：3D EDFM 两相热流耦合场景验证 (Pw, Sn, T -> Ndof = 3)
// =====================================================================
/**
 * @brief 验证 3D 网格生成及两相场景下的自由度交织存储
 */
void Verify_3D_TwoPhase_DOFMapping()
{
    std::cout << "\n========== [Test 2] 3D MeshManager: Two-Phase (Ndof = 3) ==========" << std::endl;

    // 1. 实例化 3D MeshManager (域大小 10x10x10，网格 2x2x2)
    MeshManager_3D mgr3D(10.0, 10.0, 10.0, 2, 2, 2, true, false);

    // 2. 真实生成 3D 基岩网格
    std::cout << "-> Generating 3D Matrix Grid..." << std::endl;
    mgr3D.BuildSolidMatrixGrid_3D();

    // 3. 插入 2D 裂缝并网格化 
    std::cout << "-> Adding 2D Fracture and meshing..." << std::endl;

    // 构造一个合法的 2D 裂缝顶点容器 (逆时针排布的四边形，位于 Z=5.0 的水平面)
    // 确保网格划分器 (meshAllFracturesinNetwork) 能正常获取几何信息进行 2x2 划分
    std::vector<Vector> dummyCorners = {
        Vector(2.0, 2.0, 5.0),
        Vector(8.0, 2.0, 5.0),
        Vector(8.0, 8.0, 5.0),
        Vector(2.0, 8.0, 5.0)
    };

    // 使用带参构造函数进行严谨的初始化
    Fracture_2D frac(0, dummyCorners, 1.0, 0.001);
    mgr3D.addFracturetoFractureNetwork(frac);

    // 将宏观裂缝划分为微元 (2x2)
    mgr3D.meshAllFracturesinNetwork(2, 2);

    // 4. 构建全局索引 (通知基岩与裂缝进行连续编号分配)
    mgr3D.setupGlobalIndices();

    // 5. [核心] 设置两相热流耦合自由度 (Pw, Sn, T)
    mgr3D.setNumDOFs(3);

    // 6. 获取拓扑信息
    int matrixDOFCount = static_cast<int>(mgr3D.mesh().getCells().size());
    int totalBlockCount = mgr3D.getTotalDOFCount();
    int expectedTotalEquations = totalBlockCount * 3;

    // ------------------- [严密断言测试] -------------------

    // 断言 A: 基本属性正确
    assert(mgr3D.getNumDOFs() == 3 && "DOF per block should be 3 for two-phase.");
    assert(mgr3D.getTotalEquationDOFs() == expectedTotalEquations && "Total equation size mismatch!");

    // 断言 B: 越界防护
    assert(mgr3D.getEquationIndex(-1, 0) == -1 && "Should reject negative block index.");
    assert(mgr3D.getEquationIndex(0, 3) == -1 && "Should reject out-of-bound DOF offset (>= 3).");

    // 断言 C: 基岩交织排布验证
    if (matrixDOFCount > 0)
    {
        // Block 0: Pw(0), Sn(1), T(2)
        assert(mgr3D.getEquationIndex(0, 0) == 0);
        assert(mgr3D.getEquationIndex(0, 1) == 1);
        assert(mgr3D.getEquationIndex(0, 2) == 2);

        // Block 3: Pw(9), Sn(10), T(11)
        if (matrixDOFCount > 3) {
            assert(mgr3D.getEquationIndex(3, 0) == 9);
            assert(mgr3D.getEquationIndex(3, 1) == 10);
            assert(mgr3D.getEquationIndex(3, 2) == 11);
        }
    }

    // 断言 D: 裂缝交织排布验证
    if (totalBlockCount > matrixDOFCount)
    {
        int firstFracBlock = matrixDOFCount;
        assert(mgr3D.getEquationIndex(firstFracBlock, 0) == firstFracBlock * 3);
        assert(mgr3D.getEquationIndex(firstFracBlock, 1) == firstFracBlock * 3 + 1);
        assert(mgr3D.getEquationIndex(firstFracBlock, 2) == firstFracBlock * 3 + 2);
    }

    std::cout << "[PASS] 3D Two-Phase DOF Mapping verified successfully." << std::endl;
}

// =====================================================================
// 主测试入口
// =====================================================================
/**
 * @brief 执行 DOF Mapper 的全量自动化测试
 */
void RunDOFMapperVerification()
{
    std::cout << "=====================================================================" << std::endl;
    std::cout << "            STARTING DOF MAPPER VERIFICATION SUITE                   " << std::endl;
    std::cout << "=====================================================================" << std::endl;

    try {
        Verify_2D_SinglePhase_DOFMapping();
        Verify_3D_TwoPhase_DOFMapping();

        std::cout << "\n=====================================================================" << std::endl;
        std::cout << "          [SUCCESS] ALL DOF MAPPER TESTS PASSED SUCCESSFULLY!        " << std::endl;
        std::cout << "=====================================================================" << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "\n[FAILED] Exception caught during tests: " << e.what() << std::endl;
    }
}