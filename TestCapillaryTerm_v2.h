#pragma once
#include <algorithm>
#include <string>
#include <vector>
#include <iostream>
#include "MeshManager.h"
#include "PhysicalPropertiesManager.h"
#include "Initializer.h"
#include "PostProcessor.h"
#include "IMPES_Iteration_Loop.h"

/**
 * @brief 单独测试毛细扩散算子 L[Pc] = ∇·(k λ_g ρ_g ∇Pc)
 *
 * 测试设置：
 * - kxx=kyy=kzz=1，λ_g=1，ρ_g=1，λ_w=0，ρ_w=0 ⇒ L[Pc] = ∇² Pc；
 * - 构造解析 Pc(x,y) = P0 + α (y - ymin)，满足 ∇²Pc = 0；
 * - 左右边界：Neumann 零通量；上下边界：解析 Dirichlet；
 * - 用 TwoPhaseDiffusionTemplate + assemble_COO + applyCOO 得到 y=A*Pc - b，
 *   检查 max|y|、L2|y| 是否接近 0。
 *
 * @return 0 表示流程执行成功（是否“物理上正确”看终端残差大小），非 0 表示流程错误。
 */
 /// \brief 独立测试“毛细压扩散算子”的零模（常数 Pc）性质.
 ///        期望：在均匀 k、均匀 rho_cap、常数 Pc 条件下，
 ///             离散算子 L[Pc] ≈ 0（只剩机器误差）.
int TestCapillaryTerm_v2()
{
    using namespace IMPES_Iteration;

    // ---------------- 0. 网格与场注册 ----------------
    const double lengthX = 100.0;
    const double lengthY = 100.0;
    const double lengthZ = 0.0;   // 2D
    const int    sectionNumX = 50;
    const int    sectionNumY = 50;
    const int    sectionNumZ = 0;
    const bool   usePrism = true;
    const bool   useQuadBase = false;

    std::cout << "\n===== [TEST] Capillary diffusion operator (v2, constant Pc) =====\n";

    MeshManager mgr(lengthX, lengthY, lengthZ,
        sectionNumX, sectionNumY, sectionNumZ,
        usePrism, useQuadBase);
    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection);

    Mesh& mesh = mgr.mesh();
    FieldRegistry     reg;
    FaceFieldRegistry freg;

    const auto& cells = mesh.getCells();
    const auto& id2idx = mesh.getCellId2Index();
    const size_t nCells = cells.size();

    // ---------------- 1. 定义压力与系数字段名 ----------------
    PressureAssemblyConfig cfg;
    cfg.operator_tag = "Pc_capillary_test";  // 用于 OperatorFieldNames 前缀
    cfg.Pc_field = "Pc_test";            // 当前测试用的毛细压场
    cfg.rho_capillary_field = "rho_capillary_mass"; // 可随意命名，和 cfg 一致即可
    cfg.gravity_dummy_field = "gravity_dummy_scalar"; // 这里不用重力，只做占位
    cfg.gravity = Vector{ 0.0, 0.0, 0.0 };
    cfg.enable_buoyancy = false;
    cfg.gradient_smoothing = 0;

    // ---------------- 2. 构建 Pc 场 & k 场 & rho_cap 场 ----------------
    // 2.1: 常数 Pc 场
    auto PcF = reg.getOrCreate<volScalarField>(cfg.Pc_field, nCells, 0.0);

    // 2.2: 渗透率场（各向同性 k = 1.0，实际值只影响尺度）
    auto kxx = reg.getOrCreate<volScalarField>("kxx", nCells, 1.0);
    auto kyy = reg.getOrCreate<volScalarField>("kyy", nCells, 1.0);
    auto kzz = reg.getOrCreate<volScalarField>("kzz", nCells, 1.0);

    // 2.3: “毛细系数” rho_capillary，用来承载 λ_g ρ_g 的等效系数
    auto rho_cap = reg.getOrCreate<volScalarField>(cfg.rho_capillary_field, nCells, 0.0);
    auto rho_dummy = reg.getOrCreate<volScalarField>(cfg.gravity_dummy_field, nCells, 0.0);

    const double Pc0 = 2.0e5;  // 任意常数 Pc [Pa]
    const double rhoCapVal = 1.0;    // 等效 λ_g ρ_g；只要正即可
    const double kVal = 1.0;    // 等效 k；只要正即可

    for (const auto& c : cells)
    {
        if (c.id < 0) continue;
        const size_t i = id2idx.at(c.id);

        (*PcF)[i] = Pc0;
        (*kxx)[i] = kVal;
        (*kyy)[i] = kVal;
        (*kzz)[i] = kVal;
        (*rho_cap)[i] = rhoCapVal;
        (*rho_dummy)[i] = 0.0;  // 不使用，只是占位
    }

    // ---------------- 3. 边界条件（毛细项：这里完全不设 BC，相当于自然 Neumann=0） ----------------
    PressureBC::Registry pbc_pc;  // 空 registry，不添加任何边界
    PressureBCAdapter    PbcA{ pbc_pc };

    // ---------------- 4. 构建未知量编号映射 ----------------
    int N = 0;
    auto lid_of_cell = buildUnknownMap(mesh, N);
    if (N <= 0)
    {
        std::cerr << "[CapillaryTest_v2] buildUnknownMap failed: N=" << N << "\n";
        return -1;
    }

    // ---------------- 5. 构造 OperatorFieldNames 和 mobility_tokens ----------------
    // 只需要 kxx/kyy/kzz 这三个渗透率场；rho 部分由 rho_capillary_field 承载
    auto nm = makeNames(cfg.operator_tag); // a_f_diff, s_f_diff 等名字

    std::vector<std::string> mobility_tokens{
        "kxx:kxx",
        "kyy:kyy",
        "kzz:kzz"
    };

    // ---------------- 6. 收集 Pc 场到向量 Pc_vec ----------------
    std::vector<double> Pc_vec, Lvec;
    if (!detailed::gatherFieldToVector(mesh, reg, cfg.Pc_field,
        lid_of_cell, N, Pc_vec))
    {
        std::cerr << "[CapillaryTest_v2] gatherFieldToVector failed for Pc.\n";
        return -2;
    }

    // ---------------- 7. 构建算子并作用：Lvec = A * Pc_vec - b ----------------
    if (!detailed::buildAndApplyOperator(
        mgr, reg, freg,
        PbcA,
        nm,
        mobility_tokens,
        cfg.rho_capillary_field,  // rho_coeff_field → 用毛细系数场
        cfg.gravity_dummy_field,  // rho_buoy_field  → 不用
        cfg.Pc_field,             // x_field → Pc
        cfg,
        /*enable_buoyancy=*/false,
        Pc_vec,
        Lvec))
    {
        std::cerr << "[CapillaryTest_v2] buildAndApplyOperator failed.\n";
        return -3;
    }

    // ---------------- 8. 计算残差范数并输出 ----------------
    double maxAbs = 0.0;
    double l2sum = 0.0;
    int    nUsed = 0;

    int    sampleCellId = -1;
    double samplePc = 0.0;
    double sampleL = 0.0;

    for (const auto& c : cells)
    {
        if (c.id < 0) continue;
        const size_t i = id2idx.at(c.id);
        const int    lid = lid_of_cell[i];
        if (lid < 0) continue;

        const double val = std::fabs(Lvec[lid]);
        maxAbs = std::max(maxAbs, val);
        l2sum += val * val;
        ++nUsed;

        // 记录一个样本单元
        if (sampleCellId < 0)
        {
            sampleCellId = c.id;
            samplePc = (*PcF)[i];
            sampleL = Lvec[lid];
        }
    }

    const double l2norm = (nUsed > 0) ? std::sqrt(l2sum / nUsed) : 0.0;
    const double relMax = (std::fabs(Pc0) > 0.0) ? maxAbs / std::fabs(Pc0) : 0.0;

    std::cout << "[CapillaryTest_v2] N = " << nUsed
        << ", max|L[Pc]| = " << maxAbs
        << ", L2|L[Pc]| = " << l2norm
        << ", relMax = " << relMax << "\n";

    if (sampleCellId >= 0)
    {
        std::cout << "  sample cell id = " << sampleCellId
            << ", Pc = " << samplePc
            << ", L[Pc] = " << sampleL << "\n";
    }

    std::cout << "[CapillaryTest_v2] finished. "
        "Expect residuals near machine precision if capillary operator is assembled correctly.\n";

    return 0;
}