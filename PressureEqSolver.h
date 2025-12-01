#include <vector>
#include <memory>
#include <iostream>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "BCAdapter.h"
#include "TwoPhaseWells_StrictRate.h"
#include "FVM_WellCoupling_TwoPhase.h"
#include "TimeTermAssemblerandSolver.h"
#include "SolverContrlStrName.h"
#include "LinearSolver_Eigen.h"
#include "PressureEqAssemblerandSolver.h"

namespace IMPES_Iteration
{
    //================= 控制参数 & 报告结构 =================//

    /**
     * \brief 压力求解器控制参数（单次外迭代 step）
     *
     * - assembly: 里面指定 p 字段名、Pc 字段名等
     * - linear:   线性求解器选项 (tol, maxIter, 等)
     * - under_relax: 对 cell 压力场的欠松弛系数 (1.0 = 无欠松弛)
     * - tol_abs/tol_rel: 推荐用于外迭代收敛判据 (dp_inf) 的阈值（在外层 time loop 使用）
     */
	struct PressureSolveControls
	{
        PressureAssemblyConfig assembly;   ///< 包含所有字符串字段名和 VG 参数
        LinearSolverOptions    linear;     ///< 线性求解器控制参数

        double under_relax = 1.0;   ///< 欠松弛系数 (0<urf<=1, 建议 0.3~0.7)
        double tol_abs = 1e-6; ///< dp_inf 绝对阈值（外迭代用）
        double tol_rel = 1e-6; ///< 预留：相对阈值（可在 time loop 中使用）

        int    max_outer = 10;    ///< 单时间步最大外迭代次数
        bool   verbose = false; ///< 是否打印每轮外迭代的 dp_inf / linRes
	};

    /**
     * \brief 单次压力组装 + 线性求解 的报告
     */
    struct PressureStepReport
    {
        double lin_residual = 0.0;      ///< 线性求解器返回残差（如果可用）
        int    lin_iterations = 0;      ///< 线性求解器迭代步数
        double dp_inf = 0.0;            ///< 本次 outer 中压力场的无穷范数变化
        int    outer_iterations = 0;   ///< 实际使用的外迭代次数
    };

    //================= 小工具：收集/散射 cell 场 =================//

   /**
    * \brief 把 cell 场收集到未知向量（只填充 cell 部分 DOF）
    *
    * @param reg          字段注册表
    * @param mesh         网格
    * @param name         要收集的 cell 标量场名
    * @param lid_of_cell  cell -> 全局未知编号的映射（来自 PressureAssemblyResult::cell_lid）
    * @param nUnknowns    全局未知个数（包含井 DOF）
    * @param out          输出向量（长度 nUnknowns，先置零）
    */
    inline void gatherFieldToVector(
        const FieldRegistry& reg,
        Mesh& mesh,
        const std::string& name,
        const std::vector<int>& lid_of_cell,
        int nUnknowns,
        std::vector<double>& out)
    {
        out.assign(nUnknowns, 0.0);
        auto fld = reg.get<volScalarField>(name);
        if (!fld) return;

        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();
        for (const auto& c : cells)
        {
            if (c.id < 0) continue;
            const size_t i = id2idx.at(c.id);
            const int    r = lid_of_cell[i];
            if (r >= 0 && r < nUnknowns)
            {
                out[r] = (*fld)[i];
            }
        }
    }
    /**
     * \brief 把未知向量中的 cell 部分散射回指定 cell 场
     *
     * @param reg          字段注册表（目标场会自动 getOrCreate）
     * @param mesh         网格
     * @param name         要写回的 cell 标量场名
     * @param lid_of_cell  cell -> 全局未知编号的映射
     * @param vec          解向量（长度 nUnknowns）
     */
    inline void scatterVectorToField(
        FieldRegistry& reg,
        Mesh& mesh,
        const std::string& name,
        const std::vector<int>& lid_of_cell,
        const std::vector<double>& vec)
    {
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();
        auto fld = reg.getOrCreate<volScalarField>(name, cells.size(), 0.0);

        for (const auto& c : cells)
        {
            if (c.id < 0) continue;
            const size_t i = id2idx.at(c.id);
            const int    r = lid_of_cell[i];
            if (r >= 0 && static_cast<size_t>(r) < vec.size())
            {
                (*fld)[i] = vec[r];
            }
        }
    }

    //================= 单次压力方程求解器（单 outer step） =================//

    /**
     * \brief
     * 单次 “组装两相流总压力方程 + 线性求解器 + 欠松弛 + 井底压力回写”。
     *
     * 假设：
     *  - 时间步 dt 已给定；
     *  - 本时间步的物性（lambda, rho, c_t 等）已经在外面更新好，并存入 reg；
     *  - `assemblePressureTwoPhase` 已经实现完整的时间项 + 扩散 + 毛细 + 重力 + 井耦合，并返回
     *    - `PressureAssemblyResult::system` 线性系统
     *    - `PressureAssemblyResult::cell_lid` cell -> DOF 映射
     *    - 面上的总质量/体积/速度通量（若在 cfg 中配置了名字）。:contentReference[oaicite:2]{index=2}
     *
     * 本函数：
     *  - 只做 “一次 outer 迭代”：不在内部循环；
     *  - 适合在外层时间步/外迭代 driver 中反复调用；
     *  - 通过 PressureSolveReport 返回 dp_inf 和通量信息。
     */
    inline bool solver_IMPES_Iteration_PressureEq(
        MeshManager& mgr,
        FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const PressureBCAdapter& Pbc,
        std::vector<WellDOF_TwoPhase>& wells,
        double                    dt,
        const PressureSolveControls& ctrl,
        // per-outer 输出量：
        double& dp_inf,
        double& lin_residual,
        int& lin_iterations
    )
    {
        if (dt <= 0.0)
        {
            std::cerr << "[IMPES_Iteration][Pressure] invalid dt.\n";
            return false;
        }
        // ===== 1) 组装压力方程（包括时间项、扩散、毛细、重力、井耦合、总通量） ===== //
        PressureAssemblyResult asmb;  //创建一个组装结果的实体
        if (!assemblePressureTwoPhase(mgr, reg, freg, Pbc, wells, dt, ctrl.assembly, asmb))
        {
            std::cerr << "[IMPES_Iteration][Pressure] assemblePressureTwoPhase failed.\n";
            return false;
        }

        const int nUnknowns = asmb.system.n;

        // ===== 2) 用当前压力场作为初始猜测向量 ===== //
        std::vector<double> pvec;
        gatherFieldToVector(
            reg,
            mgr.mesh(),
            ctrl.assembly.pressure_field,
            asmb.cell_lid,
            nUnknowns,
            pvec);
        
        // ===== 3) 调用线性求解器 ===== //
        double linRes = 0.0;
        int    linIters = 0;
        if (!solveCOO_Eigen(asmb.system, pvec, ctrl.linear, &linIters, &linRes))
        {
            std::cerr << "[IMPES_Iteration][Pressure] linear solver failed.\n";
            return false;
        }

        // ===== 4) 把 cell 部分解写回 pressure_field ===== //
        scatterVectorToField(
            reg,
            mgr.mesh(),
            ctrl.assembly.pressure_field,
            asmb.cell_lid,
            pvec);

        // ===== 5) 更新井底压力：p_bh = 解向量中对应 DOF ===== //
        for (auto& w : wells)
        {
            if (w.lid >= 0 && w.lid < nUnknowns)
            {
                w.p_bh = pvec[w.lid];
            }
        }

        // ===== 6) 对 cell 压力做欠松弛（相对于 *_prev） ===== //
        double dpInf = 0.0;
        if (ctrl.under_relax > 0.0 && ctrl.under_relax < 1.0 &&
            !ctrl.assembly.pressure_prev_field.empty())
        {
            // 欠松弛：p_new = urf * p_new + (1-urf) * p_prev
            underRelaxInPlace(
                reg,
                ctrl.assembly.pressure_field,
                ctrl.assembly.pressure_prev_field,
                ctrl.under_relax);
        }

        // ===== 7) 计算 dp_inf，并更新 prev 场 ===== //
        if (!ctrl.assembly.pressure_prev_field.empty())
        {
            dpInf = maxAbsDiff(
                reg,
                ctrl.assembly.pressure_field,
                ctrl.assembly.pressure_prev_field);

            updatePrevIterates(
                reg,
                {
                    { ctrl.assembly.pressure_field,
                      ctrl.assembly.pressure_prev_field }
                });
        }

        // ===== 8) 填充报告 ===== //
        dp_inf = dpInf;
        lin_residual = linRes;
        lin_iterations = linIters;
        // 总质量/体积/速度通量已经写在 freg 里（由 assemblePressureTwoPhase 按字段名存放）
        // 后续两相 flux 分配 / 饱和度方程直接从 FaceFieldRegistry 读取即可，无需在此返回指针。


        return true;
    }

}// namespace IMPES_Iteration

