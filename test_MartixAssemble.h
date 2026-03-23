#pragma once
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

// 引入你的真实底层头文件
#include "ADVar.hpp"
#include "FIM_BlockSparseMatrix.h"
#include "FIM_GlobalAssembler.h"

namespace MatrixAssemblerTest {
	void Test_BlockSparseMatrix() {

        // 设置自由度 (例如：单相两组分，或者压力+饱和度)
        const int N = 2;
        // 设置网格数量 (测试 3 个网格：0, 1, 2)
        const int total_blocks = 3;

        // 初始化全局块状稀疏矩阵
        FIM_BlockSparseMatrix<N> global_mat(total_blocks);

        std::cout << "=== 阶段 1: 组装累积项 (Accumulation) ===" << std::endl;
        // 假设网格 0 的累积项：
        // 方程 0: val = 10.0, d/dX0 = 1.0, d/dX1 = 0.1
        // 方程 1: val = 20.0, d/dX0 = 0.2, d/dX1 = 2.0
        std::vector<ADVar<N>> acc0(N);
        acc0[0].val = 10.0; acc0[0].grad(0) = 1.0; acc0[0].grad(1) = 0.1;
        acc0[1].val = 20.0; acc0[1].grad(0) = 0.2; acc0[1].grad(1) = 2.0;

        FIM_GlobalAssembler<N, ADVar<N>>::AssembleAccumulation(0, acc0, global_mat);
        std::cout << "网格 0 累积项装配完成。" << std::endl;

        std::cout << "\n=== 阶段 2: 组装跨面通量 (Flux) ===" << std::endl;
        // 模拟网格 0 到网格 1 的通量 (物理约定：算出的正值为流入网格 0 的量)
        std::vector<ADVar<N>> f_wrt_0(N), f_wrt_1(N);

        // 方程 0 (Eq 0) 跨面通量
        f_wrt_0[0].val = 4.0;
        f_wrt_0[0].grad(0) = 0.5;   // 对 网格 0 主元的导数
        f_wrt_0[0].grad(1) = 0.05;

        // 注意：AD 对网格 1 的求导，绝对值通常相等，导数符号通常相反
        f_wrt_1[0].val = 4.0;
        f_wrt_1[0].grad(0) = -0.5;  // 对 网格 1 主元的导数
        f_wrt_1[0].grad(1) = -0.05;

        // 方程 1 (Eq 1) 假设当前界面没有流体流动，保持默认 0.0
        // (ADVar 的默认构造函数已经将 val 和 grad 设为 0)

        // 调用装配器，装配网格 0 和 1 的连接
        FIM_GlobalAssembler<N, ADVar<N>>::AssembleFlux(0, 1, f_wrt_0, f_wrt_1, global_mat);
        std::cout << "网格 0 和 1 的通量装配完成 (验证 'i减j加' 与四角求导)。" << std::endl;

        std::cout << "\n=== 阶段 3: 组装源汇项 (Source/Well) ===" << std::endl;
        // 假设在网格 1 有一口生产井 (物理约定：流出系统为正)
        // 产水量 3.0，对网格 1 第 0 个主元(如压力)的偏导为 0.3
        std::vector<ADVar<N>> src1(N);
        src1[0].val = 2.0;
        src1[0].grad(0) = 0.3; src1[0].grad(1) = 0.0;

        FIM_GlobalAssembler<N, ADVar<N>>::AssembleSource(1, src1, global_mat);
        std::cout << "网格 1 生产井装配完成。" << std::endl;

        std::cout << "\n=== 阶段 4: 拓扑冻结与矩阵导出 ===" << std::endl;
        // 锁定非零结构，构建底层 CSR 直写缓存
        global_mat.FreezePattern();

        // 获取导出的 Eigen 稀疏矩阵和右端项 (Residual)
        Eigen::SparseMatrix<double, Eigen::RowMajor, int> Jacobian = global_mat.GetFrozenMatrix();
        Eigen::VectorXd Residual = global_mat.ExportEigenResidual();

        std::cout << "\n[最终残差向量 (Residual)]:" << std::endl;
        std::cout << "Block 0 (Eq0, Eq1): " << Residual(0) << ", " << Residual(1) << std::endl;
        std::cout << "Block 1 (Eq0, Eq1): " << Residual(2) << ", " << Residual(3) << std::endl;
        std::cout << "Block 2 (Eq0, Eq1): " << Residual(4) << ", " << Residual(5) << std::endl;

        std::cout << "\n[最终雅可比稀疏矩阵 (Jacobian)]:" << std::endl;
        // 使用 Eigen 自带的稠密打印，直接查看 6x6 全矩阵结构
        std::cout << Eigen::MatrixXd(Jacobian) << std::endl;

	}

}