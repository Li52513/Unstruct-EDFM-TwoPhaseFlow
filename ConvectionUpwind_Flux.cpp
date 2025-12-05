#include "ConvectionUpwind_Flux.h"
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "FieldAcessForDiscre.h"  // cellScalar(...)
#include <algorithm>
#include <iostream>
#include "BCAdapter.h"
#include "TemperatureBCAdapter.h"

namespace FVM 
{ 
	namespace Convection
	{
		inline bool isAdiabaticZeroFlux(const TemperatureBCAdapter* bc, int fid)
		{
			if (!bc) return false;
			double a = 0.0, b = 0.0, c = 0.0;
			if (!bc->getABC(fid, a, b, c)) return false;
			const double eps = 1e-30;
			return (std::abs(a) <= eps && std::abs(b) > eps && std::abs(c) <= eps);
		}


		inline bool isPureNeumann(const PressureBCAdapter* bc, int faceId) {
			if (!bc) return false;
			double a = 0.0, b = 0.0, c = 0.0;
			if (!bc->getABC(faceId, a, b, c)) return false;
			const double eps = 1e-30;
			return (std::abs(a) <= eps && std::abs(b) > eps);
		}
		inline void faceAreaVec(const Face& F, Vector& A, double& Aabs) {
		A = F.vectorE + F.vectorT;
		Aabs = A.Mag();
		if (Aabs <= 0.0) { A = F.normal; Aabs = 0.0; } // 兜底
		}

		inline bool ensure_constant_U_field
		(
			MeshManager& mgr, FieldRegistry& reg,
			const std::string& U_name, const Vector& U_const)
		{
			auto Uv = reg.get<volVectorField>(U_name);
			if (Uv) return true;  // 已有就不动
			const size_t n = mgr.mesh().getCells().size();
			auto Nu = reg.getOrCreate<volVectorField>(U_name.c_str(), n, Vector{ 0,0,0 });
			if (!Nu) return false;
			std::fill(Nu->data.begin(), Nu->data.end(), U_const);
			return true;
		}

		// 路线 A1：达西（质量通量形态，复用 a_f/s_f） //复用a_f_Diff_T s_f_Diff_T 即包含密度的达西定律离散系数与面向量
		// 约定：mf = a_f * Δp − s_f（内部 Δp=pP−pN；边界 Δp=pP）
		/**
		 * @brief 根据 a_f / s_f 面系数构建达西质量通量、体积通量和法向速度。
		 *
		 * - 内部面：使用 m_f = a_f (p_P - p_N) - s_f
		 * - Dirichlet / Robin 边界：使用 m_f = a_f p_P - s_f（p_B 已编码在 s_f 中）
		 * - 纯 Neumann 零梯度边界 (a≈0, b≠0, c≈0)：在此函数中解释为“无流边界”，
		 *   直接强制 mf = 0, Qf = 0, ufn = 0，即便 s_f 中有重力/毛细项也一并屏蔽。
		 */
		bool buildFlux_Darcy_Mass
		(
			MeshManager& mgr, const FieldRegistry& reg, FaceFieldRegistry& freg,
			const std::string& a_name, const std::string& s_name,
			const std::string& p_name, const std::string& rho_name,
			const std::string& mf_name, const std::string& Qf_name,
			const std::string& ufn_name,
			const PressureBCAdapter* bc,
			bool   clampDirichletBackflow,
			double dirichlet_zero_flux_tol
		)
		{
			Mesh& mesh = const_cast<Mesh&>(mgr.mesh());
			auto& faces = const_cast<std::vector<Face>&>(mesh.getFaces());
			const auto& id2idx = mesh.getCellId2Index();
			const auto& cells = mesh.getCells();

			auto a_f = freg.get<faceScalarField>(a_name.c_str());
			auto s_f = freg.get<faceScalarField>(s_name.c_str());

			if (!a_f || !s_f) {
				std::cerr << "[buildFlux_Darcy_Mass] face fields '" << a_name
					<< "' or '" << s_name << "' not found.\n";
				return false;
			}

			auto mf = freg.getOrCreate<faceScalarField>(mf_name.c_str(), faces.size(), 0.0);
			auto Qf = freg.getOrCreate<faceScalarField>(Qf_name.c_str(), faces.size(), 0.0);
			auto ufn = freg.getOrCreate<faceScalarField>(ufn_name.c_str(), faces.size(), 0.0);

			std::fill(mf->data.begin(), mf->data.end(), 0.0);
			std::fill(Qf->data.begin(), Qf->data.end(), 0.0);
			std::fill(ufn->data.begin(), ufn->data.end(), 0.0);

			const double eps = 1e-30;

			for (const auto& F : faces) {
				const int idx = F.id - 1;
				Vector A; double Aabs; faceAreaVec(F, A, Aabs);

				const int    P = F.ownerCell;
				const double pP = cellScalar(reg, mesh, p_name.c_str(), P, 0.0);
				const double rhoP = cellScalar(reg, mesh, rho_name.c_str(), P, 1.0);

				double mface = 0.0;
				double qface = 0.0;

				if (!F.isBoundary())
				{
					// -------------------------
					// 内部面：标准达西通量
					// m_f = a_f (pP - pN) - s_f
					// -------------------------
					const int    N = F.neighborCell;
					const double pN = cellScalar(reg, mesh, p_name.c_str(), N, 0.0);
					const double rhoN = cellScalar(reg, mesh, rho_name.c_str(), N, 1.0);

					mface = (*a_f)[idx] * (pP - pN) - (*s_f)[idx];

					const double gamma = F.f_linearInterpolationCoef;
					const double rho_f = (1.0 - gamma) * rhoP + gamma * rhoN;
					const double rho_safe = std::max(rho_f, eps);
					qface = mface / rho_safe;
				}
				else
				{
					// =======================================================
					// 边界面：先识别“无流边界”，再处理 Dirichlet / Robin。
					// 约定：a≈0, b≠0, c≈0 在本函数中表示 No-flow 边界，
					//       即总达西质量通量 mf = 0 （包括毛细 + 重力）。
					// =======================================================
					bool noFlowBoundary = false;
					bool dirichletOutflow = false;

					if (bc) {
						double a_bc = 0.0, b_bc = 0.0, c_bc = 0.0;
						if (bc->getABC(F.id, a_bc, b_bc, c_bc)) {
							const bool isDirichlet = (std::abs(a_bc) > eps && std::abs(b_bc) <= eps);
							const bool isZeroGrad = (std::abs(a_bc) <= eps && std::abs(b_bc) > eps);

							// —— 无流边界判据：纯 Neumann 且 c≈0 —— //
							if (isZeroGrad && std::abs(c_bc) <= dirichlet_zero_flux_tol) {
								noFlowBoundary = true;
							}

							// —— Dirichlet 边界：用于回流截断判断 —— //
							if (isDirichlet && clampDirichletBackflow) {
								const double pBC = c_bc / a_bc;
								// 若单元内压力 >= 边界压力，则“倾向于出流边界”
								dirichletOutflow = (pP >= pBC);
							}
						}
					}

					if (noFlowBoundary)
					{
						// -----------------------------------------------------
						// 真·无流边界：直接强制 mf=0，qf=0，忽略 a_f / s_f 中
						// 预先计算的毛细项和重力项，确保总达西通量为 0。
						// -----------------------------------------------------
						mface = 0.0;
						qface = 0.0;
					}
					else
					{
						// -----------------------------------------------------
						// 其余边界（Dirichlet / Robin）：沿用
						// m_f = a_f pP - s_f
						// 其中 s_f 已经包含 p_B / 重力 / 毛细信息。
						// -----------------------------------------------------
						mface = (*a_f)[idx] * (pP)-(*s_f)[idx];
						qface = (rhoP > eps) ? (mface / rhoP) : 0.0;

						// ---- 对“出流型 Dirichlet 边界”可选截断回流 ---- //
						if (dirichletOutflow && mface < 0.0) {
							mface = 0.0;
							qface = 0.0;
						}

						// ---- 对纯 Dirichlet 边界，当 pP≈pB 时，可选零通量 ---- //
						if (dirichlet_zero_flux_tol > 0.0 && bc) {
							double a_bc = 0.0, b_bc = 0.0, c_bc = 0.0;
							if (bc->getABC(F.id, a_bc, b_bc, c_bc)
								&& std::abs(a_bc) > eps && std::abs(b_bc) <= eps)
							{
								const double pBC = c_bc / a_bc;
								if (std::abs(pP - pBC) <= dirichlet_zero_flux_tol) {
									mface = 0.0;
									qface = 0.0;
								}
							}
						}
					}
				}

				(*mf)[idx] = mface;
				(*Qf)[idx] = qface;
				(*ufn)[idx] = (Aabs > 0.0 ? qface / Aabs : 0.0);
			}

			return true;
		}

		bool buildFlux_Darcy_Vol  //通过FVM::Diffusion::buildFaceCoeff_Cenrtal计算得到的a_f和s_f 进而计算质量通量 面通量 和 达西速度
		(
			MeshManager& mgr, const FieldRegistry& reg, FaceFieldRegistry& freg,
			const std::string& a_vol,    // 输入：a_f_diffusion
			const std::string& s_vol,    // 输入：s_f_diffusion
			const std::string& p_name,   // 输入：压力场 
			const std::string& rho_name, // 输出：密度体场名（用于 ρ_up）
			const std::string& mf_name,  // 输出：质量通量面场
			const std::string& Qf_name,  // 输出：体积通量面场
			const std::string& ufn_name,
			const PressureBCAdapter* Pbc,
			const TemperatureBCAdapter* tBC)
		{
			Mesh& mesh = const_cast<Mesh&>(mgr.mesh());
			auto& faces = const_cast<std::vector<Face>&>(mesh.getFaces());
			const auto& id2idx = mesh.getCellId2Index();

			// 输入面场：体积形 a_vol / s_vol
			auto aFv = freg.get<faceScalarField>(a_vol);
			auto sFv = freg.get<faceScalarField>(s_vol);
			// 输出面场
			auto mf = freg.getOrCreate<faceScalarField>(mf_name.c_str(), faces.size(), 0.0); //质量通量面场
			auto Qf = freg.getOrCreate<faceScalarField>(Qf_name.c_str(), faces.size(), 0.0); //体积通量面场
			auto ufn = freg.getOrCreate<faceScalarField>(ufn_name.c_str(), faces.size(), 0.0); //达西速度面场

			std::fill(mf->data.begin(), mf->data.end(), 0.0);
			std::fill(Qf->data.begin(), Qf->data.end(), 0.0);
			std::fill(ufn->data.begin(), ufn->data.end(), 0.0);

			const double eps_rho = 1e-30;
			const double eps_A = 1e-30;

			for (const auto& F : faces) {
				const int idx = F.id - 1;
				// 面向量：A = E + T（边界 T=0）
				Vector A; double Aabs; faceAreaVec(F, A, Aabs);

				const int    P = F.ownerCell;
				const size_t iP = id2idx.at(P);

				const double pP = cellScalar(reg, mesh, p_name.c_str(), P, 0.0);
				const double rhoP = cellScalar(reg, mesh, rho_name.c_str(), P, 0.0);

				double mface = 0.0;
				double qface = 0.0;

				if (!F.isBoundary())
				{
					const int    N = F.neighborCell;
					const double pN = cellScalar(reg, mesh, p_name.c_str(), N, 0.0);
					const double rhoN = cellScalar(reg, mesh, rho_name.c_str(), N, 0.0);

					// —— RHS 记账 → 通量重建 “aΔp − s” —— //
					qface = (*aFv)[idx] * (pP - pN) - (*sFv)[idx];

					// 上风密度（以体积通量方向判定）
					const double rho_up = (qface >= 0.0) ? rhoP : rhoN;
					mface = rho_up * qface;
				}
				else {
					if (isPureNeumann(Pbc, F.id)) 
					{
						qface = 0.0;
						mface = 0.0;
					}
					else {
						if (isPureNeumann(Pbc, F.id) || isAdiabaticZeroFlux(tBC, F.id)) {
							qface = 0.0;
							mface = 0.0;
						}
						else {
							qface = (*aFv)[idx] * pP - (*sFv)[idx];
							mface = rhoP * qface;
						}
					}
				}
				(*Qf)[idx] = qface;
				(*mf)[idx] = mface;
				(*ufn)[idx] = (Aabs > eps_A) ? (qface / Aabs) : 0.0;
			}

			return true;
		}

		/**
		 * 路线 B：外给速度（常量或体场 U），直接构造通量
		 * 逻辑：U_f 线性插值到面 → u_n = U_f·n → Qf = u_n * |A| → mf = ρ_up * Qf。
		 * 若 U_name 为空或未找到对应体场，则使用常量 U_const。
		 */

		bool buildFlux_FromVelocity(
			MeshManager& mgr, const FieldRegistry& reg, FaceFieldRegistry& freg,
			const std::string& U_name,
			Vector U_const,
			const std::string& rho_name,
			const std::string& mf_name,
			const std::string& Qf_name,
			const std::string& ufn_name)
		{
			Mesh& mesh = const_cast<Mesh&>(mgr.mesh());
			auto& faces = const_cast<std::vector<Face>&>(mesh.getFaces());
			const auto& id2idx = mesh.getCellId2Index();

			auto rho = reg.get<volScalarField>(rho_name);

			auto mf = freg.getOrCreate<faceScalarField>(mf_name.c_str(), faces.size(), 0.0);
			auto Qf = freg.getOrCreate<faceScalarField>(Qf_name.c_str(), faces.size(), 0.0);
			auto ufn = freg.getOrCreate<faceScalarField>(ufn_name.c_str(), faces.size(), 0.0);
			std::fill(mf->data.begin(), mf->data.end(), 0.0);
			std::fill(Qf->data.begin(), Qf->data.end(), 0.0);
			std::fill(ufn->data.begin(), ufn->data.end(), 0.0);

			const double epsA = 1e-30;

			for (const auto& F : faces) {
				const int idx = F.id - 1;

				// 面几何：A = E + T（边界 T=0），n = A/|A|
				Vector A; double Aabs; faceAreaVec(F, A, Aabs);
				const Vector n = (Aabs > 0.0 ? A / Aabs : F.normal);


				// owner 索引与量
				const int    P = F.ownerCell;
				const double rhoP = cellScalar(reg, mesh, rho_name.c_str(), P, 0.0);

				Vector	Uf{ 0.0,0.0,0.0 }; // 面速度
				// 面速度插值
				if (!F.isBoundary()) {
					// 内部面：按 gamma 做线性插值（若 gamma 无效则按投影求）
					const int    N = F.neighborCell;
					const size_t iN = id2idx.at(N);

					double gamma = F.f_linearInterpolationCoef;

					const Vector UP = cellVector(reg, mesh, U_name.c_str(), P, { 0.0,0.0,0.0 });
					const Vector UN = cellVector(reg, mesh, U_name.c_str(), N, { 0.0,0.0,0.0 });
					Uf = (1.0 - gamma) * UP + gamma * UN;
				}
				else
				{
					Uf = cellVector(reg, mesh, U_name.c_str(), P, { 0.0,0.0,0.0 });
				}
				// 法向速度与体/质量通量
				const double un = Uf * n;         // u_n = U_f · n
				const double qf = un * Aabs;      // Qf = u_n * |A|
				double mface = 0.0;               // mf = ρ_up * Qf

				if (!F.isBoundary()) {
					const int    N = F.neighborCell;
					const double rhoN = cellScalar(reg, mesh, rho_name.c_str(), N, 0.0);

					const double rho_up = (qf >= 0.0) ? rhoP : rhoN;
					mface = rho_up * qf;
				}
				else {
					// 边界：owner→外为正，上风退化为 P
					mface = rhoP * qf;
				}

				// 写回
				(*Qf)[idx] = qf;
				(*mf)[idx] = mface;
				(*ufn)[idx] = (Aabs > epsA ? qf / Aabs : 0.0); // ufn = Qf / |A|

			}
			return true;
		}

	}

}
