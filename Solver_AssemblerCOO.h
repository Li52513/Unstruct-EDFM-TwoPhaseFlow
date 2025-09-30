#pragma once
#pragma once
#include <vector>
#include <cassert>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"

//********************COO三元组稀疏矩阵系统组装器********************//
struct Triplet { int r, c; double v; };  //稀疏矩阵三元组 (row, col, value)

struct SparseSystemCOO 
{
    int n = 0;
    std::vector<Triplet> A;
    std::vector<double>  b;

    void reset(int nUnknowns, size_t reserveA = 0)
    {
        n = nUnknowns; A.clear(); b.assign(n, 0.0);
        if (reserveA) A.reserve(reserveA);
    }

    inline void addA(int r, int c, double v) 
    {
        if (v == 0.0) return;
        assert(r >= 0 && r < n && c >= 0 && c < n);
        A.push_back({ r,c,v });
    }
    inline void addb(int r, double v) 
    {
        if (v == 0.0) return;
        assert(r >= 0 && r < n);
        b[r] += v;                 // RHS 已经是 +=
    }

    // —— 新增：COO -> COO（原位）压缩合并（相同 (r,c) 累加）
    void compressInPlace(double drop_tol = 0.0) 
    {
        if (A.empty()) return;
        // 排序 by (r,c)
        std::sort(A.begin(), A.end(), [](const Triplet& a, const Triplet& b)
        {
            return (a.r < b.r) || (a.r == b.r && a.c < b.c);
        });
        // 合并
        size_t write = 0;
        for (size_t k = 0; k < A.size(); ) 
        {
            int r = A[k].r, c = A[k].c;
            double acc = 0.0;
            while ( k < A.size() && A[k].r == r && A[k].c == c ) { acc += A[k].v; ++k; }
            if (std::abs(acc) > drop_tol) { A[write++] = { r,c,acc }; }
        }
        A.resize(write);
    }

    // （可选）导出 CSR，便于自写/第三方求解器
    void toCSR(std::vector<int>& rowPtr, std::vector<int>& colInd, std::vector<double>& val,
        double drop_tol = 0.0)
    {
        compressInPlace(drop_tol);
        rowPtr.assign(n + 1, 0);
        for (auto& t : A) rowPtr[t.r + 1]++;
        for (int i = 0; i < n; i++) rowPtr[i + 1] += rowPtr[i];
        colInd.resize(A.size()); val.resize(A.size());
        auto next = rowPtr;
        for (auto& t : A) { int p = next[t.r]++; colInd[p] = t.c; val[p] = t.v; }
        // 行内已按列排序（因为整体已 sort）
    }

#ifdef USE_ARMADILLO
    // （可选）导入 Armadillo
    arma::sp_mat toArma(double drop_tol = 0.0) {
        compressInPlace(drop_tol);
        arma::umat locs(2, A.size());
        arma::vec  vals(A.size());
        for (size_t k = 0; k < A.size(); ++k) {
            locs(0, k) = A[k].r; locs(1, k) = A[k].c; vals(k) = A[k].v;
        }
        return arma::sp_mat(locs, vals, n, n, /*sort_duplicates=*/false, /*check_for_zeros=*/false);
    }
#endif
};

//===================小工具：建立 “cellId -> 方程自由度编号”===================//
inline std::vector<int> buildUnknownMap(Mesh& mesh, int& nUnknowns)
{
	const auto& cells = mesh.getCells();
	const auto& id2idx = mesh.getCellId2Index();
	std::vector<int> lid_of_cell(cells.size(), -1); // -1 表示该 cellId 不在方程中（如 ghost）
	nUnknowns = 0;
	for (const auto& c : cells) {
		if (c.id < 0) continue; // 忽略 ghost（你目前没有 ghost，但保留更健壮）
		const size_t i = id2idx.at(c.id);
		//cout << "Cell ID: " << c.id << "  index: " << i << endl;
		lid_of_cell[i] = nUnknowns++;
		//cout << "  方程自由度编号: " << lid_of_cell[i] << endl;
	}
	return lid_of_cell;
}

// ======================================================================
//  压力方程（CO2，单相）：装配  A p = b
//   - 面账本：a_f（"a_f_Diff_p_g"），s_f（"s_f_Diff_p_g"）
//   - 时间项：aC_time_p，bC_time_p
// 约定（你已经确认）：
//   内部面：
//     行P: A_PP+=a_f, A_PN+=-a_f, b_P+=+s_f
//     行N: A_NN+=a_f, A_NP+=-a_f, b_N+=-s_f
//   边界面：行P: A_PP+=a_f, b_P+=s_f
//   时间项：cell i: A_ii+=aC[i], b_i+=bC[i]
// ======================================================================

inline bool assemblePressure_singlePhase_COO
(
	MeshManager& mgr,
	const FieldRegistry& reg,
	const FaceFieldRegistry& freg,
	const std::string& a_face ,
	const std::string& s_face ,
	const std::string& a_time ,
	const std::string& b_time ,
	SparseSystemCOO* out = nullptr
)
{
	Mesh& mesh = const_cast<Mesh&>(mgr.mesh());
	const auto& faces = mesh.getFaces();
	const auto& cells = mesh.getCells();
	const auto& id2idx = mesh.getCellId2Index();

	auto aF = freg.get<faceScalarField>(a_face);
	auto sF = freg.get<faceScalarField>(s_face);
	auto aC = reg.get<volScalarField>(a_time);
	auto bC = reg.get<volScalarField>(b_time);
	if (!aF || !sF || !aC || !bC) {
		std::cerr << "[Assemble-Pressure] missing face or cell fields.\n";
		return false;
	}

	//建立未知量编号
	int N = 0;
	const auto lid_of_cell = buildUnknownMap(mesh, N);
	if (!out) return true;
	out->reset(N, /*reserveA*/ faces.size() * 4 + cells.size() * 2); // 预分配  为什么要乘以4 

	// 1) 面贡献（扩散/浮力/交叉）→ 稀疏矩阵 & RHS
	for (const auto& F : faces)
	{
		const int iF = F.id - 1;
		const double af = (*aF)[iF];
		const double sf = (*sF)[iF];
		const int Pid = F.ownerCell;
		const int Nid = F.neighborCell;
		const int iP = (Pid >= 0 ? id2idx.at(Pid) : -1);
		const int iN = (Nid >= 0 ? id2idx.at(Nid) : -1);
		const int rP = (iP >= 0 ? lid_of_cell[iP] : -1);
		if (F.isBoundary()) {
			if (rP >= 0) {
				out->addA(rP, rP, af);
				out->addb(rP, sf);
			}
			continue;
		}
		const int rN = (iN >= 0 ? lid_of_cell[iN] : -1);
		if (rP >= 0 && rN >= 0) 
		{
			// 行P
			out->addA(rP, rP, af);
			out->addA(rP, rN, -af);
			out->addb(rP, sf);
			// 行N
			out->addA(rN, rN, af);
			out->addA(rN, rP, -af);
			out->addb(rN, -sf);
		}

	
	}

	//2 ) 时间项贡献 → 对角 & RHS
	for (const auto& c : cells) 
	{
		if (c.id < 0) continue;
		const size_t i = id2idx.at(c.id);
		const int r = lid_of_cell[i];
		if (r < 0) continue;
		out->addA(r, r, (*aC)[i]);
		out->addb(r, (*bC)[i]);
	}
	return true;

}


// 水相/CO2 两个薄封装（字段名统一在这儿）
inline bool assemblePressure_water_singlePhase_COO
(
    MeshManager& mgr, const FieldRegistry& reg, const FaceFieldRegistry& freg, SparseSystemCOO* out
)
{
    return assemblePressure_singlePhase_COO
    (
        mgr, reg, freg,
        /*a_face*/"a_f_Diff_p_w", /*s_face*/"s_f_Diff_p_w",
        /*a_time*/"aC_time_p",    /*b_time*/"bC_time_p",
        out
    );
}

inline bool assemblePressure_CO2_singlePhase_COO
(
    MeshManager& mgr, const FieldRegistry& reg, const FaceFieldRegistry& freg, SparseSystemCOO* out
)
{
    return assemblePressure_singlePhase_COO
    (
        mgr, reg, freg,
        /*a_face*/"a_f_Diff_p_g", /*s_face*/"s_f_Diff_p_g",
        /*a_time*/"aC_time_p",    /*b_time*/"bC_time_p",
        out
    );
}

// ======================================================================
//  温度方程（单相，基于达西速度）：装配  A T = b
//   - 扩散面账本：a_f_Diff_T / s_f_Diff_T
//   - 对流面账本：aPP_conv / aPN_conv / bP_conv   （你的对流函数已写入这些面场）
//   - 时间项：aC_time_T / bC_time_T
// 双侧守恒装配（你确认同意）：
//   内部面对流：
//     若 aPP_conv>0（owner出流） → 行P: A_PP+=aPP；行N: A_NP+=-aPP
//     若 aPN_conv<0（owner入流） → 行P: A_PN+=aPN；行N: A_NN+=-aPN
//   边界面对流：行P: A_PP+=aPP_conv(F)，b_P+=bP_conv(F)
// ======================================================================
inline bool assembleTemperature_singlePhase_COO
(
    MeshManager& mgr,
    const FieldRegistry& reg,
    const FaceFieldRegistry& freg,
    // 扩散面账本
    const std::string& a_face_diff = "a_f_Diff_T",
    const std::string& s_face_diff = "s_f_Diff_T",
    // 对面账本
    const std::string& aPP_conv = "aPP_conv",
    const std::string& aPN_conv = "aPN_conv",
    const std::string& bP_conv = "bP_conv",
    // 时间项
    const std::string& a_time = "aC_time_T",
    const std::string& b_time = "bC_time_T",
    SparseSystemCOO* out = nullptr
) {
    Mesh& mesh = const_cast<Mesh&>(mgr.mesh());
    const auto& faces = mesh.getFaces();
    const auto& cells = mesh.getCells();
    const auto& id2idx = mesh.getCellId2Index();

    auto aF = freg.get<faceScalarField>(a_face_diff);
    auto sF = freg.get<faceScalarField>(s_face_diff);
    auto aPP = freg.get<faceScalarField>(aPP_conv);
    auto aPN = freg.get<faceScalarField>(aPN_conv);
    auto bPc = freg.get<faceScalarField>(bP_conv);
    auto aC = reg.get<volScalarField>(a_time);
    auto bC = reg.get<volScalarField>(b_time);

    if (!aF || !sF || !aPP || !aPN || !bPc || !aC || !bC) {
        std::cerr << "[assembleTemperature] missing diffusion/conv/time fields.\n";
        return false;
    }

    int N = 0;
    const auto lid_of_cell = buildUnknownMap(mesh, N);
    if (!out) return true;
    out->reset(N, /*reserveA*/ faces.size() * 6 + cells.size() * 2);

    // 1) 扩散贡献（与压力相同的装配法则）
    for (const auto& F : faces) {
        const int iF = F.id - 1;
        const double af = (*aF)[iF];
        const double sf = (*sF)[iF];
        const int Pid = F.ownerCell;
        const int Nid = F.neighborCell;
        const int iP = (Pid >= 0 ? id2idx.at(Pid) : -1);
        const int iN = (Nid >= 0 ? id2idx.at(Nid) : -1);
        const int rP = (iP >= 0 ? lid_of_cell[iP] : -1);

        if (F.isBoundary()) {
            if (rP >= 0) { out->addA(rP, rP, af); out->addb(rP, sf); }
            continue;
        }
        const int rN = (iN >= 0 ? lid_of_cell[iN] : -1);
        if (rP >= 0 && rN >= 0) {
            out->addA(rP, rP, af);
            out->addA(rP, rN, -af);
            out->addb(rP, sf);
            out->addA(rN, rN, af);
            out->addA(rN, rP, -af);
            out->addb(rN, -sf);
        }
    }

    // 2) 对流贡献（双侧守恒）
    for (const auto& F : faces) {
        const int iF = F.id - 1;
        const int Pid = F.ownerCell;
        const int Nid = F.neighborCell;
        const int iP = (Pid >= 0 ? id2idx.at(Pid) : -1);
        const int iN = (Nid >= 0 ? id2idx.at(Nid) : -1);
        const int rP = (iP >= 0 ? lid_of_cell[iP] : -1);

        const double aPPf = (*aPP)[iF];
        const double aPNf = (*aPN)[iF];
        const double bPf = (*bPc)[iF];

        if (F.isBoundary()) {
            if (rP >= 0) {
                if (aPPf != 0.0) out->addA(rP, rP, aPPf); // 出流
                if (bPf != 0.0) out->addb(rP, bPf);   // 入流
            }
            continue;
        }
        const int rN = (iN >= 0 ? lid_of_cell[iN] : -1);
        if (rP < 0 || rN < 0) continue;

        if (aPPf > 0.0) {
            // owner 出流
            out->addA(rP, rP, aPPf);
            out->addA(rN, rP, -aPPf);
        }
        else if (aPNf != 0.0) {
            // owner 入流（aPN <= 0）
            out->addA(rP, rN, aPNf);
            out->addA(rN, rN, -aPNf);
        }
        // 内部面对流没有 b 源（你已在内部面把 bP_conv 置 0）
    }

    // 3) 时间项
    for (const auto& c : cells) {
        if (c.id < 0) continue;
        const size_t i = id2idx.at(c.id);
        const int r = lid_of_cell[i];
        if (r < 0) continue;
        out->addA(r, r, (*aC)[i]);
        out->addb(r, (*bC)[i]);
    }
    return true;
}