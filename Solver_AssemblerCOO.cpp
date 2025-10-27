#include "Solver_AssemblerCOO.h"
#include <sstream>
#include <algorithm>
#include <cctype>
#include <memory>
#include <iostream>

// ===== 字符串 -> 位掩码解析(将所有大小写统一转换为同样的格式) =====
static inline std::string normToken(std::string s) {
    // 去除所有空白
    s.erase(std::remove_if(s.begin(), s.end(),
        [](unsigned char ch) { return std::isspace(ch); }),
        s.end());
    // 转小写
    std::transform(s.begin(), s.end(), s.begin(),
        [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    return s;
}


unsigned parse_ops(const std::string& expr, std::string* errMsg)
{
    if (expr.empty()) return OP_NONE;

    std::string s = expr;
    // 常见分隔符统一为 '+'
    for (char& c : s) { if (c == '|' || c == ',' || c == ';') c = '+'; }

    std::istringstream ss(s);
    std::string tok;
    unsigned mask = OP_NONE;

    while (std::getline(ss, tok, '+')) {
        auto t = normToken(tok);
        if (t.empty()) continue;

        if (t == "ddt" || t == "time" || t == "transient") { mask |= OP_TIME;       continue; }
        if (t == "diffusion" || t == "laplacian" || t == "laplace") { mask |= OP_DIFFUSION;  continue; }
        if (t == "convection" || t == "advection" || t == "convective") { mask |= OP_CONVECTION; continue; }

        if (errMsg) *errMsg = "unknown operator token: '" + tok + "'";
        return OP_NONE;
    }

    if (mask == OP_NONE && errMsg)
        *errMsg = "empty or invalid op expression: '" + expr + "'";
    return mask;
}

// 三个算子分别进行装配函数的封装
// 1)扩散项
void assemble_diffusion_faces(const AssembleCtx& Ctx, SparseSystemCOO& out)
{
    for (const auto& F : Ctx.faces) {
        const int iF = F.id - 1;
        const double af = (*Ctx.aF)[iF];
        const double sf = (*Ctx.sF)[iF];

        const int Pid = F.ownerCell, Nid = F.neighborCell;
        const int iP = (Pid >= 0 ? Ctx.id2idx.at(Pid) : -1);
        const int iN = (Nid >= 0 ? Ctx.id2idx.at(Nid) : -1);
        const int rP = (iP >= 0 ? Ctx.lid_of_cell[iP] : -1);

        if (F.isBoundary()) {
            if (rP >= 0) { out.addA(rP, rP, af); out.addb(rP, sf); }
            continue;
        }
        const int rN = (iN >= 0 ? Ctx.lid_of_cell[iN] : -1);
        if (rP >= 0 && rN >= 0) {
            out.addA(rP, rP, af);  out.addA(rP, rN, -af); out.addb(rP, sf);
            out.addA(rN, rN, af);  out.addA(rN, rP, -af); out.addb(rN, -sf);
        }
    }
}

// 2） 对流项
void assemble_convection_faces(const AssembleCtx& Ctx, SparseSystemCOO& out)
{
    for (const auto& F : Ctx.faces) {
        const int iF = F.id - 1;
        const double aPPf = (*Ctx.aPP)[iF];
        const double aPNf = (*Ctx.aPN)[iF];
        const double bPf = (*Ctx.bPc)[iF];

        const int Pid = F.ownerCell, Nid = F.neighborCell;
        const int iP = (Pid >= 0 ? Ctx.id2idx.at(Pid) : -1);
        const int iN = (Nid >= 0 ? Ctx.id2idx.at(Nid) : -1);
        const int rP = (iP >= 0 ? Ctx.lid_of_cell[iP] : -1);

        if (F.isBoundary()) {
            if (rP >= 0) {
                if (aPPf != 0.0) out.addA(rP, rP, aPPf); // 出流：对角增量
                if (bPf != 0.0) out.addb(rP, bPf);      // 入流：RHS 源（内部面应为0）
            }
            continue;
        }
        const int rN = (iN >= 0 ? Ctx.lid_of_cell[iN] : -1);
        if (rP < 0 || rN < 0) continue;

        if (aPPf > 0.0) {              // owner 出流
            out.addA(rP, rP, aPPf);
            out.addA(rN, rP, -aPPf);
        }
        else if (aPNf != 0.0) {      // owner 入流（aPN<=0）
            out.addA(rP, rN, aPNf);
            out.addA(rN, rN, -aPNf);
        }
    }
}

// 3） 时间项
void assemble_time_cells(const AssembleCtx& Ctx, SparseSystemCOO& out)
{
    for (const auto& c : Ctx.cells) {
        if (c.id < 0) continue;
        const size_t i = Ctx.id2idx.at(c.id);
        const int r = Ctx.lid_of_cell[i];
        if (r < 0) continue;
        out.addA(r, r, (*Ctx.aC)[i]);
        out.addb(r, (*Ctx.bC)[i]);
    }
}

//==================== 顶层调度（位掩码版） ====================//

// 预留 nnz 的经验系数（如采用更致密格式可上调）
static constexpr size_t kReservePerFaceDiff = 4; // 若 MPFA/非正交较强可调为 8
static constexpr size_t kReservePerFaceConv = 2; // 若中心差分/对称对流可调为 4
static constexpr size_t kReservePerCellTime = 1;

bool assemble_sparse_system_coo_item
(
    MeshManager& mgr,
    const FieldRegistry& reg,
    const FaceFieldRegistry& freg,
    unsigned ops, 
    const OperatorFieldNames& nm,
    SparseSystemCOO* out
)
{
    Mesh& mesh = const_cast<Mesh&>(mgr.mesh());
    const auto& faces = mesh.getFaces();
    const auto& cells = mesh.getCells();
    const auto& id2idx = mesh.getCellId2Index();

    // —— 取字段（注意：get 返回 shared_ptr，这里先接住，再 .get()） —— //
    std::shared_ptr<faceScalarField> aF_sp, sF_sp, aPP_sp, aPN_sp, bPc_sp;
    std::shared_ptr<volScalarField>  aC_sp, bC_sp;

    if (has(ops, OP_DIFFUSION)) {
        aF_sp = freg.get<faceScalarField>(nm.a_f_diff);
        sF_sp = freg.get<faceScalarField>(nm.s_f_diff);
        if (!aF_sp || !sF_sp) { std::cerr << "[assemble] missing diffusion fields.\n"; return false; }
    }
    if (has(ops, OP_CONVECTION)) {
        aPP_sp = freg.get<faceScalarField>(nm.aPP_conv);
        aPN_sp = freg.get<faceScalarField>(nm.aPN_conv);
        bPc_sp = freg.get<faceScalarField>(nm.bP_conv);
        if (!aPP_sp || !aPN_sp || !bPc_sp) { std::cerr << "[assemble] missing convection fields.\n"; return false; }
    }
    if (has(ops, OP_TIME)) {
        aC_sp = reg.get<volScalarField>(nm.a_time);
        bC_sp = reg.get<volScalarField>(nm.b_time);
        if (!aC_sp || !bC_sp) { std::cerr << "[assemble] missing time fields.\n"; return false; }
    }

    // —— 建未知量编号 —— //
    int N = 0;
    const auto lid_of_cell = buildUnknownMap(mesh, N);
    if (!out) return true;

    // —— 估算 nnz 并 reset —— //
    size_t reserveA = 0;
    if (has(ops, OP_DIFFUSION))  reserveA += faces.size() * kReservePerFaceDiff;
    if (has(ops, OP_CONVECTION)) reserveA += faces.size() * kReservePerFaceConv;
    if (has(ops, OP_TIME))       reserveA += cells.size() * kReservePerCellTime;
    out->reset(N, reserveA);

    // —— 只读上下文（使用原生指针，作用域内 shared_ptr 保活） —— //
    const faceScalarField* aF = aF_sp ? aF_sp.get() : nullptr;
    const faceScalarField* sF = sF_sp ? sF_sp.get() : nullptr;
    const faceScalarField* aPP = aPP_sp ? aPP_sp.get() : nullptr;
    const faceScalarField* aPN = aPN_sp ? aPN_sp.get() : nullptr;
    const faceScalarField* bPc = bPc_sp ? bPc_sp.get() : nullptr;
    const volScalarField* aC = aC_sp ? aC_sp.get() : nullptr;
    const volScalarField* bC = bC_sp ? bC_sp.get() : nullptr;


    AssembleCtx Ctx{
        mesh, faces, cells, id2idx, lid_of_cell,
        aF, sF, aPP, aPN, bPc, aC, bC
    };

    // —— 分支装配 —— //
    if (has(ops, OP_DIFFUSION))  assemble_diffusion_faces(Ctx, *out);
    if (has(ops, OP_CONVECTION)) assemble_convection_faces(Ctx, *out);
    if (has(ops, OP_TIME))       assemble_time_cells(Ctx, *out);

    return true;


}