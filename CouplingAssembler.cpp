#include "CouplingAssembler.h"
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>   
#include <iostream>
#include <cstdio>
#include <sstream>
#include "FracIndex.h"         // buildFracElemIndex(...)
#include "FractureNetwork.h"   // FractureNetwork / globalFFPts
#include "Fracture.h"          // Fracture / FractureElement::FFF_Exchange
#include "MeshManager.h"       // MeshManager / Mesh
#include "FieldRegistry.h"     // FieldRegistry
#include "VolField.h"          // volScalarField
#include "UserDefineVarType.h" // Vector
#include <iomanip>
// 你已有的物性与相对渗透率工具
// pc_vG(...), kr_Mualem_vG(...)

using std::size_t;
using std::vector;
using volScalarField = VolField<double>;


// 按 id 在 globalFFPts 里查交点；找不到则返回 nullptr
static const FractureNetwork::GlobalFFPoint* findGlobalFFByID(const FractureNetwork& net, int id) {
    for (const auto& P : net.globalFFPts) {
        if (P.id == id) return &P;
    }
    return nullptr;
}

// ------------- 工具：点到线段投影、单位向量 -------------
static Vector projectPointSegment(const Vector& p, const Vector& a, const Vector& b) {
    Vector ab = b - a;
    double L2 = ab.m_x * ab.m_x + ab.m_y * ab.m_y + ab.m_z * ab.m_z;
    if (L2 <= 1e-30) return a;
    double t = ((p - a) * ab) / L2;
    if (t < 0.0) t = 0.0; else if (t > 1.0) t = 1.0;
    return a + ab * t;
}
static Vector unit(const Vector& v) {
    double m = v.Mag();
    return (m > 1e-30) ? (v * (1.0 / m)) : Vector{ 0,0,0 };
}

// ------------- 工具：相位势（可选重力） -------------
static inline double phi(double p, double rho, double z, bool use_g, double g) {
    return use_g ? (p - rho * g * z) : p;
}

// ------------- 工具：相对渗透率（裂缝侧临时计算） -------------
static inline void kr_vG(const VGParams& vg, const RelPermParams& rp, double Sw,
    double& krw, double& krg) {
    kr_Mualem_vG(Sw, vg, rp, krw, krg);
}

// ------------- 工具：矩阵法向投影 k_n = n^T K n -------------
static double kn_project(const FieldRegistry& R, size_t iCell, const Vector& n) {
    auto kxx = R.get<volScalarField>("kxx");
    auto kyy = R.get<volScalarField>("kyy");
    double kxy_val = 0.0;
    if (R.has("kxy")) {
        auto kxy = R.get<volScalarField>("kxy");
        kxy_val = (*kxy)[iCell];
    }
    if (!kxx || !kyy) return 0.0;
    double nx = n.m_x, ny = n.m_y;
    return nx * nx * (*kxx)[iCell] + 2.0 * nx * ny * kxy_val + ny * ny * (*kyy)[iCell];
}

// ------------- 工具：裂缝切向渗透率 k_t（支持 fr_k_t；否则回退 b²/12） -------------
static double kt_project(const FieldRegistry& Rf, size_t gSeg, const Vector& /*t*/) {
    if (Rf.has("fr_k_t")) {
        auto kt = Rf.get<volScalarField>("fr_k_t");
        return (*kt)[gSeg];
    }
    if (Rf.has("fr_aperture")) {
        auto b = Rf.get<volScalarField>("fr_aperture");
        double ap = (*b)[gSeg];
        return (ap * ap) / 12.0; // 平行板立方定律（单位厚度 2D）
    }
    return 0.0;
}


//用于外部检查和诊断
static inline void ensureFracCIDiagFields(FieldRegistry& reg_fr, size_t nSeg)
{
    if (!reg_fr.has("diag_knM"))      reg_fr.create<volScalarField>("diag_knM", nSeg, 0.0);
    if (!reg_fr.has("diag_lamW_up"))  reg_fr.create<volScalarField>("diag_lamW_up", nSeg, 0.0);
    if (!reg_fr.has("diag_lamG_up"))  reg_fr.create<volScalarField>("diag_lamG_up", nSeg, 0.0);
    if (!reg_fr.has("alpha_geom"))    reg_fr.create<volScalarField>("alpha_geom", nSeg, 0.0); // 方便导出
}

static inline void ensureFracTIDiagFields(FieldRegistry& reg_fr, size_t nSeg) {
    if (!reg_fr.has("diag_Gw_arm"))   reg_fr.create<volScalarField>("diag_Gw_arm", nSeg, 0.0);
    if (!reg_fr.has("diag_Gg_arm"))   reg_fr.create<volScalarField>("diag_Gg_arm", nSeg, 0.0);
    if (!reg_fr.has("diag_ktb"))      reg_fr.create<volScalarField>("diag_ktb", nSeg, 0.0);
    if (!reg_fr.has("diag_kt"))      reg_fr.create<volScalarField>("diag_kt", nSeg, 0.0);
}







// ------------- 创建/确保 CI 场 -------------
void ensureFracCouplingFields(FieldRegistry& reg_fr, size_t nSeg) {
    if (!reg_fr.has("CIw")) reg_fr.create<volScalarField>("CIw", nSeg, 0.0);
    if (!reg_fr.has("CIg")) reg_fr.create<volScalarField>("CIg", nSeg, 0.0);
}

// ------------- 基岩-裂缝耦合 CI（两相） -------------
void updateMatrixFractureCI(MeshManager& mgr,
    const FieldRegistry& Rm, FieldRegistry& Rf,
    const VGParams& vg,
    bool upwind,
    bool include_gravity, double g)
{
    auto& mesh = mgr.mesh();
    auto& cells = mesh.getCells();
    auto& net = mgr.fracture_network();

    // 全局段索引
    auto idx = buildFracElemIndex(net);
    // 确保输出场存在
    ensureFracCouplingFields(Rf, idx.total);


    //确保诊断量存在
    ensureFracCIDiagFields (Rf, idx.total);

    
    // —— 取矩阵侧场
    auto Sw_m = Rm.get<volScalarField>("S_w");
    auto pw_m = Rm.get<volScalarField>("p_w");
    auto pg_m = Rm.get<volScalarField>("p_g");
    auto muw_m = Rm.get<volScalarField>("mu_w");
    auto mug_m = Rm.get<volScalarField>("mu_g");
    auto rhow_m = Rm.get<volScalarField>("rho_w");
    auto rhog_m = Rm.get<volScalarField>("rho_g");

    auto kxx_m = Rm.get<volScalarField>("kxx");
    auto kyy_m = Rm.get<volScalarField>("kyy");
    auto kxy_m = Rm.has("kxy") ? Rm.get<volScalarField>("kxy") : nullptr;

    // —— 取裂缝侧场
    auto Sw_f = Rf.get<volScalarField>("Sf_w");
    auto pfw_f = Rf.get<volScalarField>("pf_w");
    auto Tf_f = Rf.get<volScalarField>("Tf");
    auto muw_f = Rf.get<volScalarField>("mu_w");
    auto mug_f = Rf.get<volScalarField>("mu_g");
    auto rhow_f = Rf.get<volScalarField>("rho_w");
    auto rhog_f = Rf.get<volScalarField>("rho_g");

    auto CIw = Rf.get<volScalarField>("CIw");
    auto CIg = Rf.get<volScalarField>("CIg");


    // --取诊断量
    auto knM = Rf.get<volScalarField>("diag_knM");
    auto lamWu = Rf.get<volScalarField>("diag_lamW_up");
    auto lamGu = Rf.get<volScalarField>("diag_lamG_up");
    auto aGeom = Rf.get<volScalarField>("alpha_geom");


    if (!Sw_m || !pw_m || !pg_m || !muw_m || !mug_m || !rhow_m || !rhog_m ||
        !kxx_m || !kyy_m || !Sw_f || !pfw_f || !Tf_f || !CIw || !CIg) {
        std::cerr << "[Coupling] CI: missing fields; update skipped.\n";
        return;
    }

    // —— 预计算裂缝侧 kr（若你未在 reg_fr 中存 kr_w/kr_g，则临时算一遍）
    RelPermParams rp_def; // 用你的默认参数；如需不同区域可传参或查场
    vector<double> krw_f_loc(idx.total), krg_f_loc(idx.total);
    for (size_t gseg = 0; gseg < idx.total; ++gseg) {
        double kw, kg;
        kr_vG(vg, rp_def, (*Sw_f)[gseg], kw, kg);
        krw_f_loc[gseg] = kw;
        krg_f_loc[gseg] = kg;
    }

    // —— 主循环：每段裂缝
    for (size_t f = 0; f < net.fractures.size(); ++f) {
        const auto& F = net.fractures[f];
        const auto& I = F.intersections;
        size_t base = idx.offset[f];

        for (size_t e = 0; e < F.elements.size(); ++e) {
            const auto& E = F.elements[e];
            size_t gseg = base + e;

            // 找宿主 cell
            auto it = mesh.getCellId2Index().find(E.cellID);
            if (it == mesh.getCellId2Index().end()) {
                (*CIw)[gseg] = (*CIg)[gseg] = 0.0;
                continue;
            }
            size_t im = it->second;

            // 方向：cell center -> 段最近点 的单位法向
            const auto& p1 = I[E.id - 1].point;
            const auto& p2 = I[E.id].point;
            Vector q = projectPointSegment(cells[im].center, p1, p2);
            Vector n = unit(q - cells[im].center);

            // 几何：CI_geom 由你在 mgr.ComputeFractureGeometryCouplingCoefficient() 中已填
            double geomCI = E.geomCI;

            // 法向等效渗透率（只用矩阵侧控制；若要调和平均，可在此扩展）
            double knm = kn_project(Rm, im, n);
            (*knM)[gseg] = knm; //存进场

            // 相位势（上风判据）
            double zM = cells[im].center.m_z, zF = 0.0; // 2D
            double pcF = pc_vG((*Sw_f)[gseg], vg);
            // 水相
            double phiMw = phi((*pw_m)[im], (rhow_m ? (*rhow_m)[im] : 0.0), zM, include_gravity, g);
            double phiFw = phi((*pfw_f)[gseg], (rhow_f ? (*rhow_f)[gseg] : 0.0), zF, include_gravity, g);
            bool M_up_w = (phiMw >= phiFw);
            // 气相
            double pgF = (*pfw_f)[gseg] + pcF;
            double phiMg = phi((*pg_m)[im], (rhog_m ? (*rhog_m)[im] : 0.0), zM, include_gravity, g);
            double phiFg = phi(pgF, (rhog_f ? (*rhog_f)[gseg] : 0.0), zF, include_gravity, g);
            bool M_up_g = (phiMg >= phiFg);

            // 可动度（迎风 or 中点）
            double lamw_up = 0.0, lamg_up = 0.0;
            if (!upwind) {
                lamw_up = std::max(0.0, Rm.get<volScalarField>("kr_w")->operator[](im)) /
                    std::max(1e-30, muw_m->operator[](im));
                (*lamWu)[gseg] = lamw_up;

                lamg_up = std::max(0.0, Rm.get<volScalarField>("kr_g")->operator[](im)) /
                    std::max(1e-30, mug_m->operator[](im));
                (*lamGu)[gseg] = lamg_up;


            }
            else {
                if (M_up_w) {
                    lamw_up = std::max(0.0, Rm.get<volScalarField>("kr_w")->operator[](im)) /
                        std::max(1e-30, muw_m->operator[](im));
                    (*lamWu)[gseg] = lamw_up;
                }
                else {
                    lamw_up = std::max(0.0, krw_f_loc[gseg]) /
                        std::max(1e-30, (muw_f ? (*muw_f)[gseg] : muw_m->operator[](im)));
                    (*lamWu)[gseg] = lamw_up;
                }
                if (M_up_g) {
                    lamg_up = std::max(0.0, Rm.get<volScalarField>("kr_g")->operator[](im)) /
                        std::max(1e-30, mug_m->operator[](im));
                    (*lamGu)[gseg] = lamg_up;

                }
                else {
                    lamg_up = std::max(0.0, krg_f_loc[gseg]) /
                        std::max(1e-30, (mug_f ? (*mug_f)[gseg] : mug_m->operator[](im)));
                    (*lamGu)[gseg] = lamg_up;

                }
            }

            // 最终 CI
            (*CIw)[gseg] = geomCI * knm * lamw_up;
            (*CIg)[gseg] = geomCI * knm * lamg_up;
        }
    }
}

// ------------- 裂缝-裂缝 TI（交点 Star–Delta，多臂通用） -------------
void updateFractureFractureTI(MeshManager& mgr,
    FieldRegistry& Rf,
    const VGParams& vg, const RelPermParams& rp,
    bool /*include_gravity*/, double /*g*/)
{
    auto& net = mgr.fracture_network();
    auto idx = buildFracElemIndex(net);
    ensureFracTIDiagFields(Rf, idx.total);

    auto Sw_f = Rf.get<volScalarField>("Sf_w");
    auto pfw_f = Rf.get<volScalarField>("pf_w");
    auto muw_f = Rf.get<volScalarField>("fr_mu_w");
    auto mug_f = Rf.get<volScalarField>("fr_mu_g");
    auto bf = Rf.get<volScalarField>("fr_aperture");

    //用于检查
    auto Gwarm = Rf.get<volScalarField>("diag_Gw_arm");
    auto Ggarm = Rf.get<volScalarField>("diag_Gg_arm");
    auto ktb = Rf.get<volScalarField>("diag_ktb");
    auto kt_diag = Rf.get<volScalarField>("diag_kt");
    auto aGeom = Rf.has("alpha_geom") ? Rf.get<volScalarField>("alpha_geom") : nullptr; // 已由上一步创建

    if (!Sw_f || !pfw_f || !muw_f || !mug_f || !bf) {
        std::cerr << "[Coupling] TI: fields missing.\n";
        return;
    }

    // 1) 每段“臂导纳” G（相分开算）
    vector<double> Gw(idx.total, 0.0), Gg(idx.total, 0.0);
    for (size_t f = 0; f < net.fractures.size(); ++f) {
        const auto& F = net.fractures[f];
        const auto& I = F.intersections;
        size_t base = idx.offset[f];

        for (size_t e = 0; e < F.elements.size(); ++e) {
            const auto& E = F.elements[e];
            size_t gseg = base + e;

            Vector t = unit(I[E.id].point - I[E.id - 1].point);
            double kt = kt_project(Rf, gseg, t);  // 切向渗透率（已带兜底）
            (*kt_diag)[gseg] = kt ;

            double alpha = E.geomAlpha;              // 几何 alpha（=2/L 或你的定义）
            if (aGeom) (*aGeom)[gseg] = E.geomAlpha; // 再保险


            double b = (*bf)[gseg];              // 开度

            (*ktb)[gseg] = kt * b;

            double krw, krg; kr_vG(vg, rp, (*Sw_f)[gseg], krw, krg);
            double lamw = krw / std::max(1e-30, (*muw_f)[gseg]);
            double lamg = krg / std::max(1e-30, (*mug_f)[gseg]);

            // 注意乘以开度 b（2D 单位厚度）
            Gw[gseg] = alpha * (kt * b) * lamw;
            (*Gwarm)[gseg] = Gw[gseg];
            Gg[gseg] = alpha * (kt * b) * lamg;
            (*Ggarm)[gseg] = Gg[gseg];
        }
    }

    // 2) 收集每个“交点ID”的入射半臂
    int maxID = 0;
    for (const auto& gp : net.globalFFPts) if (gp.id > maxID) maxID = gp.id;
    std::vector<std::vector<std::pair<size_t, bool>>> arms(maxID + 1);

    for (size_t f = 0; f < net.fractures.size(); ++f) {
        const auto& F = net.fractures[f];
        size_t base = idx.offset[f];
        for (size_t e = 0; e < F.elements.size(); ++e) {
            const auto& E = F.elements[e];
            size_t gseg = base + e;
            if (E.isFFatStart) arms[E.gIDstart].push_back({ gseg, true });
            if (E.isFFatEnd)   arms[E.gIDend].push_back({ gseg, false });
        }
    }

    // 3) 清零旧 TI
    for (auto& F : net.fractures) {
        for (auto& E : F.elements) {
            for (auto& ex : E.ffEx) { ex.TIw = 0.0; ex.TIg = 0.0; }
        }
    }

    // 4) 逐交点做 Star–Delta： Tij = Gi*Gj / sum(G)
    for (size_t gid = 0; gid < arms.size(); ++gid) {
        const auto& A = arms[gid];
        if (A.size() < 2) continue;

        double sumGw = 0.0, sumGg = 0.0;
        for (size_t ii = 0; ii < A.size(); ++ii) {
            size_t gs = A[ii].first;   // bool tail = A[ii].second; // 如需
            sumGw += Gw[gs];
            sumGg += Gg[gs];
        }
        if (sumGw <= 1e-30 && sumGg <= 1e-30) continue;

        // 写回到两两“边”
        auto locate_and_set = [&](size_t gA, size_t gB, double Tw, double Tg, int gid) {
            // 找 gA 属于哪条裂缝/哪个段
            size_t fA = 0, eA = 0;
            for (; fA < net.fractures.size(); ++fA) {
                size_t baseA = idx.offset[fA];
                size_t nA = net.fractures[fA].elements.size();
                if (gA < baseA + nA) { eA = gA - baseA; break; }
            }
            auto& EL = net.fractures[fA].elements[eA];
            bool hit = false;
            for (auto& ex : EL.ffEx) {
                if ((size_t)ex.peerGlobalSeg == gB) { ex.TIw = Tw; ex.TIg = Tg;  ex.atGlobalFF = gid; hit = true; break;}
            }
            if (!hit) {
                FractureElement::FFF_Exchange ex;
                ex.peerGlobalSeg = (int)gB;
                ex.TIw = Tw; ex.TIg = Tg;
                ex.atGlobalFF = gid;
                // 如需同时维护 peerFracID / peerSegID，可在此反查赋值
                EL.ffEx.push_back(ex);
            }
            };

        for (size_t u = 0; u < A.size(); ++u) {
            for (size_t v = u + 1; v < A.size(); ++v) {
                size_t gi = A[u].first, gj = A[v].first;
                double T_w = (sumGw > 0.0) ? (Gw[gi] * Gw[gj] / sumGw) : 0.0;
                double T_g = (sumGg > 0.0) ? (Gg[gi] * Gg[gj] / sumGg) : 0.0;
                locate_and_set(gi, gj, T_w, T_g, (int)gid);
                locate_and_set(gj, gi, T_w, T_g, (int)gid);
            }
        }
    }
}

void printCI_Diagnostics( MeshManager& mgr,
    const FieldRegistry& Rm,      // 这里只是为了接口统一，不读取
    const FieldRegistry& Rf,
    size_t maxFracs,
    size_t maxSegsPerFrac)
{
    const auto& net = mgr.fracture_network();
    FracElemIndex idx = buildFracElemIndex(net);

    auto CIw = Rf.get<volScalarField>("CIw");
    auto CIg = Rf.get<volScalarField>("CIg");
    auto knM = Rf.get<volScalarField>("diag_knM");
    auto lamWu = Rf.get<volScalarField>("diag_lamW_up");
    auto lamGu = Rf.get<volScalarField>("diag_lamG_up");
    auto aGeom = Rf.get<volScalarField>("alpha_geom");

    std::cout << "\n==== CI Diagnostics (READ-ONLY) ====\n";
    for (size_t f = 0; f < net.fractures.size(); ++f) {
        if (maxFracs != (size_t)-1 && f >= maxFracs) break;
        const auto& F = net.fractures[f];
        size_t base = idx.offset[f];
        std::cout << "Fracture #" << (f + 1) << " (" << F.start.m_x << "," << F.start.m_y << ") → ("
            << F.end.m_x << "," << F.end.m_y << ")\n";
        for (size_t e = 0; e < F.elements.size(); ++e) {
            if (maxSegsPerFrac != (size_t)-1 && e >= maxSegsPerFrac) break;
            const auto& E = F.elements[e];
            size_t g = base + e;
            std::cout << " Seg " << std::setw(2) << E.id
                << " | cell=" << std::setw(3) << E.cellID
                << " | L=" << std::fixed << std::setprecision(6) << E.length
                << " d=" << E.avgDistance
                << " | geomCI=" << std::scientific << E.geomCI
                << " | kn=" << (knM ? (*knM)[g] : 0.0)
                << " lamW=" << (lamWu ? (*lamWu)[g] : 0.0)
                << " lamG=" << (lamGu ? (*lamGu)[g] : 0.0)
                << " | CIw=" << (CIw ? (*CIw)[g] : 0.0)
                << " CIg=" << (CIg ? (*CIg)[g] : 0.0)
                << " | alpha=" << (aGeom ? (*aGeom)[g] : E.geomAlpha)
                << "\n";
        }
    }
}

//void exportCI_CSV( MeshManager& mgr,
//    const FieldRegistry& Rf,
//    const std::string& outdir,
//    int step, double time)
//{
//    const auto& net = mgr.fracture_network();
//    FracElemIndex idx = buildFracElemIndex(net);
//
//    auto CIw = Rf.get<volScalarField>("CIw");
//    auto CIg = Rf.get<volScalarField>("CIg");
//    auto knM = Rf.get<volScalarField>("diag_knM");
//    auto lamWu = Rf.get<volScalarField>("diag_lamW_up");
//    auto lamGu = Rf.get<volScalarField>("diag_lamG_up");
//    auto aGeom = Rf.get<volScalarField>("alpha_geom");
//
//    char sn[32]; std::snprintf(sn, sizeof(sn), "%05d", step);
//    std::ofstream f(outdir + "/CI_diag_step" + sn + ".csv");
//    f << "time,fid,sid,cellID,L,d,geomCI,kn,lamW,lamG,CIw,CIg,alpha\n";
//    f << std::scientific << std::setprecision(10);
//
//    for (size_t fid = 0; fid < net.fractures.size(); ++fid) {
//        size_t base = idx.offset[fid];
//        const auto& F = net.fractures[fid];
//        for (size_t e = 0; e < F.elements.size(); ++e) {
//            const auto& E = F.elements[e];
//            size_t g = base + e;
//            f << time << "," << (fid + 1) << "," << E.id << "," << E.cellID << ","
//                << E.length << "," << E.avgDistance << "," << E.geomCI << ","
//                << (knM ? (*knM)[g] : 0.0) << ","
//                << (lamWu ? (*lamWu)[g] : 0.0) << ","
//                << (lamGu ? (*lamGu)[g] : 0.0) << ","
//                << (CIw ? (*CIw)[g] : 0.0) << ","
//                << (CIg ? (*CIg)[g] : 0.0) << ","
//                << (aGeom ? (*aGeom)[g] : E.geomAlpha)
//                << "\n";
//        }
//    }
//    std::cout << "[Diag] wrote CI CSV.\n";
//}


static inline void g2fe(const FractureNetwork& net, const FracElemIndex& idx,
    size_t g, size_t& fid, size_t& eid) {
    fid = 0;
    while (fid + 1 < idx.offset.size() && idx.offset[fid + 1] <= g) ++fid;
    eid = g - idx.offset[fid];
}

void printTI_Diagnostics( MeshManager& mgr,
    const FieldRegistry& Rf,
    size_t maxFracs)
{
    const auto& net = mgr.fracture_network();
    FracElemIndex idx = buildFracElemIndex(net);
    auto aGeom = Rf.has("alpha_geom") ? Rf.get<volScalarField>("alpha_geom") : nullptr;

    std::cout << "\n==== TI Diagnostics (READ-ONLY) ====\n";
    for (size_t f = 0; f < net.fractures.size(); ++f) {
        if (maxFracs != (size_t)-1 && f >= maxFracs) break;
        const auto& F = net.fractures[f];
        size_t base = idx.offset[f];

        for (size_t e = 0; e < F.elements.size(); ++e) {
            const auto& E = F.elements[e];
            if (E.ffEx.empty()) continue;
            size_t gA = base + e;
            double aA = aGeom ? (*aGeom)[gA] : E.geomAlpha;

            for (const auto& ex : E.ffEx) {
                size_t fidB = 0, eidB = 0; g2fe(net, idx, (size_t)ex.peerGlobalSeg, fidB, eidB);
                const auto& EB = net.fractures[fidB].elements[eidB];
                size_t gB = idx.offset[fidB] + eidB;
                double aB = aGeom ? (*aGeom)[gB] : EB.geomAlpha;

                // 交点信息（若 atGlobalFF 有效）
                std::string atFF = "N/A";
                if (ex.atGlobalFF >= 0) {
                    if (const auto* G = findGlobalFFByID(net, ex.atGlobalFF)) {
                        std::ostringstream oss;
                        oss << "ID " << G->id << " @(" << G->point.m_x << "," << G->point.m_y << ")";
                        atFF = oss.str();
                    }
                    else {
                        atFF = std::string("ID ") + std::to_string(ex.atGlobalFF) + " (not found)";
                    }
                }

                std::cout << " [FF] Frac#" << (f + 1) << " Seg" << E.id
                    << "  <->  Frac#" << (fidB + 1) << " Seg" << EB.id
                    << " | at " << atFF
                    << " | alphaA=" << aA << " alphaB=" << aB
                    << " | TIw=" << ex.TIw << " TIg=" << ex.TIg
                    << "\n";
            }
        }
    }
}