#pragma once
// ================== PostProcess_.h ==================
// 声明：cell→face（含边界ABC + 梯度缓冲）、face→node、Tecplot写出

#include <string>
#include <vector>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "TemperatureBCAdapter.h"
#include "BCAdapter.h"

// ------------------ 1) cell → face：模板声明（定义在 .cpp，且做显式实例化） ------------------
template<class BCAdapter>
bool getFaceValueFromCellValue_BC(
    MeshManager& mgr,
    const FieldRegistry& reg,
    FaceFieldRegistry& freg,
    const BCAdapter& bc,
    const std::string& cellFieldName,
    const std::string& faceOutName,
    const std::vector<Vector>* gradBuf
)
{
    Mesh& mesh = const_cast<Mesh&>(mgr.mesh());
    const auto& faces = mesh.getFaces();
    const auto& cells = mesh.getCells();
    const auto& id2idx = mesh.getCellId2Index();

    auto phiC = reg.get<volScalarField>(cellFieldName);
    if (!phiC) { std::cerr << "[getFaceValueFromCellValue_BC] missing cell field: " << cellFieldName << "\n"; return false; }
    if (phiC->size != cells.size()) {
        std::cerr << "[getFaceValueFromCellValue_BC] size mismatch: " << cellFieldName << " vs cells\n"; return false;
    }

    const bool haveGrad = (gradBuf && gradBuf->size() == cells.size());
    auto phiF = freg.getOrCreate<faceScalarField>(faceOutName, faces.size(), 0.0);
    for (const auto& F : faces)
    {
        const int iF = F.id - 1;

        if (!F.isBoundary())
        {
            // —— 内部面：线性插值 —— //
            const int Pid = F.ownerCell, Nid = F.neighborCell;
            const int iP = id2idx.at(Pid), iN = id2idx.at(Nid);

            double gamma = F.f_linearInterpolationCoef;
            if (!(gamma > 0.0 && gamma < 1.0)) gamma = 0.5;

            const double phiP = (*phiC)[iP];
            const double phiN = (*phiC)[iN];
            (*phiF)[iF] = (1.0 - gamma) * phiP + gamma * phiN;
        }
        else
        {
            // —— 边界面：Adapter 取 a,b,c；未设置则自然边界 φf=φP —— //
            double a = 0.0, b = 0.0, c = 0.0;
            const bool hasBC = bc.getABC(F.id, a, b, c);

            const int  Pid = F.ownerCell;
            const int  iP = id2idx.at(Pid);
            const double phiP = (*phiC)[iP];

            if (!hasBC) { (*phiF)[iF] = phiP; continue; }

            const double dperp = std::max(F.ownerToNeighbor.Mag(), 1e-14);
            const double Eabs = F.vectorE.Mag();
            const double Sabs = F.length;

            double grad_dot_T = 0.0;
            if (haveGrad) grad_dot_T = ((*gradBuf)[iP] * F.vectorT);

            const double denom = a * Sabs + b * Eabs / dperp;
            double alpha = 0.0, beta = 0.0;
            if (std::abs(denom) > 1e-30) {
                alpha = (b * Eabs / dperp) / denom;
                beta = (c * Sabs - grad_dot_T * b) / denom;
            }
            else {
                alpha = 1.0; beta = 0.0; // 极限退化
            }
            (*phiF)[iF] = alpha * phiP + beta;

        }

    }

    return true;
}

// ------------------ 2) 便捷封装声明（定义在 .cpp） ------------------
bool getFaceValueFromCellValue_T(
    MeshManager& mgr, const FieldRegistry& reg, FaceFieldRegistry& freg,
    const TemperatureBCAdapter& Tbc,
    const std::string& T_cell, const std::string& T_face_out,
    const std::vector<Vector>* gradT
);

bool getFaceValueFromCellValue_P(
    MeshManager& mgr, const FieldRegistry& reg, FaceFieldRegistry& freg,
    const PressureBCAdapter& Pbc,
    const std::string& p_cell, const std::string& p_face_out,
    const std::vector<Vector>* gradp
);

// ------------------ 3) 无 BC 情况：直接依据单元值插值得到面值 ------------------
bool getFaceValueFromCellValue_plain(
    MeshManager& mgr,
    const FieldRegistry& reg,
    FaceFieldRegistry& freg,
    const std::string& cellFieldName,
    const std::string& faceOutName);

// ------------------ 4) face → node 声明（定义在 .cpp） ------------------
bool getNodeValueFromFaceValue(
    MeshManager& mgr,
    FaceFieldRegistry& freg,
    const std::string& faceFieldName,
    std::vector<double>& nodeValsOut
);

// ------------------ 4) 一条龙 Tecplot 输出 声明（定义在 .cpp） ------------------
bool outputTecplot_cellToFaceToNode_BC(
    MeshManager& mgr,
    const FieldRegistry& reg,
    FaceFieldRegistry& freg,
    const TemperatureBCAdapter* Tbc,
    const PressureBCAdapter* Pbc,
    const std::string& cellFieldName,
    const std::string& faceFieldName,
    const std::vector<Vector>* gradBuf,
    const std::string& outFilename
);


