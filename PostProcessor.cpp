// PostProcessor.cpp  —— 无 <filesystem> 版本
#include "PostProcessor.h"
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <string>

#include "Mesh.h"
#include "FractureNetwork.h"
#include "FieldRegistry.h"
#include "VolField.h"
#include "FieldAccess.h"

using std::ofstream;
using ScalarField = VolField<double>;

// ========== 目录创建的跨平台实现（不用 std::filesystem） ==========
#ifdef _WIN32
#include <direct.h>     // _mkdir
#include <sys/stat.h>   // _stat
static inline bool dirExists(const std::string& p) {
    struct _stat st; return (_stat(p.c_str(), &st) == 0) && (st.st_mode & _S_IFDIR);
}
static inline int makeDir(const std::string& p) {
    return _mkdir(p.c_str()); // 已存在返回 -1，无妨
}
#define SEP1 '\\'
#define SEP2 '/'
#else
#include <sys/stat.h>   // mkdir, stat
#include <sys/types.h>
#include <unistd.h>
static inline bool dirExists(const std::string& p) {
    struct stat st; return (stat(p.c_str(), &st) == 0) && S_ISDIR(st.st_mode);
}
static inline int makeDir(const std::string& p) {
    return mkdir(p.c_str(), 0755); // 已存在返回 -1，无妨
}
#define SEP1 '/'
#define SEP2 '/'
#endif

// 递归创建多级目录：e.g. "results/fracture/run01"
static void ensureDir(const std::string& dir) {
    if (dir.empty()) return;

    std::string path;
    size_t i = 0;

#ifdef _WIN32
    // 处理盘符前缀，比如 "C:\..."
    if (dir.size() >= 2 && dir[1] == ':') { path = dir.substr(0, 2); i = 2; }
#endif
    // 处理以根分隔符开头的情况（如 "/tmp"）
    if (i < dir.size() && (dir[i] == SEP1 || dir[i] == SEP2)) {
        path.push_back(SEP1);
        ++i;
    }

    while (i < dir.size()) {
        // 收集一个路径片段，直到遇到分隔符
        size_t j = i;
        while (j < dir.size() && dir[j] != SEP1 && dir[j] != SEP2) ++j;
        std::string part = dir.substr(i, j - i);
        if (!part.empty()) {
            if (!path.empty() && path.back() != SEP1) path.push_back(SEP1);
            path += part;
            if (!dirExists(path)) makeDir(path);
        }
        // 跳过连续分隔符
        i = j + 1;
    }
}

// ============================================================

PostProcessor::PostProcessor(std::string outdir) : outdir_(std::move(outdir)) {
    ensureDir(outdir_);
}

std::string PostProcessor::stepName_(int step) {
    char buf[32]; std::snprintf(buf, sizeof(buf), "%05d", step);
    return std::string(buf);
}

void PostProcessor::appendTimeIndex_(double time, int step) const {
    ofstream ft(outdir_ + "/times.csv", std::ios::app);
    ft << std::fixed << std::setprecision(8) << step << "," << time << "\n";
}

void PostProcessor::exportMatrixValue(const Mesh& mesh,
    const FieldRegistry& reg,
    double time, int step) const {
    exportMatrix_(mesh, reg, nullptr, 0.0, time, step);
    appendTimeIndex_(time, step);
}

void PostProcessor::exportFractureValue(const Mesh& mesh,
    const FractureNetwork& frNet,
    const FieldRegistry& reg_fr,
    double time, int step) const {
    exportFracture_(mesh, frNet, reg_fr, nullptr, 0.0, time, step);
    appendTimeIndex_(time, step);
}


void PostProcessor::exportMatrixValueInterpolated(const Mesh& mesh,
    const FieldRegistry& reg_prev,
    const FieldRegistry& reg_curr,
    double time, int step, double alpha) const {
    exportMatrix_(mesh, reg_prev, &reg_curr, alpha, time, step);
    appendTimeIndex_(time, step);
}

void PostProcessor::exportFractureValueInterpolated(const Mesh& mesh,
    const FractureNetwork& frNet,
    const FieldRegistry& reg_prev_fr,
    const FieldRegistry& reg_curr_fr,
    double time, int step, double alpha) const {
    exportFracture_(mesh, frNet, reg_prev_fr, &reg_curr_fr, alpha, time, step);
    appendTimeIndex_(time, step);
}

// —— 工具：统计全局裂缝段总数 —— //
static size_t countAllSegments(const FractureNetwork& frNet) {
    size_t n = 0;
    for (const auto& F : frNet.fractures) n += F.elements.size();
    return n;
}

// —— 工具：列出“尺寸匹配”的全部标量场名（自动导出所有字段） —— //
static std::vector<std::string> listScalarFieldNamesBySize(
    const FieldRegistry& R, size_t wantedSize)
{
    std::vector<std::string> names;
    names.reserve(R.fields.size());
    for (const auto& kv : R.fields) {
        auto fld = std::dynamic_pointer_cast<ScalarField>(kv.second);
        if (!fld) continue;
        if (fld->data.size() == wantedSize) names.push_back(kv.first);
    }
    std::sort(names.begin(), names.end());
    return names;
}

// ============ 写矩阵（自动列） ============
void PostProcessor::exportMatrix_(const Mesh& mesh,
    const FieldRegistry& regA,
    const FieldRegistry* regB, double alpha,
    double time, int step) const
{
    const auto& cells = const_cast<Mesh&>(mesh).getCells();
    const size_t N = cells.size();

    // 自动列 + 常用列优先
    std::vector<std::string> autoCols = listScalarFieldNamesBySize(regA, N);
    auto takeIfExists = [&](const char* name, std::vector<std::string>& out) {
        auto it = std::find(autoCols.begin(), autoCols.end(), name);
        if (it != autoCols.end()) { out.push_back(*it); autoCols.erase(it); }
        };
    std::vector<std::string> orderedCols;
    const char* preferred[] = { "p_w","S_w","T","p_c","p_g","kr_w","kr_g","C_eff","lambda_eff" };
    for (auto* n : preferred) takeIfExists(n, orderedCols);
    orderedCols.insert(orderedCols.end(), autoCols.begin(), autoCols.end());

    // 写文件
    ensureDir(outdir_); // 保险
    std::string fn = outdir_ + "/state_matrix_step" + stepName_(step) + ".csv";
    ofstream f(fn);
    f << "time,cellID,cx,cy,cz,region,location";
    for (const auto& n : orderedCols) f << "," << n;
    f << "\n";
    f << std::scientific << std::setprecision(10);

    for (size_t i = 0; i < N; ++i) {
        const auto& c = cells[i];
        if (c.id < 0) continue; // 保险

        f << time << ","
            << c.id << ","
            << c.center.m_x << "," << c.center.m_y << "," << c.center.m_z << ","
            << static_cast<int>(c.region) << ","
            << static_cast<int>(c.location);

        for (const auto& name : orderedCols) {
            double val = (!regB)
                ? FieldAccess::getScalar(regA, name, i)
                : FieldAccess::getScalarLerp(regA, *regB, name, i, alpha);
            f << "," << val;
        }
        f << "\n";
    }
}

// ============ 写裂缝（自动列） ============
// 说明：这里的 getG 就是用“全局段索引 g”从裂缝 FieldRegistry 中取值（或插值后取值）
//       g 的定义：按 frac 顺序串接每条裂缝的所有元素（段），与你的 initFracturePrimaries 顺序一致。
void PostProcessor::exportFracture_(const Mesh&,
    const FractureNetwork& frNet,
    const FieldRegistry& regA,
    const FieldRegistry* regB, double alpha,
    double time, int step) const
{
    const size_t Nseg = countAllSegments(frNet);

    // 自动列 + 常用列优先
    std::vector<std::string> autoCols = listScalarFieldNamesBySize(regA, Nseg);
    auto take = [&](const char* name, std::vector<std::string>& out) {
        auto it = std::find(autoCols.begin(), autoCols.end(), name);
        if (it != autoCols.end()) { out.push_back(*it); autoCols.erase(it); }
        };
    std::vector<std::string> orderedCols;
    const char* preferred[] = { "pf_w","Sf_w","Tf" };
    for (auto* n : preferred) take(n, orderedCols);
    orderedCols.insert(orderedCols.end(), autoCols.begin(), autoCols.end());

    ensureDir(outdir_); // 保险
    std::string fn = outdir_ + "/state_fracture_step" + stepName_(step) + ".csv";
    ofstream f(fn);
    f << "time,fid,segid,hostCell,x1,y1,x2,y2,mx,my";
    for (const auto& n : orderedCols) f << "," << n;
    f << "\n";
    f << std::scientific << std::setprecision(10);

    size_t g = 0; // 全局段索引（getG 使用它来从 reg 中取值）
    for (size_t fid = 0; fid < frNet.fractures.size(); ++fid) {
        const auto& F = frNet.fractures[fid];
        const auto& I = F.intersections;

        for (const auto& elem : F.elements) {
            int sid = elem.id;
            const auto& p1 = I[sid - 1].point;
            const auto& p2 = I[sid].point;
            auto mid = 0.5 * (p1 + p2);

            f << time << ","
                << (fid + 1) << "," << sid << "," << elem.cellID << ","
                << p1.m_x << "," << p1.m_y << ","
                << p2.m_x << "," << p2.m_y << ","
                << mid.m_x << "," << mid.m_y;

            for (const auto& name : orderedCols) {
                double val = (!regB)
                    ? FieldAccess::getScalar(regA, name, g)
                    : FieldAccess::getScalarLerp(regA, *regB, name, g, alpha);
                f << "," << val;
            }
            f << "\n";
            ++g; // 下一段
        }
    }
}