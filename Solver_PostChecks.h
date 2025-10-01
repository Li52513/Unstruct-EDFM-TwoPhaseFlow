#pragma once
// Solver_PostChecks.h
//
// 依赖：MeshManager.h / FieldRegistry.h / Solver_AssemblerCOO.h
// 功能：装配体检、COO->MatrixMarket 导出、每步 p/T 导出(CSV/TXT)，并在 MSVC 下也保证可创建目录。

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "MeshManager.h"
#include "FieldRegistry.h"
#include "Solver_AssemblerCOO.h"  // Triplet / SparseSystemCOO

// ----------- 目录创建（同时兼容 MSVC 没有 __cplusplus==201703 的情况） -----------
#if __cplusplus >= 201703L
#include <filesystem>
namespace fs = std::filesystem;
inline void ensureParentDirExists(const std::string& filepath) {
    try { fs::create_directories(fs::path(filepath).parent_path()); }
    catch (...) {}
}
inline void ensureDirExists(const std::string& dirpath) {
    try { fs::create_directories(dirpath); }
    catch (...) {}
}
#else
#ifdef _WIN32
#include <direct.h>
inline void ensureParentDirExists(const std::string& filepath) {
    std::string path = filepath;
    // 去掉文件名，保留目录
    size_t pos = path.find_last_of("/\\");
    if (pos == std::string::npos) return;
    path = path.substr(0, pos);
    std::string acc;
    for (size_t i = 0; i < path.size(); ++i) {
        char ch = path[i];
        acc.push_back(ch);
        if (ch == '\\' || ch == '/') {
            _mkdir(acc.c_str()); // 已存在则忽略
        }
    }
    _mkdir(path.c_str());
}
inline void ensureDirExists(const std::string& dirpath) {
    // 递归创建
    std::string acc;
    for (char ch : dirpath) {
        acc.push_back(ch);
        if (ch == '\\' || ch == '/') _mkdir(acc.c_str());
    }
    _mkdir(dirpath.c_str());
}
#else
#include <sys/stat.h>
#include <sys/types.h>
inline void ensureParentDirExists(const std::string& filepath) {
    std::string path = filepath;
    size_t pos = path.find_last_of("/\\");
    if (pos == std::string::npos) return;
    path = path.substr(0, pos);
    std::string acc;
    for (size_t i = 0; i < path.size(); ++i) {
        char ch = path[i];
        acc.push_back(ch);
        if (ch == '/') mkdir(acc.c_str(), 0755);
    }
    mkdir(path.c_str(), 0755);
}
inline void ensureDirExists(const std::string& dirpath) {
    std::string acc;
    for (char ch : dirpath) {
        acc.push_back(ch);
        if (ch == '/') mkdir(acc.c_str(), 0755);
    }
    mkdir(dirpath.c_str(), 0755);
}
#endif
#endif

namespace PostChecks {

    // ======================== 组装诊断报告 ========================
    struct AssemblyReport {
        int rows = 0, nnz = 0, strongDiag = 0, badOff = 0;
        double minDiag = +1e300, maxDiag = -1e300;
        double bMin = +1e300, bMax = -1e300, bSum = 0.0;
    };

    // ——装配体检：强对角性/非对角正值/对角范围/右端——
    inline AssemblyReport reportAssembly(SparseSystemCOO sys, bool force_compress = true) {
        if (force_compress) sys.compressInPlace(0.0);

        AssemblyReport R;
        R.rows = sys.n;
        R.nnz = static_cast<int>(sys.A.size());

        std::vector<double> diag(sys.n, 0.0), sumOff(sys.n, 0.0);
        for (const auto& t : sys.A) {
            if (t.r == t.c) diag[t.r] += t.v;
            else {
                sumOff[t.r] += std::abs(t.v);
                if (t.v > 0.0) ++R.badOff; // TPFA 期望非对角 ≤ 0
            }
        }
        for (int i = 0; i < sys.n; ++i) {
            R.minDiag = std::min(R.minDiag, diag[i]);
            R.maxDiag = std::max(R.maxDiag, diag[i]);
            if (diag[i] >= sumOff[i]) ++R.strongDiag;
        }
        if (!sys.b.empty()) {
            for (double bi : sys.b) {
                R.bMin = std::min(R.bMin, bi);
                R.bMax = std::max(R.bMax, bi);
                R.bSum += bi;
            }
        }
        return R;
    }

    inline void printAssemblyReport(const AssemblyReport& R, const char* tag = "sys") {
        std::cout << "[Asm:" << tag << "] rows=" << R.rows
            << " nnz=" << R.nnz
            << " diag[min,max]=[" << R.minDiag << "," << R.maxDiag << "]"
            << " strongDiag=" << R.strongDiag << "/" << R.rows
            << " badOff(>0)=" << R.badOff
            << "  b[min,max,sum]=[" << R.bMin << "," << R.bMax << "," << R.bSum << "]\n";
    }

    // ——压缩一致性自检（单元测试用，可选）——
    inline bool selfCheckCompression(const SparseSystemCOO& sys, double tol = 1e-12) {
        struct KeyHash {
            size_t operator()(const std::pair<int, int>& k) const noexcept {
                return (size_t(uint32_t(k.first)) << 32) ^ uint32_t(k.second);
            }
        };
        std::unordered_map<std::pair<int, int>, double, KeyHash> acc;
        acc.reserve(sys.A.size() * 2);
        for (const auto& t : sys.A) acc[{t.r, t.c}] += t.v;

        SparseSystemCOO tmp = sys;
        tmp.compressInPlace(0.0);

        if (tmp.A.size() != acc.size()) {
            std::cerr << "[compress-check] nnz mismatch: compressed=" << tmp.A.size()
                << " unique=" << acc.size() << "\n";
            return false;
        }
        for (const auto& t : tmp.A) {
            auto it = acc.find({ t.r,t.c });
            if (it == acc.end()) {
                std::cerr << "[compress-check] missing key (" << t.r << "," << t.c << ")\n";
                return false;
            }
            double ref = it->second;
            double err = std::abs(t.v - ref);
            double scale = std::max(1.0, std::abs(ref));
            if (err > tol * scale) {
                std::cerr << "[compress-check] value mismatch at (" << t.r << "," << t.c
                    << "): got=" << t.v << " ref=" << ref << " relerr=" << err / scale << "\n";
                return false;
            }
        }
        return true;
    }

    // ======================== Matrix Market 导出 ========================
    inline bool dumpCOO_to_matrix_market(SparseSystemCOO sys,
        const std::string& matrixPath,
        const std::string& rhsPath = "",
        bool compress = true,
        double drop_tol = 0.0)
    {
        if (compress) sys.compressInPlace(drop_tol);

        ensureParentDirExists(matrixPath);
        std::ofstream fm(matrixPath);
        if (!fm) { std::cerr << "[MM] cannot open " << matrixPath << "\n"; return false; }

        fm << "%%MatrixMarket matrix coordinate real general\n";
        fm << sys.n << " " << sys.n << " " << sys.A.size() << "\n";
        fm << std::scientific << std::setprecision(16);
        for (const auto& t : sys.A) {
            // MatrixMarket 1-based 索引
            fm << (t.r + 1) << " " << (t.c + 1) << " " << t.v << "\n";
        }
        fm.close();

        if (!rhsPath.empty()) {
            ensureParentDirExists(rhsPath);
            std::ofstream fb(rhsPath);
            if (!fb) { std::cerr << "[MM] cannot open " << rhsPath << "\n"; return false; }
            fb << std::scientific << std::setprecision(16);
            for (double bi : sys.b) fb << bi << "\n";
            fb.close();
        }
        return true;
    }

    // ======================== p/T 导出（CSV or TXT） ========================

    // ——CSV：用于表格/Excel——
    inline bool export_pT_csv_after_step(MeshManager& mgr,
        const FieldRegistry& reg,
        const std::string& p_field,
        const std::string& T_field,
        int step, double t,
        const std::string& folder = "out")
    {
        // 保证目录存在（MSVC 下也可用）
        ensureDirExists(folder);

        std::ostringstream fn;
        fn << folder << "/pT_step_" << std::setw(6) << std::setfill('0') << step << ".csv";
        ensureParentDirExists(fn.str());

        std::ofstream fo(fn.str());
        if (!fo) { std::cerr << "[CSV] cannot open " << fn.str() << "\n"; return false; }

        auto& mesh = mgr.mesh();
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();
        auto p = reg.get<volScalarField>(p_field.c_str());
        auto T = reg.get<volScalarField>(T_field.c_str());
        if (!p || !T) { std::cerr << "[CSV] missing p/T fields.\n"; return false; }

        fo << "# step," << step << ", time," << std::setprecision(16) << std::scientific << t << "\n";
        fo << "cell_id,cx,cy,cz,volume,p,T\n";
        fo << std::scientific << std::setprecision(16);

        for (const auto& c : cells) {
            if (c.id < 0) continue;
            size_t i = id2idx.at(c.id);
            fo << c.id << ","
                << c.center.m_x << "," << c.center.m_y << "," << c.center.m_z << ","
                << c.volume << ","
                << (*p)[i] << ","
                << (*T)[i] << "\n";
        }
        fo.close();
        return true;
    }

    // ——TXT：便于你说的“简单每步一份 txt”——
    inline bool export_pT_txt_after_step(MeshManager& mgr,
        const FieldRegistry& reg,
        const std::string& p_field,
        const std::string& T_field,
        int step, double t,
        const std::string& folder = "out_txt")
    {
        ensureDirExists(folder);

        std::ostringstream fn;
        fn << folder << "/pT_step_" << std::setw(6) << std::setfill('0') << step << ".txt";
        ensureParentDirExists(fn.str());

        std::ofstream fo(fn.str());
        if (!fo) { std::cerr << "[TXT] cannot open " << fn.str() << "\n"; return false; }

        auto& mesh = mgr.mesh();
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();
        auto p = reg.get<volScalarField>(p_field.c_str());
        auto T = reg.get<volScalarField>(T_field.c_str());
        if (!p || !T) { std::cerr << "[TXT] missing p/T fields.\n"; return false; }

        fo << "# step " << step << "  time " << std::setprecision(16) << std::scientific << t << "\n";
        fo << "# columns: cell_id  cx  cy  cz  volume  p  T\n";
        fo << std::scientific << std::setprecision(16);

        for (const auto& c : cells) {
            if (c.id < 0) continue;
            size_t i = id2idx.at(c.id);
            fo << c.id << "  "
                << c.center.m_x << "  " << c.center.m_y << "  " << c.center.m_z << "  "
                << c.volume << "  "
                << (*p)[i] << "  "
                << (*T)[i] << "\n";
        }
        fo.close();
        return true;
    }

} // namespace PostChecks
