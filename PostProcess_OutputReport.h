#pragma once
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>

#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FVM_WellDOF_TwoPhase.h"
#include "CapRelPerm.h"            
#include "TwoPhaseWells_StrictRate.h" 

namespace IMPES_Iteration {

    /// \brief 单个时间步的井 + 边界质量流量记录（一行CSV）
    struct WellBoundaryCSVRow
    {
        double time = 0.0;

        // 注入井
        double pbh_inj = 0.0;
        double pcell_inj = 0.0;
        double Mw_inj = 0.0;
        double Mg_inj = 0.0;
        double fw_inj = 0.0;
        double Tin_inj = 0.0;

        // 采出井
        double pbh_prod = 0.0;
        double pcell_prod = 0.0;
        double Mw_prod = 0.0;
        double Mg_prod = 0.0;
        double fw_prod = 0.0;
        double Tout_prod = 0.0; ///< 先占位，可填 NaN 或 -1

        // 边界流出
        double Mw_bdry_out = 0.0;
        double Mg_bdry_out = 0.0;
        double fw_bdry_out = 0.0;
    };

    /// \brief 如果文件不存在，则写入一行表头；如果已存在则不动。
    inline void ensureWellBoundaryCSVHeader(const std::string& filename)
    {
        std::ifstream fin(filename);
        if (fin.good()) return; // 已存在，不重复写

        std::ofstream fout(filename);
        if (!fout.good()) {
            std::cerr << "[WellCSV] Cannot open file for header: " << filename << "\n";
            return;
        }

        fout << "time,"
            // inj
            << "pbh_inj,pcell_inj,Mw_inj,Mg_inj,fw_inj,Tin_inj,"
            // prod
            << "pbh_prod,pcell_prod,Mw_prod,Mg_prod,fw_prod,Tout_prod,"
            // boundary
            << "Mw_bdry_out,Mg_bdry_out,fw_bdry_out"
            << "\n";
    }

    /// \brief 追加一行记录到 CSV 文件（假定表头已经存在或稍后会手动调用 ensureHeader）
    inline void appendWellBoundaryCSVRow(
        const std::string& filename,
        const WellBoundaryCSVRow& row,
        int time_precision = 8,
        int value_precision = 8
    )
    {
        std::ofstream fout(filename, std::ios::app);
        if (!fout.good()) {
            std::cerr << "[WellCSV] Cannot open file for append: " << filename << "\n";
            return;
        }

        fout << std::setprecision(time_precision) << std::scientific
            << row.time << ",";

        auto printVal = [&](double v, bool last = false)
            {
                fout << std::setprecision(value_precision) << std::scientific << v;
                if (!last) fout << ",";
            };

        // inj
        printVal(row.pbh_inj);
        printVal(row.pcell_inj);
        printVal(row.Mw_inj);
        printVal(row.Mg_inj);
        printVal(row.fw_inj);
        printVal(row.Tin_inj);

        // prod
        printVal(row.pbh_prod);
        printVal(row.pcell_prod);
        printVal(row.Mw_prod);
        printVal(row.Mg_prod);
        printVal(row.fw_prod);
        printVal(row.Tout_prod);

        // boundary
        printVal(row.Mw_bdry_out);
        printVal(row.Mg_bdry_out, /*last=*/false);
        fout << std::setprecision(value_precision) << std::scientific << row.fw_bdry_out;

        fout << "\n";
    }
}