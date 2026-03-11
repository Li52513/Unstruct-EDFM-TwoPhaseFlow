/**
 * @file FIM_TransientSupport.hpp
 * @brief 瞬态全隐式(FIM)仿真测试的通用验证与支持工具
 */

#pragma once

#include "SolverContrlStrName_op.h"

#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>

namespace FIM_Engine {

    /**
     * @brief 验证导出的 VTK 文件是否存在、非空，并且包含必需的标量字段
     * @param filename VTK 文件的绝对或相对路径
     * @param check_sw 如果是两相流，是否必须校验饱和度标量
     */
    inline void VerifyVtkExport(const std::string& filename, bool check_sw) {
        const auto pCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
        const auto tCfg = PhysicalProperties_string_op::TemperatureEquation_String::FIM();
        const auto satCfg = PhysicalProperties_string_op::SaturationEquation_String::FIM();

        std::ifstream file(filename, std::ios::ate | std::ios::binary);
        if (!file.is_open() || file.tellg() == 0) {
            throw std::runtime_error("[FAIL] VTK file export failed or is empty: " + filename);
        }
        file.seekg(0, std::ios::beg);
        std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

        auto hasScalar = [&](const std::string& scalarName) {
            return content.find("SCALARS " + scalarName + " ") != std::string::npos;
            };

        if (!hasScalar("P") && !hasScalar(pCfg.pressure_field) && !hasScalar("p_w")) {
            throw std::runtime_error("[FAIL] VTK missing pressure scalar");
        }
        if (!hasScalar("T") && !hasScalar(tCfg.temperatue_field)) { // 注意底层是 temperatue_field
            throw std::runtime_error("[FAIL] VTK missing temperature scalar");
        }
        if (check_sw && !hasScalar("S_w") && !hasScalar(satCfg.saturation)) {
            throw std::runtime_error("[FAIL] VTK missing saturation scalar");
        }
    }

} // namespace FIM_Engine