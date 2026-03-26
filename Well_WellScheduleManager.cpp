/**
 * @file Well_WellScheduleManager.cpp
 * @brief WAG 调度与井控配置管理器实现
 */
#include "Well_WellScheduleManager.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <cctype>

std::string WellScheduleManager::Trim(const std::string& str) {
    if (str.empty()) return str;
    size_t first = str.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, (last - first + 1));
}

std::string WellScheduleManager::ToLower(std::string str) {
    // [Patch 5] 安全转换，防非ASCII越界UB
    std::transform(str.begin(), str.end(), str.begin(), [](unsigned char c) { return std::tolower(c); });
    return str;
}

bool WellScheduleManager::LoadFromCsv(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "[WellScheduleManager] Error: Cannot open CSV file: " << filepath << std::endl;
        return false;
    }

    std::string line;
    if (!std::getline(file, line)) {
        std::cerr << "[WellScheduleManager] Error: CSV file is empty." << std::endl;
        return false;
    }

    schedules_.clear();
    int lineNumber = 1;

    while (std::getline(file, line)) {
        lineNumber++;
        std::string trimmedLine = Trim(line);
        if (trimmedLine.empty() || trimmedLine[0] == '#') continue;

        std::stringstream ss(trimmedLine);
        std::string token;
        std::vector<std::string> tokens;

        while (std::getline(ss, token, ',')) {
            tokens.push_back(Trim(token));
        }

        // [Patch 6] 严格要求 15 列
        if (tokens.size() != 15) {
            std::cerr << "[WellScheduleManager] Error: Invalid column count at line " << lineNumber
                << ". Expected exactly 15, got " << tokens.size() << ". Fail-fast." << std::endl;
            return false;
        }

        try {
            WellScheduleStep step;
            step.t_start = std::stod(tokens[0]);
            step.t_end = std::stod(tokens[1]);
            step.well_name = tokens[2];
            if (step.well_name.empty()) {
                std::cerr << "[WellScheduleManager] Error: Empty well_name at line " << lineNumber << ". Fail-fast." << std::endl;
                return false;
            }
            step.domain = ParseDomain(tokens[3]);
            step.control_mode = ParseControlMode(tokens[4]);
            step.target_value = std::stod(tokens[5]);
            step.component_mode = ParseComponentMode(tokens[6]);
            step.rw = std::stod(tokens[7]);
            step.skin = std::stod(tokens[8]);
            step.well_axis = ParseAxis(tokens[9]);
            step.completion_id = std::stoi(tokens[10]);
            step.wi_override = std::stod(tokens[11]);
            step.L_override = std::stod(tokens[12]);
            step.frac_w = std::stod(tokens[13]);
            step.frac_g = std::stod(tokens[14]);

            std::string err;
            if (!ValidateStep(step, err)) {
                std::cerr << "[WellScheduleManager] Validation Error at line " << lineNumber << ": " << err << std::endl;
                return false;
            }
            schedules_.push_back(step);
        }
        catch (const std::exception& e) {
            std::cerr << "[WellScheduleManager] Exception parsing data at line " << lineNumber << ": " << e.what() << std::endl;
            return false;
        }
    }
    return Validate();
}

bool WellScheduleManager::ValidateStep(const WellScheduleStep& step, std::string& err) const {
    if (step.t_start >= step.t_end) { err = "t_start >= t_end"; return false; }
    if (step.rw <= 0.0) { err = "rw <= 0.0"; return false; }
    if (step.completions.empty()) {
        if (step.completion_id < 0) { err = "completion_id < 0"; return false; }
    }
    else {
        for (const auto& c : step.completions) {
            if (c.completion_id < 0) { err = "completion list contains completion_id < 0"; return false; }
        }
    }

    // [Patch 7] 物理有效性及无用参数豁免
    if (step.control_mode == WellControlMode::BHP) {
        if (std::isnan(step.target_value) || std::isinf(step.target_value)) {
            err = "target_value for BHP is NaN or Inf"; return false;
        }
    }

    if (step.component_mode == WellComponentMode::Total) {
        if (step.frac_w < 0.0 || step.frac_w > 1.0 || step.frac_g < 0.0 || step.frac_g > 1.0) {
            err = "frac_w and frac_g must be in [0, 1]"; return false;
        }
        if (std::abs(step.frac_w + step.frac_g - 1.0) > 1e-8) {
            err = "frac_w + frac_g != 1.0 in Total mode"; return false;
        }
    }

    if (std::abs(step.wi_override) < 1e-12) {
        std::cerr << "[WellScheduleManager] Warning: wi_override is 0.0 for well " << step.well_name
            << ". Usually lacks physical meaning." << std::endl;
    }

    return true;
}

bool WellScheduleManager::Validate() const {
    auto sorted_scheds = schedules_;
    std::sort(sorted_scheds.begin(), sorted_scheds.end(), [](const WellScheduleStep& a, const WellScheduleStep& b) {
        if (a.well_name != b.well_name) return a.well_name < b.well_name;
        return a.t_start < b.t_start;
        });

    for (size_t i = 1; i < sorted_scheds.size(); ++i) {
        if (sorted_scheds[i].well_name == sorted_scheds[i - 1].well_name) {
            if (sorted_scheds[i].t_start < sorted_scheds[i - 1].t_end) {
                std::cerr << "[WellScheduleManager] Time overlap error for well " << sorted_scheds[i].well_name
                    << " around time " << sorted_scheds[i].t_start << std::endl;
                return false;
            }
        }
    }
    return true;
}

std::vector<WellScheduleStep> WellScheduleManager::GetActiveSteps(double current_time) const {
    std::vector<WellScheduleStep> active_steps;
    for (const auto& step : schedules_) {
        if (current_time >= step.t_start && current_time < step.t_end) {
            active_steps.push_back(step);
        }
    }
    return active_steps;
}

WellTargetDomain WellScheduleManager::ParseDomain(const std::string& str) const {
    std::string s = ToLower(str);
    if (s == "matrix") return WellTargetDomain::Matrix;
    if (s == "fracture") return WellTargetDomain::Fracture;
    throw std::invalid_argument("Unknown domain: " + str);
}

WellControlMode WellScheduleManager::ParseControlMode(const std::string& str) const {
    std::string s = ToLower(str);
    if (s == "bhp") return WellControlMode::BHP;
    if (s == "rate") return WellControlMode::Rate;
    throw std::invalid_argument("Unknown control mode: " + str);
}

WellComponentMode WellScheduleManager::ParseComponentMode(const std::string& str) const {
    std::string s = ToLower(str);
    if (s == "water") return WellComponentMode::Water;
    if (s == "gas") return WellComponentMode::Gas;
    if (s == "total") return WellComponentMode::Total;
    throw std::invalid_argument("Unknown component mode: " + str);
}

WellAxis WellScheduleManager::ParseAxis(const std::string& str) const {
    std::string s = ToLower(str);
    if (s == "x") return WellAxis::X;
    if (s == "y") return WellAxis::Y;
    if (s == "z") return WellAxis::Z;
    if (s == "none") return WellAxis::None;
    throw std::invalid_argument("Unknown axis: " + str);
}

