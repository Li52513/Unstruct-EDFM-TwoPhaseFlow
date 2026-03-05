/**
 * @file Well_WellScheduleManager.h
 * @brief WAG 调度与井控配置管理器
 */
#ifndef WELL_WELLSCHEDULEMANAGER_H
#define WELL_WELLSCHEDULEMANAGER_H

#include <vector>
#include <string>
#include "Well_WellControlTypes.h"

class WellScheduleManager {
public:
    WellScheduleManager() = default;
    ~WellScheduleManager() = default;

    bool LoadFromCsv(const std::string& filepath);
    bool Validate() const;
    std::vector<WellScheduleStep> GetActiveSteps(double current_time) const;

private:
    std::vector<WellScheduleStep> schedules_;

    static std::string Trim(const std::string& str);
    static std::string ToLower(std::string str);
    bool ValidateStep(const WellScheduleStep& step, std::string& err) const;

    WellTargetDomain ParseDomain(const std::string& str) const;
    WellControlMode ParseControlMode(const std::string& str) const;
    WellComponentMode ParseComponentMode(const std::string& str) const;
    WellAxis ParseAxis(const std::string& str) const;
};

#endif // WELL_WELLSCHEDULEMANAGER_H