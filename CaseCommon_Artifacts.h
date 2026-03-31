#pragma once

#include "Well_WellControlTypes.h"

#include <string>
#include <vector>

namespace CaseCommon {

enum class CaseStage {
    SolveOnly,
    PrepareReference,
    ValidateOnly,
    FullWorkflow
};

struct CaseArtifactPaths {
    std::string root_dir;
    std::string case_dir;
    std::string studies_dir;
    std::string figures_dir;
    std::string engineering_dir;
    std::string reference_dir;
    std::string report_dir;
    std::string report_scripts_dir;
    std::string engineering_stage_manifest_path;
    std::string reference_contract_path;
    std::string report_status_markdown_path;
};

struct ValidationVariableSpec {
    std::string variable_key;
    std::string display_name;
};

struct ValidationTimeSnapshotSpec {
    std::string snapshot_id;
    double time = 0.0;
};

struct ValidationLineSpec2D {
    std::string line_id;
    double x0 = 0.0;
    double y0 = 0.0;
    double x1 = 0.0;
    double y1 = 0.0;
    std::vector<std::string> variable_keys;
};

struct ValidationPointSpec2D {
    std::string point_id;
    double x = 0.0;
    double y = 0.0;
    std::vector<std::string> variable_keys;
};

struct ValidationLineSpec3D {
    std::string line_id;
    double x0 = 0.0;
    double y0 = 0.0;
    double z0 = 0.0;
    double x1 = 0.0;
    double y1 = 0.0;
    double z1 = 0.0;
    std::vector<std::string> variable_keys;
};

struct ValidationSliceSpec3D {
    std::string slice_id;
    std::string axis;
    double coordinate = 0.0;
    std::vector<std::string> variable_keys;
};

struct ValidationPointSpec3D {
    std::string point_id;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    std::vector<std::string> variable_keys;
};

struct ValidationSpec2D {
    std::vector<ValidationVariableSpec> variables;
    std::vector<ValidationTimeSnapshotSpec> profile_snapshots;
    std::vector<ValidationLineSpec2D> characteristic_lines;
    std::vector<ValidationPointSpec2D> observation_points;
    std::vector<std::string> grid_study_tags;
    std::vector<double> time_step_study_values;
    std::vector<std::string> required_well_series;
};

struct ValidationSpec3D {
    std::vector<ValidationVariableSpec> variables;
    std::vector<ValidationTimeSnapshotSpec> profile_snapshots;
    std::vector<ValidationLineSpec3D> characteristic_lines;
    std::vector<ValidationSliceSpec3D> characteristic_slices;
    std::vector<ValidationPointSpec3D> observation_points;
    std::vector<std::string> grid_study_tags;
    std::vector<double> time_step_study_values;
    std::vector<std::string> required_well_series;
};

struct WellTemplateLocation {
    std::string label;
    double x_fraction = 0.0;
    double y_fraction = 0.5;
    double z_fraction = 0.5;
    WellTargetDomain domain = WellTargetDomain::Matrix;
    WellAxis axis = WellAxis::None;
};

struct WellTemplateSpec {
    std::string well_mode = "injprod";
    std::string control_policy = "injector_rate_producer_bhp";
    std::string injection_fluid = "co2";
    std::string thermal_policy = "none";
    std::string injector_name = "INJ";
    std::string producer_name = "PROD";
    WellTemplateLocation injector_location = {"injector", 0.25, 0.5, 0.5, WellTargetDomain::Matrix, WellAxis::None};
    WellTemplateLocation producer_location = {"producer", 0.75, 0.5, 0.5, WellTargetDomain::Matrix, WellAxis::None};
    WellControlMode injector_control_mode = WellControlMode::Rate;
    WellControlMode producer_control_mode = WellControlMode::BHP;
    WellComponentMode injector_component_mode = WellComponentMode::Gas;
    WellComponentMode producer_component_mode = WellComponentMode::Total;
    double injector_target_value = 0.0;
    double producer_target_value = 0.0;
    double injector_temperature = -1.0;
    bool avoid_direct_fracture_completion = true;
};

struct ValidationToleranceSpec {
    double relative_l2_max = -1.0;
    double relative_linf_max = -1.0;
    bool is_placeholder = true;
};

struct FamilyAcceptancePolicy {
    std::string family_group_key;
    std::vector<std::string> required_profile_variables;
    std::vector<std::string> required_monitor_variables;
    std::vector<std::string> required_well_series_when_present;
    bool require_reference_contour_overlay = true;
    bool require_profile_line_overlay = true;
    bool require_monitor_time_series_overlay = true;
    ValidationToleranceSpec pressure_tolerance;
    ValidationToleranceSpec temperature_tolerance;
    ValidationToleranceSpec saturation_tolerance;
};

void EnsureDirRecursive(const std::string& rawPath);
CaseArtifactPaths BuildArtifactPaths(const std::string& rootDir,
                                     const std::string& caseCode,
                                     const std::string& caseSlug);
const char* ToString(CaseStage stage);
bool ParseCaseStage(const std::string& text, CaseStage* stage);

} // namespace CaseCommon
