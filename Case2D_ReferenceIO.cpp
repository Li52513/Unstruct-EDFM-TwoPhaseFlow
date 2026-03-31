#include "Case2D_ReferenceIO.h"

#include <chrono>
#include <functional>
#include <system_error>
#include <thread>

#ifdef _WIN32
#include <Windows.h>
#endif

namespace Case2DReferenceIO {

std::string Trim(std::string value) {
    static const char* kWhitespace = " \t\r\n";
    const std::size_t begin = value.find_first_not_of(kWhitespace);
    if (begin == std::string::npos) return "";
    const std::size_t end = value.find_last_not_of(kWhitespace);
    return value.substr(begin, end - begin + 1);
}

std::vector<std::string> SplitCsvLine(const std::string& line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string token;
    while (std::getline(ss, token, ',')) out.push_back(Trim(token));
    return out;
}

CsvTable ReadCsvTable(const std::string& path) {
    std::ifstream in(std::filesystem::path(path), std::ios::in);
    if (!in.good()) throw std::runtime_error("[Case2DReferenceIO] failed to read CSV: " + path);
    CsvTable table;
    std::string line;
    if (!std::getline(in, line)) {
        throw std::runtime_error("[Case2DReferenceIO] empty CSV: " + path);
    }
    table.headers = SplitCsvLine(line);
    for (std::size_t i = 0; i < table.headers.size(); ++i) {
        table.index_by_header[table.headers[i]] = i;
    }
    while (std::getline(in, line)) {
        if (!Trim(line).empty()) table.rows.push_back(SplitCsvLine(line));
    }
    return table;
}

double CsvGetDouble(const CsvTable& table, std::size_t row, const std::string& column) {
    const auto it = table.index_by_header.find(column);
    if (it == table.index_by_header.end()) {
        throw std::runtime_error("[Case2DReferenceIO] missing CSV column: " + column);
    }
    if (row >= table.rows.size() || it->second >= table.rows[row].size()) {
        throw std::runtime_error("[Case2DReferenceIO] malformed CSV row for column: " + column);
    }
    return std::stod(table.rows[row][it->second]);
}

std::string BoolString(bool value) {
    return value ? "true" : "false";
}

std::string PathToGenericString(const std::filesystem::path& path) {
    return path.generic_string();
}

std::string BuildOutputOpenDiagnostics(const std::filesystem::path& targetPath) {
    std::error_code ec;
    std::ostringstream oss;
    oss << "target=" << PathToGenericString(targetPath);
    if (!targetPath.parent_path().empty()) {
        oss << ", parent=" << PathToGenericString(targetPath.parent_path())
            << ", parent_exists=" << BoolString(std::filesystem::exists(targetPath.parent_path(), ec));
        ec.clear();
    }
    oss << ", target_exists=" << BoolString(std::filesystem::exists(targetPath, ec));
    ec.clear();
    const std::filesystem::path cwd = std::filesystem::current_path(ec);
    if (!ec) oss << ", cwd=" << PathToGenericString(cwd);
    return oss.str();
}

#ifdef _WIN32
std::wstring ToExtendedLengthPath(const std::filesystem::path& path) {
    std::wstring native = path.native();
    if (native.rfind(L"\\\\?\\", 0) == 0) return native;
    if (native.rfind(L"\\\\", 0) == 0) return L"\\\\?\\UNC\\" + native.substr(2);
    return L"\\\\?\\" + native;
}
#endif

std::filesystem::path BuildShortAsciiStagingPath(const std::string& targetPath) {
    const std::filesystem::path stagingRoot = std::filesystem::temp_directory_path() / "edfm_case2d_refio";
    std::error_code ec;
    std::filesystem::create_directories(stagingRoot, ec);

    std::ostringstream oss;
    oss << "stage_" << std::hex << std::hash<std::string>{}(targetPath) << ".csv";
    return stagingRoot / oss.str();
}

std::ofstream OpenAsciiStagingStream(const std::string& targetPath,
                                     const char* context,
                                     std::filesystem::path& stagingPath) {
    stagingPath = BuildShortAsciiStagingPath(targetPath);
    std::ofstream out;
    for (int attempt = 0; attempt < 8; ++attempt) {
        std::error_code ec;
        std::filesystem::remove(stagingPath, ec);
        out = std::ofstream(stagingPath, std::ios::out | std::ios::trunc);
        if (out.good()) return out;
        out.clear();
        std::this_thread::sleep_for(std::chrono::milliseconds(40));
    }
    throw std::runtime_error(
        std::string(context) + ": failed to open ASCII staging file for " + targetPath +
        " | staging=" + PathToGenericString(stagingPath) +
        " | " + BuildOutputOpenDiagnostics(std::filesystem::path(targetPath)));
}

void CommitAsciiStagingFile(std::ofstream& out,
                            const std::filesystem::path& stagingPath,
                            const std::string& targetPath,
                            const char* context) {
    out.flush();
    if (!out.good()) {
        out.close();
        std::error_code ec;
        std::filesystem::remove(stagingPath, ec);
        throw std::runtime_error(
            std::string(context) + ": failed while flushing staging file for " + targetPath +
            " | staging=" + PathToGenericString(stagingPath));
    }
    out.close();

    std::error_code ec;
    std::filesystem::path target(targetPath);
    std::filesystem::path absoluteTarget = std::filesystem::absolute(target, ec);
    if (ec) {
        ec.clear();
        absoluteTarget = std::filesystem::current_path(ec) / target;
    }
    if (!absoluteTarget.parent_path().empty()) {
        std::filesystem::create_directories(absoluteTarget.parent_path(), ec);
        ec.clear();
    }

    bool copied = false;
    for (int attempt = 0; attempt < 8; ++attempt) {
        std::error_code removeEc;
        if (std::filesystem::exists(absoluteTarget, removeEc)) {
            std::filesystem::remove(absoluteTarget, removeEc);
        }
        std::error_code copyEc;
        std::filesystem::copy_file(stagingPath, absoluteTarget, std::filesystem::copy_options::none, copyEc);
        if (!copyEc) {
            copied = true;
            break;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(40));
    }
#ifdef _WIN32
    if (!copied) {
        const std::wstring stagingW = ToExtendedLengthPath(stagingPath);
        const std::wstring targetW = ToExtendedLengthPath(absoluteTarget);
        for (int attempt = 0; attempt < 8; ++attempt) {
            DeleteFileW(targetW.c_str());
            if (CopyFileW(stagingW.c_str(), targetW.c_str(), FALSE) != 0) {
                copied = true;
                break;
            }
            Sleep(40);
        }
    }
#endif
    if (!copied) {
        std::error_code cleanupEc;
        std::filesystem::remove(stagingPath, cleanupEc);
        throw std::runtime_error(
            std::string(context) + ": failed to commit staging file to " + targetPath +
            " | staging=" + PathToGenericString(stagingPath) +
            " | abs_target=" + PathToGenericString(absoluteTarget) +
            " | " + BuildOutputOpenDiagnostics(absoluteTarget));
    }
    std::filesystem::remove(stagingPath, ec);
}

ProfileReferenceTable LoadProfileReference(const std::string& path) {
    const CsvTable table = ReadCsvTable(path);
    ProfileReferenceTable ref;
    for (std::size_t row = 0; row < table.rows.size(); ++row) {
        ProfileReferenceRow item;
        item.station_id = static_cast<int>(std::llround(CsvGetDouble(table, row, "station_id")));
        item.target_axis_m = CsvGetDouble(table, row, "target_axis_m");
        item.target_x = CsvGetDouble(table, row, "target_x_m");
        item.target_y = CsvGetDouble(table, row, "target_y_m");
        item.target_time_s = CsvGetDouble(table, row, "target_time_s");
        item.p_ref = CsvGetDouble(table, row, "p_ref_pa");
        item.t_ref = CsvGetDouble(table, row, "t_ref_k");
        ref.rows_by_station_id[item.station_id] = item;
    }
    return ref;
}

} // namespace Case2DReferenceIO
