#include "CalibrationExporter.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <string>

int main() {
    std::map<int, hpge::CalibrationResult> results;
    hpge::CalibrationResult crystal0;
    crystal0.crystal = 0;
    crystal0.success = true;
    crystal0.p0 = -1.25;
    crystal0.p1 = 0.625;
    crystal0.p2 = 1.75e-5;
    results[0] = crystal0;
    hpge::CalibrationResult crystal63;
    crystal63.crystal = 63;
    crystal63.success = true;
    crystal63.p0 = 2.5;
    crystal63.p1 = 0.75;
    crystal63.p2 = 2.25e-5;
    results[63] = crystal63;

    const auto path = std::filesystem::temp_directory_path() /
                      "hpge_calibration_export_test.hpp";
    if (!hpge::CalibrationExporter::WriteCppCoefficientLists(path.string(), results)) {
        std::cerr << "Could not write C++ coefficient lists\n";
        return 1;
    }
    std::ifstream input(path);
    const std::string contents((std::istreambuf_iterator<char>(input)),
                               std::istreambuf_iterator<char>());
    std::filesystem::remove(path);
    std::size_t listCount = 0;
    std::size_t position = 0;
    const std::string declaration = "std::array<double, 64>";
    while ((position = contents.find(declaration, position)) != std::string::npos) {
        ++listCount;
        position += declaration.size();
    }
    if (listCount != 3 || contents.find("p0 =") == std::string::npos ||
        contents.find("p1 =") == std::string::npos ||
        contents.find("p2 =") == std::string::npos ||
        contents.find("-1.25") == std::string::npos ||
        contents.find("2.25e-05") == std::string::npos ||
        contents.find("missing") == std::string::npos) {
        std::cerr << "Export did not contain exactly three indexed coefficient lists\n";
        return 1;
    }
    std::cout << "PASS: three C++ coefficient arrays with missing-crystal placeholders\n";
    return 0;
}
