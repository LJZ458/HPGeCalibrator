#include "CalibrationExporter.h"

#include <fstream>
#include <iomanip>

namespace hpge {

bool CalibrationExporter::WriteCppCoefficientLists(
    const std::string& path, const std::map<int, CalibrationResult>& results) {
    std::ofstream output(path);
    if (!output) return false;
    output << "#pragma once\n\n#include <array>\n#include <limits>\n\n"
              "namespace hpge_calibration {\n\n"
              "inline const double missing = std::numeric_limits<double>::quiet_NaN();\n\n";
    const auto writeList = [&](const char* name, int parameter) {
        output << "inline const std::array<double, 64> " << name << " = {\n    ";
        for (int crystal = 0; crystal < 64; ++crystal) {
            const auto found = results.find(crystal);
            if (found == results.end() || !found->second.success) {
                output << "missing";
            } else {
                const double value = parameter == 0 ? found->second.p0
                                   : parameter == 1 ? found->second.p1
                                                    : found->second.p2;
                output << std::setprecision(12) << value;
            }
            if (crystal != 63) output << ", ";
            if ((crystal + 1) % 8 == 0 && crystal != 63) output << "\n    ";
        }
        output << "\n};\n\n";
    };
    writeList("p0", 0);
    writeList("p1", 1);
    writeList("p2", 2);
    output << "} // namespace hpge_calibration\n";
    return static_cast<bool>(output);
}

} // namespace hpge
