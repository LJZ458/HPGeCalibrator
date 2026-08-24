#pragma once

#include "CalibrationTypes.h"

#include <map>
#include <string>

namespace hpge {

class CalibrationExporter {
public:
    static bool WriteCppCoefficientLists(
        const std::string& path,
        const std::map<int, CalibrationResult>& results);
};

} // namespace hpge
