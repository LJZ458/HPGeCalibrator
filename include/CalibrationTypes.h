#pragma once

#include <string>
#include <vector>

namespace hpge {

enum class AxisOrientation { ChargeOnX, ChargeOnY };

struct ReferencePeak {
    std::string datasetId;
    double charge = 0.0;
    double energy = 0.0;
    std::string label;
};

struct CalibrationPoint {
    std::string datasetId;
    double charge = 0.0;
    double energy = 0.0;
    double chargeError = 0.0;
    bool manual = false;
    double residual = 0.0;
};

struct PeakMatchResult {
    bool success = false;
    std::vector<double> charges;
    std::vector<bool> matched;
    double scale = 1.0;
    double offset = 0.0;
    double score = 0.0;
};

struct CalibrationResult {
    int crystal = -1;
    bool success = false;
    bool needsReview = true;
    std::string status;
    double p0 = 0.0;
    double p1 = 0.0;
    double p2 = 0.0;
    double chi2 = 0.0;
    int ndf = 0;
    double residualRms = 0.0;
    std::vector<CalibrationPoint> points;
};

} // namespace hpge

