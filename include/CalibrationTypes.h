#pragma once

#include <string>
#include <vector>

namespace hpge {

enum class AxisOrientation { ChargeOnX, ChargeOnY };

struct PeakFitResult {
    bool success = false;
    std::string status;
    double rangeLow = 0.0;
    double rangeHigh = 0.0;
    double centroid = 0.0;
    double centroidError = 0.0;
    double sigma = 0.0;
    double height = 0.0;
    double tailFraction = 0.0;
    double beta = 0.0;
    double stepFraction = 0.0;
    double background0 = 0.0;
    double background1 = 0.0;
    double background2 = 0.0;
    double chi2 = 0.0;
    int ndf = 0;
};

struct ReferencePeak {
    std::string datasetId;
    double charge = 0.0;
    double energy = 0.0;
    std::string label;
    PeakFitResult peakFit;
};

struct CalibrationPoint {
    std::string datasetId;
    double charge = 0.0;
    double energy = 0.0;
    double chargeError = 0.0;
    bool manual = false;
    double residual = 0.0;
    PeakFitResult peakFit;
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
