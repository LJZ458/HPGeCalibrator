#pragma once

#include "CalibrationTypes.h"

#include <memory>
#include <string>
#include <vector>

class TH1;
class TH1D;

namespace hpge {

struct CalibratedSpectrumInput {
    const TH1* spectrum = nullptr;
    double p0 = 0.0;
    double p1 = 0.0;
    double p2 = 0.0;
};

class CombinedSpectrumAnalyzer {
public:
    static double CalibratedEnergy(double charge, const CalibratedSpectrumInput& input);

    static std::shared_ptr<TH1D> Combine(
        const std::string& datasetId,
        const std::vector<CalibratedSpectrumInput>& inputs,
        std::string& error);

    static CombinedPeakQuality EvaluatePeak(
        const TH1& combinedSpectrum,
        const std::string& datasetId,
        double expectedEnergy,
        double halfWindowKeV);

    static CombinedPeakQuality EvaluatePeakInRange(
        const TH1& combinedSpectrum,
        const std::string& datasetId,
        double expectedEnergy,
        double rangeLowKeV,
        double rangeHighKeV);
};

} // namespace hpge
