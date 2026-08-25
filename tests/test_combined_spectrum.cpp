#include "CombinedSpectrumAnalyzer.h"

#include <TH1D.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

void FillGaussian(TH1D& histogram, double center, double sigma, double height) {
    for (int bin = 1; bin <= histogram.GetNbinsX(); ++bin) {
        const double charge = histogram.GetXaxis()->GetBinCenter(bin);
        const double distance = (charge - center) / sigma;
        const double content = 3.0 + height * std::exp(-0.5 * distance * distance);
        histogram.SetBinContent(bin, content);
        histogram.SetBinError(bin, std::sqrt(content));
    }
}

} // namespace

int main() {
    TH1D first("combined_test_first", "", 2000, 0.0, 2000.0);
    TH1D second("combined_test_second", "", 2000, 0.0, 2000.0);
    FillGaussian(first, 1000.0, 4.0, 1200.0);
    FillGaussian(second, 1000.0, 4.0, 900.0);

    std::vector<hpge::CalibratedSpectrumInput> inputs{
        {&first, 0.0, 0.5, 0.0},
        {&second, 10.0, 0.49, 0.0}
    };
    std::string error;
    auto combined = hpge::CombinedSpectrumAnalyzer::Combine("synthetic", inputs, error);
    if (!combined || !error.empty() || combined->Integral() <= first.Integral()) {
        std::cerr << "Could not combine calibrated spectra: " << error << '\n';
        return 1;
    }
    const auto quality = hpge::CombinedSpectrumAnalyzer::EvaluatePeak(
        *combined, "synthetic", 500.0, 15.0);
    if (!quality.success || std::abs(quality.fittedEnergy - 500.0) > 0.5 ||
        std::abs(quality.residualKeV) > 0.5 || quality.fwhmKeV < 3.0 ||
        quality.fwhmKeV > 7.0 || quality.resolutionPercent <= 0.0) {
        std::cerr << "Combined peak QA failed: " << quality.status
                  << ", centroid=" << quality.fittedEnergy
                  << ", residual=" << quality.residualKeV
                  << ", FWHM=" << quality.fwhmKeV << '\n';
        return 1;
    }
    const auto selectedRangeQuality = hpge::CombinedSpectrumAnalyzer::EvaluatePeakInRange(
        *combined, "synthetic", 500.0, 493.0, 508.0);
    if (!selectedRangeQuality.success ||
        std::abs(selectedRangeQuality.peakFit.rangeLow - 493.0) > 1e-9 ||
        std::abs(selectedRangeQuality.peakFit.rangeHigh - 508.0) > 1e-9 ||
        std::abs(selectedRangeQuality.fittedEnergy - 500.0) > 0.5) {
        std::cerr << "User-selected combined peak interval was not fitted exactly\n";
        return 1;
    }

    const auto invalid = hpge::CombinedSpectrumAnalyzer::Combine("empty", {}, error);
    if (invalid || error.empty()) {
        std::cerr << "Empty combined-spectrum input was not rejected\n";
        return 1;
    }
    std::cout << "PASS: calibrated spectrum combination, energy residual and FWHM\n";
    return 0;
}
