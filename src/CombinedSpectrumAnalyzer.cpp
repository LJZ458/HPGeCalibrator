#include "CombinedSpectrumAnalyzer.h"

#include "CalibrationEngine.h"

#include <TAxis.h>
#include <TH1.h>
#include <TH1D.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>

namespace hpge {
namespace {

constexpr double kGaussianFwhm = 2.3548200450309493;

std::pair<double, double> EnergyBounds(const CalibratedSpectrumInput& input) {
    const TAxis* axis = input.spectrum->GetXaxis();
    const double low = axis->GetXmin();
    const double high = axis->GetXmax();
    double minimum = std::min(CombinedSpectrumAnalyzer::CalibratedEnergy(low, input),
                              CombinedSpectrumAnalyzer::CalibratedEnergy(high, input));
    double maximum = std::max(CombinedSpectrumAnalyzer::CalibratedEnergy(low, input),
                              CombinedSpectrumAnalyzer::CalibratedEnergy(high, input));
    if (std::abs(input.p2) > 1e-18) {
        const double vertex = -input.p1 / (2.0 * input.p2);
        if (vertex > low && vertex < high) {
            const double energy = CombinedSpectrumAnalyzer::CalibratedEnergy(vertex, input);
            minimum = std::min(minimum, energy);
            maximum = std::max(maximum, energy);
        }
    }
    return {minimum, maximum};
}

} // namespace

double CombinedSpectrumAnalyzer::CalibratedEnergy(
    double charge, const CalibratedSpectrumInput& input) {
    return input.p0 + input.p1 * charge + input.p2 * charge * charge;
}

std::shared_ptr<TH1D> CombinedSpectrumAnalyzer::Combine(
    const std::string& datasetId,
    const std::vector<CalibratedSpectrumInput>& inputs,
    std::string& error) {
    error.clear();
    std::vector<double> channelWidths;
    double minimumEnergy = std::numeric_limits<double>::infinity();
    double maximumEnergy = -std::numeric_limits<double>::infinity();
    for (const auto& input : inputs) {
        if (!input.spectrum || input.spectrum->GetNbinsX() <= 0) continue;
        const auto bounds = EnergyBounds(input);
        if (!std::isfinite(bounds.first) || !std::isfinite(bounds.second)) continue;
        minimumEnergy = std::min(minimumEnergy, bounds.first);
        maximumEnergy = std::max(maximumEnergy, bounds.second);
        const TAxis* axis = input.spectrum->GetXaxis();
        const double center = 0.5 * (axis->GetXmin() + axis->GetXmax());
        const double derivative = std::abs(input.p1 + 2.0 * input.p2 * center);
        const double width = derivative * axis->GetBinWidth(std::max(1, input.spectrum->GetNbinsX() / 2));
        if (std::isfinite(width) && width > 0.0) channelWidths.push_back(width);
    }
    if (!std::isfinite(minimumEnergy) || !std::isfinite(maximumEnergy) ||
        maximumEnergy <= minimumEnergy || channelWidths.empty()) {
        error = "No successfully calibrated spectra are available for combination";
        return {};
    }
    std::sort(channelWidths.begin(), channelWidths.end());
    const double medianWidth = channelWidths[channelWidths.size() / 2];
    double targetWidth = std::clamp(0.5 * medianWidth, 0.05, 1.0);
    constexpr int maximumBins = 200000;
    int bins = static_cast<int>(std::ceil((maximumEnergy - minimumEnergy) / targetWidth));
    bins = std::clamp(bins, 256, maximumBins);
    targetWidth = (maximumEnergy - minimumEnergy) / bins;
    (void)targetWidth;

    static std::atomic<unsigned long> sequence{0};
    const std::string name = "combined_calibrated_" +
        std::to_string(std::hash<std::string>{}(datasetId)) + "_" +
        std::to_string(sequence.fetch_add(1));
    auto combined = std::shared_ptr<TH1D>(
        new TH1D(name.c_str(), "Combined calibrated spectrum;Energy (keV);Counts",
                 bins, minimumEnergy, maximumEnergy));
    combined->SetDirectory(nullptr);
    std::vector<double> variances(static_cast<std::size_t>(bins + 2), 0.0);
    for (const auto& input : inputs) {
        if (!input.spectrum) continue;
        for (int sourceBin = 1; sourceBin <= input.spectrum->GetNbinsX(); ++sourceBin) {
            const double content = input.spectrum->GetBinContent(sourceBin);
            if (!(content > 0.0) || !std::isfinite(content)) continue;
            double errorValue = input.spectrum->GetBinError(sourceBin);
            const double sourceVariance = errorValue > 0.0 ? errorValue * errorValue : content;
            const TAxis* sourceAxis = input.spectrum->GetXaxis();
            const double chargeLow = sourceAxis->GetBinLowEdge(sourceBin);
            const double chargeHigh = sourceAxis->GetBinUpEdge(sourceBin);
            double energyLow = CalibratedEnergy(chargeLow, input);
            double energyHigh = CalibratedEnergy(chargeHigh, input);
            if (!std::isfinite(energyLow) || !std::isfinite(energyHigh)) continue;
            if (energyHigh < energyLow) std::swap(energyLow, energyHigh);
            const double energyWidth = energyHigh - energyLow;
            if (energyWidth <= std::numeric_limits<double>::epsilon()) {
                const int targetBin = combined->FindFixBin(energyLow);
                if (targetBin >= 1 && targetBin <= bins) {
                    combined->AddBinContent(targetBin, content);
                    variances[static_cast<std::size_t>(targetBin)] += sourceVariance;
                }
                continue;
            }
            const int firstTarget = std::max(1, combined->FindFixBin(energyLow));
            const int lastTarget = std::min(bins, combined->FindFixBin(
                std::nextafter(energyHigh, energyLow)));
            for (int targetBin = firstTarget; targetBin <= lastTarget; ++targetBin) {
                const double overlapLow = std::max(energyLow, combined->GetXaxis()->GetBinLowEdge(targetBin));
                const double overlapHigh = std::min(energyHigh, combined->GetXaxis()->GetBinUpEdge(targetBin));
                if (overlapHigh <= overlapLow) continue;
                const double fraction = (overlapHigh - overlapLow) / energyWidth;
                combined->AddBinContent(targetBin, content * fraction);
                variances[static_cast<std::size_t>(targetBin)] +=
                    sourceVariance * fraction * fraction;
            }
        }
    }
    for (int bin = 1; bin <= bins; ++bin) {
        combined->SetBinError(bin, std::sqrt(variances[static_cast<std::size_t>(bin)]));
    }
    return combined;
}

CombinedPeakQuality CombinedSpectrumAnalyzer::EvaluatePeak(
    const TH1& combinedSpectrum, const std::string& datasetId,
    double expectedEnergy, double halfWindowKeV) {
    if (!(expectedEnergy > 0.0) || !std::isfinite(expectedEnergy)) {
        CombinedPeakQuality quality;
        quality.datasetId = datasetId;
        quality.expectedEnergy = expectedEnergy;
        quality.status = "Expected energy must be positive";
        return quality;
    }
    const double binWidth = combinedSpectrum.GetXaxis()->GetBinWidth(1);
    const double halfWindow = std::max(std::abs(halfWindowKeV), 8.0 * binWidth);
    const double low = std::max(combinedSpectrum.GetXaxis()->GetXmin(), expectedEnergy - halfWindow);
    const double high = std::min(combinedSpectrum.GetXaxis()->GetXmax(), expectedEnergy + halfWindow);
    return EvaluatePeakInRange(combinedSpectrum, datasetId, expectedEnergy, low, high);
}

CombinedPeakQuality CombinedSpectrumAnalyzer::EvaluatePeakInRange(
    const TH1& combinedSpectrum, const std::string& datasetId,
    double expectedEnergy, double rangeLowKeV, double rangeHighKeV) {
    CombinedPeakQuality quality;
    quality.datasetId = datasetId;
    quality.expectedEnergy = expectedEnergy;
    if (!(expectedEnergy > 0.0) || !std::isfinite(expectedEnergy)) {
        quality.status = "Expected energy must be positive";
        return quality;
    }
    if (!std::isfinite(rangeLowKeV) || !std::isfinite(rangeHighKeV) ||
        rangeHighKeV <= rangeLowKeV) {
        quality.status = "Peak-fit range must have two finite increasing limits";
        return quality;
    }
    const double low = std::max(combinedSpectrum.GetXaxis()->GetXmin(), rangeLowKeV);
    const double high = std::min(combinedSpectrum.GetXaxis()->GetXmax(), rangeHighKeV);
    quality.peakFit = CalibrationEngine::FitRadwarePeak(combinedSpectrum, low, high);
    quality.success = quality.peakFit.success;
    quality.status = quality.peakFit.status;
    if (!quality.success) return quality;
    quality.fittedEnergy = quality.peakFit.centroid;
    quality.residualKeV = quality.fittedEnergy - expectedEnergy;
    quality.fwhmKeV = kGaussianFwhm * std::abs(quality.peakFit.sigma);
    quality.resolutionPercent = 100.0 * quality.fwhmKeV / expectedEnergy;
    return quality;
}

} // namespace hpge
