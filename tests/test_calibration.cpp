#include "CalibrationEngine.h"

#include <TH1D.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

bool Near(double actual, double expected, double tolerance, const std::string& description) {
    if (std::abs(actual - expected) <= tolerance) return true;
    std::cerr << description << ": expected " << expected << ", got " << actual << '\n';
    return false;
}

hpge::CalibrationPoint MakeCalibrationPoint(double charge, double energy,
                                            const std::string& dataset = "test") {
    hpge::CalibrationPoint point;
    point.datasetId = dataset;
    point.charge = charge;
    point.energy = energy;
    return point;
}

TH1D MakeSpectrum(const char* name, const std::vector<double>& centers,
                  double scale = 1.0, double offset = 0.0) {
    TH1D histogram(name, name, 4096, 0.0, 4096.0);
    for (int bin = 1; bin <= histogram.GetNbinsX(); ++bin) {
        const double x = histogram.GetBinCenter(bin);
        double value = 2.0 + 0.0002 * x;
        for (double center : centers) {
            const double mapped = offset + scale * center;
            value += 1000.0 * std::exp(-0.5 * std::pow((x - mapped) / 3.0, 2));
        }
        histogram.SetBinContent(bin, value);
    }
    histogram.SetEntries(100000.0);
    return histogram;
}

TH1D MakeVariableIntensitySpectrum(const char* name,
                                   const std::vector<double>& centers,
                                   const std::vector<double>& amplitudes,
                                   double scale, double offset,
                                   const std::vector<double>& extraPeaks = {},
                                   double quadratic = 0.0) {
    TH1D histogram(name, name, 4096, 0.0, 4096.0);
    for (int bin = 1; bin <= histogram.GetNbinsX(); ++bin) {
        const double x = histogram.GetBinCenter(bin);
        double value = 4.0 + 0.001 * x + 0.7 * std::sin(0.021 * x);
        for (std::size_t index = 0; index < centers.size(); ++index) {
            const double mapped = offset + scale * centers[index] +
                                  quadratic * centers[index] * centers[index];
            value += amplitudes[index] *
                     std::exp(-0.5 * std::pow((x - mapped) / 3.5, 2));
        }
        for (double extra : extraPeaks) {
            value += 550.0 * std::exp(-0.5 * std::pow((x - extra) / 4.0, 2));
        }
        histogram.SetBinContent(bin, std::max(value, 0.0));
    }
    histogram.SetEntries(150000.0);
    return histogram;
}

TH1D MakeLowStatisticsSpikySpectrum(const char* name,
                                    const std::vector<double>& centers,
                                    const std::vector<double>& amplitudes,
                                    double scale, double offset,
                                    const std::vector<int>& spikeBins) {
    TH1D histogram(name, name, 4096, 0.0, 4096.0);
    for (int bin = 1; bin <= histogram.GetNbinsX(); ++bin) {
        const double x = histogram.GetBinCenter(bin);
        double value = 1.4 + 0.18 * std::sin(0.019 * x) + 0.12 * std::sin(0.071 * x);
        for (std::size_t index = 0; index < centers.size(); ++index) {
            const double mapped = offset + scale * centers[index];
            value += amplitudes[index] *
                     std::exp(-0.5 * std::pow((x - mapped) / 3.2, 2));
        }
        histogram.SetBinContent(bin, std::max(value, 0.0));
    }
    for (std::size_t index = 0; index < spikeBins.size(); ++index) {
        const int bin = spikeBins[index];
        histogram.SetBinContent(bin, histogram.GetBinContent(bin) +
            32.0 + 7.0 * static_cast<double>(index % 4));
        // A small two-bin skirt makes the detector-search stage see the feature,
        // while its shape remains much narrower than a real HPGe photopeak.
        histogram.SetBinContent(bin - 1, histogram.GetBinContent(bin - 1) + 5.0);
        histogram.SetBinContent(bin + 1, histogram.GetBinContent(bin + 1) + 5.0);
    }
    histogram.SetEntries(4200.0);
    return histogram;
}

void AddGaussianPeak(TH1D& histogram, double center, double amplitude, double sigma) {
    for (int bin = 1; bin <= histogram.GetNbinsX(); ++bin) {
        const double x = histogram.GetBinCenter(bin);
        histogram.SetBinContent(bin, histogram.GetBinContent(bin) +
            amplitude * std::exp(-0.5 * std::pow((x - center) / sigma, 2)));
    }
}

double ChargeForEnergy(double energy, double p0, double p1, double p2) {
    return (-p1 + std::sqrt(p1 * p1 - 4.0 * p2 * (p0 - energy))) / (2.0 * p2);
}

bool TestFindPeaks() {
    auto histogram = MakeSpectrum("find_peaks", {500.0, 1400.0, 2800.0});
    hpge::CalibrationEngine::SearchOptions options;
    options.sigmaBins = 2.0;
    options.threshold = 0.05;
    const auto peaks = hpge::CalibrationEngine::FindPeakCandidates(histogram, options);
    if (peaks.size() != 3) {
        std::cerr << "Expected 3 peak candidates, got " << peaks.size() << '\n';
        return false;
    }
    return Near(peaks[0], 500.0, 2.0, "first candidate") &&
           Near(peaks[1], 1400.0, 2.0, "second candidate") &&
           Near(peaks[2], 2800.0, 2.0, "third candidate");
}

bool TestRefinePeak() {
    TH1D histogram("refine", "refine", 200, 0.0, 200.0);
    const double center = 83.35;
    for (int bin = 1; bin <= histogram.GetNbinsX(); ++bin) {
        const double x = histogram.GetBinCenter(bin);
        histogram.SetBinContent(bin, 2000.0 * std::exp(-0.5 * std::pow((x - center) / 2.5, 2)));
    }
    histogram.SetEntries(10000.0);
    const double refined = hpge::CalibrationEngine::RefinePeak(histogram, 80.0, 10.0);
    return Near(refined, center, 0.2, "refined peak");
}

bool TestMatchMultiplePeaks() {
    const std::vector<double> reference{600.0, 1200.0, 1850.0, 2700.0};
    const double scale = 1.08;
    const double offset = 37.0;
    auto target = MakeSpectrum("match_multiple", reference, scale, offset);
    hpge::CalibrationEngine::SearchOptions options;
    const auto match = hpge::CalibrationEngine::MatchReferencePeaks(target, reference, options);
    if (!match.success || match.matched.size() != reference.size()) {
        std::cerr << "Multiple-peak mapping did not succeed\n";
        return false;
    }
    bool ok = Near(match.scale, scale, 0.01, "mapping scale") &&
              Near(match.offset, offset, 4.0, "mapping offset");
    for (std::size_t i = 0; i < reference.size(); ++i) {
        if (!match.matched[i]) {
            std::cerr << "Reference peak " << i << " was not matched\n";
            ok = false;
        } else {
            ok &= Near(match.charges[i], offset + scale * reference[i], 2.0,
                       "mapped peak " + std::to_string(i));
        }
    }
    return ok;
}

bool TestMatchSinglePeak() {
    auto target = MakeSpectrum("match_single", {1000.0}, 1.0, 25.0);
    hpge::CalibrationEngine::SearchOptions options;
    const auto match = hpge::CalibrationEngine::MatchReferencePeaks(target, {1000.0}, options);
    return match.success && match.matched.size() == 1 && match.matched[0] &&
           Near(match.charges[0], 1025.0, 2.0, "single mapped peak") &&
           Near(match.offset, 25.0, 2.0, "single mapping offset");
}

bool TestRadwarePeakFit() {
    hpge::PeakFitResult expected;
    expected.rangeLow = 850.0;
    expected.rangeHigh = 1150.0;
    expected.centroid = 1002.4;
    expected.sigma = 7.5;
    expected.height = 1800.0;
    expected.tailFraction = 0.08;
    expected.beta = 24.0;
    expected.stepFraction = 0.015;
    expected.background0 = 35.0;
    expected.background1 = 4.0;
    expected.background2 = 2.0;

    TH1D histogram("radware_fit", "radware_fit", 600, 700.0, 1300.0);
    for (int bin = 1; bin <= histogram.GetNbinsX(); ++bin) {
        const double x = histogram.GetBinCenter(bin);
        double counts = hpge::CalibrationEngine::EvaluateRadwarePeak(x, expected);
        counts += 0.35 * std::sqrt(std::max(counts, 1.0)) * std::sin(0.37 * bin);
        histogram.SetBinContent(bin, std::max(counts, 0.0));
    }
    histogram.SetEntries(200000.0);
    const auto fitted = hpge::CalibrationEngine::FitRadwarePeak(
        histogram, expected.rangeHigh, expected.rangeLow);
    if (!fitted.success) {
        std::cerr << "RadWare peak fit failed: " << fitted.status << '\n';
        return false;
    }
    return Near(fitted.centroid, expected.centroid, 0.35, "RadWare centroid") &&
           Near(fitted.sigma, expected.sigma, 1.5, "RadWare sigma") &&
           fitted.centroidError > 0.0 && fitted.rangeLow == expected.rangeLow &&
           fitted.rangeHigh == expected.rangeHigh && fitted.ndf > 0;
}

bool TestRadwareFitValidation() {
    TH1D empty("radware_empty", "radware_empty", 100, 0.0, 100.0);
    const auto invalidRange = hpge::CalibrationEngine::FitRadwarePeak(empty, 20.0, 20.0);
    if (invalidRange.success || invalidRange.status.find("invalid") == std::string::npos) {
        std::cerr << "RadWare fit did not reject an invalid interval\n";
        return false;
    }
    const auto tooNarrow = hpge::CalibrationEngine::FitRadwarePeak(empty, 20.0, 25.0);
    if (tooNarrow.success || tooNarrow.status.find("12 bins") == std::string::npos) {
        std::cerr << "RadWare fit did not reject a narrow interval\n";
        return false;
    }
    const auto noCounts = hpge::CalibrationEngine::FitRadwarePeak(empty, 20.0, 50.0);
    if (noCounts.success || noCounts.status.find("no counts") == std::string::npos) {
        std::cerr << "RadWare fit did not reject an empty interval\n";
        return false;
    }
    return true;
}

bool TestMappedRadwareFit() {
    const std::vector<double> reference{600.0, 1200.0, 1850.0, 2700.0};
    const double scale = 1.08;
    const double offset = 37.0;
    auto target = MakeSpectrum("mapped_radware", reference, scale, offset);
    hpge::CalibrationEngine::SearchOptions options;
    const auto match = hpge::CalibrationEngine::MatchReferencePeaks(target, reference, options);
    if (!match.success) {
        std::cerr << "Could not map the target spectrum before interval fitting\n";
        return false;
    }
    for (std::size_t index = 0; index < reference.size(); ++index) {
        if (!match.matched[index]) {
            std::cerr << "Mapped interval is missing reference peak " << index << '\n';
            return false;
        }
        const double low = match.offset + match.scale * (reference[index] - 35.0);
        const double high = match.offset + match.scale * (reference[index] + 35.0);
        const auto fit = hpge::CalibrationEngine::FitRadwarePeak(target, low, high);
        if (!fit.success ||
            !Near(fit.centroid, offset + scale * reference[index], 1.0,
                  "mapped RadWare centroid " + std::to_string(index))) {
            std::cerr << "Mapped RadWare interval fit failed: " << fit.status << '\n';
            return false;
        }
    }
    return true;
}

bool TestAlignmentConstrainedCorrespondence() {
    const std::vector<double> pattern{420.0, 780.0, 1210.0, 1740.0, 2380.0, 3190.0};
    const std::vector<double> assigned{780.0, 1740.0, 3190.0};
    const double expectedScale = 1.07;
    const double expectedOffset = 29.0;
    const double expectedQuadratic = 1.4e-5;
    auto reference = MakeVariableIntensitySpectrum(
        "constrained_correspondence_reference", pattern,
        {950.0, 260.0, 780.0, 330.0, 620.0, 410.0}, 1.0, 0.0);
    // These three contaminants form a tempting identity-mapped copy of the
    // selected subset. The full six-peak pattern identifies the transformed
    // photopeaks and must override that locally plausible but wrong solution.
    auto target = MakeVariableIntensitySpectrum(
        "constrained_correspondence_target", pattern,
        {250.0, 980.0, 310.0, 760.0, 390.0, 680.0},
        expectedScale, expectedOffset, assigned, expectedQuadratic);
    hpge::CalibrationEngine::SearchOptions options;
    options.alignmentSensitivity = 0.55;
    const auto match = hpge::CalibrationEngine::FindCorrespondingPeaks(
        reference, target, assigned, options);
    if (!match.success) {
        std::cerr << "Alignment-constrained corresponding-peak search failed\n";
        return false;
    }
    for (std::size_t index = 0; index < assigned.size(); ++index) {
        const double expected = expectedOffset + expectedScale * assigned[index] +
                                expectedQuadratic * assigned[index] * assigned[index];
        if (!match.matched[index] ||
            !Near(match.charges[index], expected, 3.0,
                  "alignment-constrained corresponding peak")) {
            return false;
        }
        const auto fit = hpge::CalibrationEngine::FitRadwarePeak(
            target, match.charges[index] - 30.0, match.charges[index] + 30.0);
        if (!fit.success ||
            !Near(fit.centroid, expected, 1.0,
                  "aligned corresponding-peak fit")) {
            std::cerr << "Aligned candidate was not suitable for recentered fitting\n";
            return false;
        }
    }
    return true;
}

bool TestSpectrumAlignment() {
    const std::vector<double> peaks{470.0, 840.0, 1190.0, 1760.0, 2300.0, 3010.0};
    const double expectedScale = 1.075;
    const double expectedOffset = 31.0;
    const double expectedQuadratic = 1.2e-5;
    auto reference = MakeVariableIntensitySpectrum(
        "alignment_reference", peaks, {1100.0, 280.0, 850.0, 420.0, 730.0, 330.0},
        1.0, 0.0, {3650.0});
    auto target = MakeVariableIntensitySpectrum(
        "alignment_target", peaks, {260.0, 1050.0, 310.0, 900.0, 390.0, 780.0},
        expectedScale, expectedOffset, {330.0, 3550.0}, expectedQuadratic);
    hpge::CalibrationEngine::SearchOptions options;
    const auto match = hpge::CalibrationEngine::AlignSpectrumPatterns(reference, target, options);
    if (!match.success || !(match.scale > 0.0)) {
        std::cerr << "Pre-calibration spectrum alignment failed\n";
        return false;
    }
    const int matched = static_cast<int>(
        std::count(match.matched.begin(), match.matched.end(), true));
    return matched >= 5 &&
           Near(match.scale, expectedScale, 0.01, "intensity-independent alignment scale") &&
           Near(match.offset, expectedOffset, 8.0, "intensity-independent alignment offset") &&
           Near(match.quadratic, expectedQuadratic, 3e-6,
                "intensity-independent alignment curvature") &&
           Near(hpge::CalibrationEngine::MapTargetChargeToReference(
                    match, hpge::CalibrationEngine::MapReferenceCharge(match, 2500.0)),
                2500.0, 0.2, "alignment forward/inverse mapping");
}

bool TestTwoPeakPatternAlignment() {
    const std::vector<double> peaks{1450.0, 1690.0};
    auto reference = MakeVariableIntensitySpectrum(
        "co60_pattern_reference", peaks, {1200.0, 500.0}, 1.0, 0.0);
    auto target = MakeVariableIntensitySpectrum(
        "co60_pattern_target", peaks, {350.0, 1300.0}, 1.04, 22.0);
    hpge::CalibrationEngine::SearchOptions options;
    const auto match = hpge::CalibrationEngine::AlignSpectrumPatterns(reference, target, options);
    return match.success &&
           Near(match.scale, 1.04, 0.01, "two-peak pattern scale") &&
           Near(match.offset, 22.0, 8.0, "two-peak pattern offset");
}

bool TestCo56CrystalPatternAlignment() {
    const std::vector<double> energies{846.771, 1037.840, 1238.282,
                                       1771.351, 2598.459, 3201.962};
    std::vector<double> referenceCharges;
    std::vector<double> targetCharges;
    for (double energy : energies) {
        referenceCharges.push_back(ChargeForEnergy(energy, -1.5, 0.6800, 1.80e-5));
        targetCharges.push_back(ChargeForEnergy(energy, 0.9, 0.7240, 2.25e-5));
    }
    auto reference = MakeVariableIntensitySpectrum(
        "co56_crystal_reference", referenceCharges,
        {1200.0, 250.0, 900.0, 350.0, 700.0, 420.0}, 1.0, 0.0, {3800.0});
    auto target = MakeVariableIntensitySpectrum(
        "co56_crystal_target", targetCharges,
        {300.0, 1100.0, 270.0, 850.0, 400.0, 760.0}, 1.0, 0.0, {420.0});
    hpge::CalibrationEngine::SearchOptions options;
    const auto match = hpge::CalibrationEngine::AlignSpectrumPatterns(reference, target, options);
    if (!match.success) {
        std::cerr << "Co-56 crystal pattern alignment failed\n";
        return false;
    }
    bool accurate = true;
    for (std::size_t index = 0; index < referenceCharges.size(); ++index) {
        accurate &= Near(hpge::CalibrationEngine::MapReferenceCharge(
                             match, referenceCharges[index]),
                         targetCharges[index], 10.0,
                         "Co-56 mapped charge " + std::to_string(index));
    }
    if (!accurate) {
        std::cerr << "Co-56 alignment model="
                  << (match.quadraticModel ? "quadratic" : "affine")
                  << " candidates=" << match.referenceCharges.size()
                  << " coefficients=" << match.offset << ',' << match.scale << ','
                  << match.quadratic << '\n';
        for (std::size_t index = 0; index < match.referenceCharges.size(); ++index) {
            std::cerr << "  " << match.referenceCharges[index] << " -> "
                      << match.charges[index] << " matched=" << match.matched[index] << '\n';
        }
    }
    return accurate;
}

bool TestLowStatisticsSpikeAlignment() {
    const std::vector<double> peaks{590.0, 980.0, 1480.0, 2160.0, 2920.0};
    const double expectedScale = 1.047;
    const double expectedOffset = 24.0;
    auto reference = MakeLowStatisticsSpikySpectrum(
        "low_statistics_reference", peaks, {18.0, 10.0, 15.0, 11.0, 14.0}, 1.0, 0.0,
        {173, 418, 755, 1211, 1813, 2457, 3331, 3719});
    auto target = MakeLowStatisticsSpikySpectrum(
        "low_statistics_target", peaks, {8.0, 17.0, 9.0, 16.0, 10.0},
        expectedScale, expectedOffset,
        {95, 347, 681, 1137, 1679, 2533, 3187, 3971});

    hpge::CalibrationEngine::SearchOptions conservativeOptions;
    conservativeOptions.alignmentSensitivity = 0.15;
    conservativeOptions.autoTuneAlignmentSensitivity = false;
    const auto conservative = hpge::CalibrationEngine::AlignSpectrumPatterns(
        reference, target, conservativeOptions);
    if (!conservative.success) {
        std::cerr << "Conservative low-statistics alignment failed\n";
        return false;
    }
    const int matched = static_cast<int>(
        std::count(conservative.matched.begin(), conservative.matched.end(), true));
    if (conservative.referenceCharges.size() > peaks.size() + 1) {
        std::cerr << "Conservative sensitivity retained too many reference candidates: "
                  << conservative.referenceCharges.size() << '\n';
        return false;
    }
    for (double peak : peaks) {
        const bool retained = std::any_of(
            conservative.referenceCharges.begin(), conservative.referenceCharges.end(),
            [peak](double candidate) { return std::abs(candidate - peak) <= 3.0; });
        if (!retained) {
            std::cerr << "Conservative sensitivity rejected real low-count peak " << peak << '\n';
            return false;
        }
    }

    hpge::CalibrationEngine::SearchOptions highOptions;
    highOptions.alignmentSensitivity = 1.0;
    highOptions.autoTuneAlignmentSensitivity = false;
    const auto high = hpge::CalibrationEngine::AlignSpectrumPatterns(reference, target, highOptions);
    if (high.referenceCharges.size() <= conservative.referenceCharges.size()) {
        std::cerr << "Alignment sensitivity switch did not retain additional candidates\n";
        return false;
    }
    return matched >= 4 &&
           Near(conservative.scale, expectedScale, 0.012,
                "low-statistics alignment scale") &&
           Near(conservative.offset, expectedOffset, 9.0,
                "low-statistics alignment offset");
}

bool TestAdaptiveSpectrumSensitivity() {
    const std::vector<double> peaks{610.0, 1050.0, 1610.0, 2290.0, 3010.0};
    const double expectedScale = 1.035;
    const double expectedOffset = 18.0;
    auto reference = MakeLowStatisticsSpikySpectrum(
        "adaptive_low_count_reference", peaks, {16.0, 10.0, 14.0, 9.0, 12.0},
        1.0, 0.0, {217, 833, 1901, 2777, 3691});
    auto target = MakeLowStatisticsSpikySpectrum(
        "adaptive_high_count_target", peaks, {9.0, 17.0, 8.0, 15.0, 11.0},
        expectedScale, expectedOffset, {143, 719, 1327, 2591, 3529});
    for (int bin = 1; bin <= target.GetNbinsX(); ++bin) {
        target.SetBinContent(bin, 10.0 * target.GetBinContent(bin));
    }
    target.SetEntries(42000.0);

    hpge::CalibrationEngine::SearchOptions options;
    options.alignmentSensitivity = 0.35;
    options.autoTuneAlignmentSensitivity = true;
    options.alignmentModel = hpge::CalibrationEngine::AlignmentModel::Affine;
    const auto match = hpge::CalibrationEngine::AlignSpectrumPatterns(reference, target, options);
    return match.success && !match.quadraticModel &&
           match.referenceSensitivity > match.targetSensitivity + 0.08 &&
           Near(match.scale, expectedScale, 0.012, "per-spectrum tuned scale") &&
           Near(match.offset, expectedOffset, 9.0, "per-spectrum tuned offset");
}

bool TestQuadraticAlignmentWithIntenseXrays() {
    const std::vector<double> peaks{620.0, 1080.0, 1650.0, 2380.0, 3150.0};
    const double expectedScale = 1.025;
    const double expectedOffset = 21.0;
    const double expectedQuadratic = 3.2e-5;
    auto reference = MakeVariableIntensitySpectrum(
        "quadratic_xray_reference", peaks, {75.0, 42.0, 68.0, 38.0, 55.0},
        1.0, 0.0);
    auto target = MakeVariableIntensitySpectrum(
        "quadratic_xray_target", peaks, {35.0, 82.0, 31.0, 70.0, 40.0},
        expectedScale, expectedOffset, {}, expectedQuadratic);
    AddGaussianPeak(reference, 118.0, 22000.0, 3.0);
    AddGaussianPeak(reference, 162.0, 14000.0, 3.0);
    AddGaussianPeak(reference, 207.0, 18000.0, 3.0);
    AddGaussianPeak(target, 137.0, 26000.0, 3.0);
    AddGaussianPeak(target, 194.0, 17000.0, 3.0);

    hpge::CalibrationEngine::SearchOptions options;
    options.alignmentSensitivity = 0.5;
    options.autoTuneAlignmentSensitivity = true;
    options.alignmentModel = hpge::CalibrationEngine::AlignmentModel::Quadratic;
    const auto match = hpge::CalibrationEngine::AlignSpectrumPatterns(reference, target, options);
    if (!match.success || !match.quadraticModel) {
        std::cerr << "Explicit quadratic alignment failed with intense X-rays\n";
        return false;
    }
    auto affineOptions = options;
    affineOptions.alignmentModel = hpge::CalibrationEngine::AlignmentModel::Affine;
    const auto affine = hpge::CalibrationEngine::AlignSpectrumPatterns(
        reference, target, affineOptions);
    double affineWorstError = 0.0;
    for (double peak : peaks) {
        const double expected = expectedOffset + expectedScale * peak +
                                expectedQuadratic * peak * peak;
        affineWorstError = std::max(affineWorstError,
            std::abs(hpge::CalibrationEngine::MapReferenceCharge(affine, peak) - expected));
    }
    if (affine.quadraticModel || affineWorstError < 20.0) {
        std::cerr << "Quadratic regression is not distinct from affine alignment\n";
        return false;
    }
    const int matched = static_cast<int>(
        std::count(match.matched.begin(), match.matched.end(), true));
    bool mapped = matched >= 4 &&
        Near(match.scale, expectedScale, 0.02, "quadratic X-ray scale") &&
        Near(match.offset, expectedOffset, 15.0, "quadratic X-ray offset") &&
        Near(match.quadratic, expectedQuadratic, 6e-6, "quadratic X-ray curvature");
    for (double peak : peaks) {
        mapped &= Near(hpge::CalibrationEngine::MapReferenceCharge(match, peak),
                       expectedOffset + expectedScale * peak +
                           expectedQuadratic * peak * peak,
                       9.0, "quadratic X-ray mapped peak");
    }
    if (!mapped) {
        std::cerr << "Quadratic X-ray coefficients=" << match.offset << ',' << match.scale
                  << ',' << match.quadratic << " candidates="
                  << match.referenceCharges.size() << '\n';
        for (std::size_t index = 0; index < match.referenceCharges.size(); ++index) {
            std::cerr << "  " << match.referenceCharges[index] << " -> "
                      << match.charges[index] << " matched=" << match.matched[index] << '\n';
        }
    }
    return mapped;
}

bool TestPatternAlignmentValidation() {
    TH1D emptyReference("empty_pattern_reference", "empty", 100, 0.0, 100.0);
    TH1D emptyTarget("empty_pattern_target", "empty", 100, 0.0, 100.0);
    hpge::CalibrationEngine::SearchOptions options;
    const auto empty = hpge::CalibrationEngine::AlignSpectrumPatterns(
        emptyReference, emptyTarget, options);
    if (empty.success || !empty.referenceCharges.empty()) {
        std::cerr << "Empty spectra were accepted for pattern alignment\n";
        return false;
    }
    auto oneReference = MakeSpectrum("one_pattern_reference", {1000.0});
    auto oneTarget = MakeSpectrum("one_pattern_target", {1000.0}, 1.1, 20.0);
    const auto single = hpge::CalibrationEngine::AlignSpectrumPatterns(
        oneReference, oneTarget, options);
    if (single.success) {
        std::cerr << "Single peaks were accepted as a spectrum pattern\n";
        return false;
    }
    hpge::PeakMatchResult linear;
    linear.offset = 17.0;
    linear.scale = 1.06;
    const double mapped = hpge::CalibrationEngine::MapReferenceCharge(linear, 1234.0);
    return Near(hpge::CalibrationEngine::MapTargetChargeToReference(linear, mapped),
                1234.0, 1e-9, "linear alignment forward/inverse mapping");
}

bool TestFitQuadratic() {
    const double p0 = -2.3;
    const double p1 = 0.71;
    const double p2 = 2.1e-5;
    std::vector<hpge::CalibrationPoint> points;
    for (double charge : {550.0, 1050.0, 1700.0, 2400.0, 3300.0}) {
        const double energy = p0 + p1 * charge + p2 * charge * charge;
        points.push_back(MakeCalibrationPoint(charge, energy, "synthetic"));
    }
    const auto fit = hpge::CalibrationEngine::FitSecondOrder(7, points, 0.1);
    if (!fit.success || fit.crystal != 7 || fit.needsReview || fit.ndf != 2) {
        std::cerr << "Exact quadratic fit has incorrect status: " << fit.status << '\n';
        return false;
    }
    return Near(fit.p0, p0, 1e-7, "p0") &&
           Near(fit.p1, p1, 1e-10, "p1") &&
           Near(fit.p2, p2, 1e-12, "p2") &&
           Near(fit.residualRms, 0.0, 1e-8, "residual RMS");
}

bool TestFitValidation() {
    std::vector<hpge::CalibrationPoint> twoPoints{
        MakeCalibrationPoint(100.0, 70.0), MakeCalibrationPoint(200.0, 140.0)};
    const auto insufficient = hpge::CalibrationEngine::FitSecondOrder(0, twoPoints, 1.0);
    if (insufficient.success || insufficient.status.find("fewer than 3") == std::string::npos) {
        std::cerr << "Fit did not reject fewer than three points\n";
        return false;
    }

    std::vector<hpge::CalibrationPoint> duplicateCharges{
        MakeCalibrationPoint(100.0, 60.0), MakeCalibrationPoint(100.0, 70.0),
        MakeCalibrationPoint(100.0, 80.0)};
    const auto singular = hpge::CalibrationEngine::FitSecondOrder(0, duplicateCharges, 1.0);
    if (singular.success || singular.status.find("not distinct") == std::string::npos) {
        std::cerr << "Fit did not reject identical charges\n";
        return false;
    }

    std::vector<hpge::CalibrationPoint> exactlyThree{
        MakeCalibrationPoint(100.0, 70.0), MakeCalibrationPoint(200.0, 142.0),
        MakeCalibrationPoint(300.0, 216.0)};
    const auto noRedundancy = hpge::CalibrationEngine::FitSecondOrder(0, exactlyThree, 1.0);
    if (!noRedundancy.success || !noRedundancy.needsReview || noRedundancy.ndf != 0) {
        std::cerr << "Exactly determined quadratic was not flagged for review\n";
        return false;
    }

    std::vector<hpge::CalibrationPoint> noisy{
        MakeCalibrationPoint(100.0, 70.0), MakeCalibrationPoint(200.0, 160.0),
        MakeCalibrationPoint(300.0, 205.0), MakeCalibrationPoint(400.0, 330.0)};
    const auto review = hpge::CalibrationEngine::FitSecondOrder(0, noisy, 0.1);
    if (!review.success || !review.needsReview || review.residualRms <= 0.1) {
        std::cerr << "High-residual fit was not flagged for review\n";
        return false;
    }
    return true;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Expected one test-case name\n";
        return 2;
    }
    const std::string test = argv[1];
    bool passed = false;
    if (test == "find-peaks") passed = TestFindPeaks();
    else if (test == "refine-peak") passed = TestRefinePeak();
    else if (test == "match-multiple-peaks") passed = TestMatchMultiplePeaks();
    else if (test == "match-single-peak") passed = TestMatchSinglePeak();
    else if (test == "radware-peak-fit") passed = TestRadwarePeakFit();
    else if (test == "radware-fit-validation") passed = TestRadwareFitValidation();
    else if (test == "mapped-radware-fit") passed = TestMappedRadwareFit();
    else if (test == "alignment-constrained-correspondence") passed = TestAlignmentConstrainedCorrespondence();
    else if (test == "spectrum-alignment") passed = TestSpectrumAlignment();
    else if (test == "two-peak-pattern-alignment") passed = TestTwoPeakPatternAlignment();
    else if (test == "co56-crystal-pattern-alignment") passed = TestCo56CrystalPatternAlignment();
    else if (test == "low-statistics-spike-alignment") passed = TestLowStatisticsSpikeAlignment();
    else if (test == "adaptive-spectrum-sensitivity") passed = TestAdaptiveSpectrumSensitivity();
    else if (test == "quadratic-xray-alignment") passed = TestQuadraticAlignmentWithIntenseXrays();
    else if (test == "pattern-alignment-validation") passed = TestPatternAlignmentValidation();
    else if (test == "fit-quadratic") passed = TestFitQuadratic();
    else if (test == "fit-validation") passed = TestFitValidation();
    else {
        std::cerr << "Unknown test case: " << test << '\n';
        return 2;
    }
    if (passed) std::cout << "PASS: " << test << '\n';
    return passed ? 0 : 1;
}
