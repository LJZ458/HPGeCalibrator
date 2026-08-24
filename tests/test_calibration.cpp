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
    return accurate;
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
    else if (test == "spectrum-alignment") passed = TestSpectrumAlignment();
    else if (test == "two-peak-pattern-alignment") passed = TestTwoPeakPatternAlignment();
    else if (test == "co56-crystal-pattern-alignment") passed = TestCo56CrystalPatternAlignment();
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
