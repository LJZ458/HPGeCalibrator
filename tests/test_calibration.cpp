#include "CalibrationEngine.h"

#include <TH1D.h>

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

bool TestFitQuadratic() {
    const double p0 = -2.3;
    const double p1 = 0.71;
    const double p2 = 2.1e-5;
    std::vector<hpge::CalibrationPoint> points;
    for (double charge : {550.0, 1050.0, 1700.0, 2400.0, 3300.0}) {
        const double energy = p0 + p1 * charge + p2 * charge * charge;
        points.push_back({"synthetic", charge, energy, 0.0, false, 0.0});
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
        {"test", 100.0, 70.0}, {"test", 200.0, 140.0}};
    const auto insufficient = hpge::CalibrationEngine::FitSecondOrder(0, twoPoints, 1.0);
    if (insufficient.success || insufficient.status.find("fewer than 3") == std::string::npos) {
        std::cerr << "Fit did not reject fewer than three points\n";
        return false;
    }

    std::vector<hpge::CalibrationPoint> duplicateCharges{
        {"test", 100.0, 60.0}, {"test", 100.0, 70.0}, {"test", 100.0, 80.0}};
    const auto singular = hpge::CalibrationEngine::FitSecondOrder(0, duplicateCharges, 1.0);
    if (singular.success || singular.status.find("not distinct") == std::string::npos) {
        std::cerr << "Fit did not reject identical charges\n";
        return false;
    }

    std::vector<hpge::CalibrationPoint> exactlyThree{
        {"test", 100.0, 70.0}, {"test", 200.0, 142.0}, {"test", 300.0, 216.0}};
    const auto noRedundancy = hpge::CalibrationEngine::FitSecondOrder(0, exactlyThree, 1.0);
    if (!noRedundancy.success || !noRedundancy.needsReview || noRedundancy.ndf != 0) {
        std::cerr << "Exactly determined quadratic was not flagged for review\n";
        return false;
    }

    std::vector<hpge::CalibrationPoint> noisy{
        {"test", 100.0, 70.0}, {"test", 200.0, 160.0},
        {"test", 300.0, 205.0}, {"test", 400.0, 330.0}};
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
    else if (test == "fit-quadratic") passed = TestFitQuadratic();
    else if (test == "fit-validation") passed = TestFitValidation();
    else {
        std::cerr << "Unknown test case: " << test << '\n';
        return 2;
    }
    if (passed) std::cout << "PASS: " << test << '\n';
    return passed ? 0 : 1;
}

