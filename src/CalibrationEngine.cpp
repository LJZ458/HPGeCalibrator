#include "CalibrationEngine.h"

#include <TAxis.h>
#include <TH1.h>
#include <TSpectrum.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>

namespace hpge {
namespace {

bool SolveThreeByThree(double matrix[3][4], double solution[3]) {
    for (int column = 0; column < 3; ++column) {
        int pivot = column;
        for (int row = column + 1; row < 3; ++row) {
            if (std::abs(matrix[row][column]) > std::abs(matrix[pivot][column])) pivot = row;
        }
        if (std::abs(matrix[pivot][column]) < 1e-14) return false;
        if (pivot != column) {
            for (int item = column; item < 4; ++item) {
                std::swap(matrix[column][item], matrix[pivot][item]);
            }
        }
        const double divisor = matrix[column][column];
        for (int item = column; item < 4; ++item) matrix[column][item] /= divisor;
        for (int row = 0; row < 3; ++row) {
            if (row == column) continue;
            const double factor = matrix[row][column];
            for (int item = column; item < 4; ++item) {
                matrix[row][item] -= factor * matrix[column][item];
            }
        }
    }
    for (int i = 0; i < 3; ++i) solution[i] = matrix[i][3];
    return true;
}

double NearestCandidate(const std::vector<double>& candidates, double expected,
                        double tolerance, bool& found) {
    found = false;
    double best = 0.0;
    double distance = tolerance;
    const auto it = std::lower_bound(candidates.begin(), candidates.end(), expected);
    if (it != candidates.end() && std::abs(*it - expected) <= distance) {
        best = *it;
        distance = std::abs(*it - expected);
        found = true;
    }
    if (it != candidates.begin()) {
        const double candidate = *std::prev(it);
        if (std::abs(candidate - expected) <= distance) {
            best = candidate;
            found = true;
        }
    }
    return best;
}

} // namespace

std::vector<double> CalibrationEngine::FindPeakCandidates(
    const TH1& histogram, const SearchOptions& options) {
    auto copy = std::unique_ptr<TH1>(dynamic_cast<TH1*>(histogram.Clone()));
    if (!copy || copy->GetEntries() <= 0) return {};
    copy->SetDirectory(nullptr);

    TSpectrum search(std::max(2, options.maxPeaks));
    const int count = search.Search(copy.get(), std::max(0.5, options.sigmaBins),
                                    "goff nodraw", std::clamp(options.threshold, 0.001, 0.99));
    std::vector<double> peaks;
    peaks.reserve(count);
    const double* positions = search.GetPositionX();
    for (int i = 0; i < count; ++i) peaks.push_back(positions[i]);
    std::sort(peaks.begin(), peaks.end());
    return peaks;
}

double CalibrationEngine::RefinePeak(const TH1& histogram, double approximateCharge,
                                     double halfWindow) {
    const TAxis* axis = histogram.GetXaxis();
    const int first = std::max(1, axis->FindBin(approximateCharge - halfWindow));
    const int last = std::min(histogram.GetNbinsX(),
                              axis->FindBin(approximateCharge + halfWindow));
    if (first > last) return approximateCharge;

    int maximumBin = first;
    for (int bin = first + 1; bin <= last; ++bin) {
        if (histogram.GetBinContent(bin) > histogram.GetBinContent(maximumBin)) {
            maximumBin = bin;
        }
    }

    // A three-bin parabolic interpolation avoids a fit failure on sparse spectra.
    if (maximumBin > 1 && maximumBin < histogram.GetNbinsX()) {
        const double left = histogram.GetBinContent(maximumBin - 1);
        const double middle = histogram.GetBinContent(maximumBin);
        const double right = histogram.GetBinContent(maximumBin + 1);
        const double denominator = left - 2.0 * middle + right;
        if (std::abs(denominator) > 1e-12) {
            const double shift = std::clamp(0.5 * (left - right) / denominator, -1.0, 1.0);
            return axis->GetBinCenter(maximumBin) + shift * axis->GetBinWidth(maximumBin);
        }
    }
    return axis->GetBinCenter(maximumBin);
}

PeakMatchResult CalibrationEngine::MatchReferencePeaks(
    const TH1& target, const std::vector<double>& referenceCharges,
    const SearchOptions& options) {
    PeakMatchResult result;
    result.charges.assign(referenceCharges.size(), 0.0);
    result.matched.assign(referenceCharges.size(), false);
    if (referenceCharges.empty()) return result;

    auto candidates = FindPeakCandidates(target, options);
    if (candidates.empty()) return result;

    const double targetRange = target.GetXaxis()->GetXmax() - target.GetXaxis()->GetXmin();
    const double tolerance = std::max(2.5 * target.GetXaxis()->GetBinWidth(1),
                                      options.matchToleranceFraction * targetRange);

    if (referenceCharges.size() == 1) {
        bool found = false;
        const double matched = NearestCandidate(candidates, referenceCharges.front(),
                                                4.0 * tolerance, found);
        if (found) {
            result.charges[0] = matched;
            result.matched[0] = true;
            result.offset = matched - referenceCharges.front();
            result.score = 1.0;
            result.success = true;
        }
        return result;
    }

    std::vector<std::size_t> order(referenceCharges.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return referenceCharges[a] < referenceCharges[b];
    });
    const double refLow = referenceCharges[order.front()];
    const double refHigh = referenceCharges[order.back()];
    if (refHigh <= refLow) return result;

    double bestScore = -std::numeric_limits<double>::infinity();
    double bestScale = 1.0;
    double bestOffset = 0.0;
    std::vector<double> bestCharges(referenceCharges.size(), 0.0);
    std::vector<bool> bestMatched(referenceCharges.size(), false);

    for (std::size_t i = 0; i + 1 < candidates.size(); ++i) {
        for (std::size_t j = i + 1; j < candidates.size(); ++j) {
            const double scale = (candidates[j] - candidates[i]) / (refHigh - refLow);
            if (scale < 0.2 || scale > 5.0) continue;
            const double offset = candidates[i] - scale * refLow;
            std::vector<double> trialCharges(referenceCharges.size(), 0.0);
            std::vector<bool> trialMatched(referenceCharges.size(), false);
            int count = 0;
            double normalizedError = 0.0;
            for (std::size_t r = 0; r < referenceCharges.size(); ++r) {
                const double expected = offset + scale * referenceCharges[r];
                bool found = false;
                const double match = NearestCandidate(candidates, expected, tolerance, found);
                if (found) {
                    trialCharges[r] = match;
                    trialMatched[r] = true;
                    ++count;
                    normalizedError += std::abs(match - expected) / tolerance;
                }
            }
            // Count dominates; error and extreme scale changes break ties.
            const double score = 100.0 * count - normalizedError - 0.5 * std::abs(std::log(scale));
            if (score > bestScore) {
                bestScore = score;
                bestScale = scale;
                bestOffset = offset;
                bestCharges = std::move(trialCharges);
                bestMatched = std::move(trialMatched);
            }
        }
    }

    int matchedCount = 0;
    for (bool value : bestMatched) matchedCount += value ? 1 : 0;
    result.success = matchedCount >= 2;
    result.charges = std::move(bestCharges);
    result.matched = std::move(bestMatched);
    result.scale = bestScale;
    result.offset = bestOffset;
    result.score = bestScore;
    return result;
}

CalibrationResult CalibrationEngine::FitSecondOrder(
    int crystal, std::vector<CalibrationPoint> points, double residualRmsLimitKeV) {
    CalibrationResult result;
    result.crystal = crystal;
    result.points = std::move(points);
    if (result.points.size() < 3) {
        result.status = "fewer than 3 calibration peaks";
        return result;
    }

    std::sort(result.points.begin(), result.points.end(),
              [](const CalibrationPoint& a, const CalibrationPoint& b) {
                  return a.charge < b.charge;
              });
    const double meanCharge = std::accumulate(
        result.points.begin(), result.points.end(), 0.0,
        [](double sum, const CalibrationPoint& point) { return sum + point.charge; }) /
        static_cast<double>(result.points.size());
    double chargeScale = 0.0;
    for (const auto& point : result.points) {
        chargeScale = std::max(chargeScale, std::abs(point.charge - meanCharge));
    }
    if (chargeScale <= 0.0) {
        result.status = "calibration charges are not distinct";
        return result;
    }

    // Fit in z=(charge-mean)/scale coordinates to keep the normal equations
    // well conditioned even for large ADC values.
    double sums[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double rhs[3] = {0.0, 0.0, 0.0};
    for (const auto& point : result.points) {
        const double z = (point.charge - meanCharge) / chargeScale;
        double power = 1.0;
        for (double& sum : sums) {
            sum += power;
            power *= z;
        }
        rhs[0] += point.energy;
        rhs[1] += point.energy * z;
        rhs[2] += point.energy * z * z;
    }
    double system[3][4] = {
        {sums[0], sums[1], sums[2], rhs[0]},
        {sums[1], sums[2], sums[3], rhs[1]},
        {sums[2], sums[3], sums[4], rhs[2]}
    };
    double scaled[3] = {0.0, 0.0, 0.0};
    if (!SolveThreeByThree(system, scaled)) {
        result.status = "quadratic least-squares system is singular";
        return result;
    }
    result.p2 = scaled[2] / (chargeScale * chargeScale);
    result.p1 = scaled[1] / chargeScale - 2.0 * meanCharge * result.p2;
    result.p0 = scaled[0] - scaled[1] * meanCharge / chargeScale +
                scaled[2] * meanCharge * meanCharge / (chargeScale * chargeScale);
    result.ndf = static_cast<int>(result.points.size()) - 3;
    double sumSquares = 0.0;
    for (auto& point : result.points) {
        const double fitted = result.p0 + result.p1 * point.charge +
                              result.p2 * point.charge * point.charge;
        point.residual = point.energy - fitted;
        sumSquares += point.residual * point.residual;
    }
    result.chi2 = sumSquares;
    result.residualRms = std::sqrt(sumSquares / static_cast<double>(result.points.size()));
    result.success = true;
    result.needsReview = result.residualRms > residualRmsLimitKeV || result.ndf == 0;
    if (result.ndf == 0) {
        result.status = "fit complete; exactly 3 peaks (no fit redundancy)";
    } else if (result.needsReview) {
        result.status = "fit complete; residual RMS exceeds review limit";
    } else {
        result.status = "fit complete";
    }
    return result;
}

} // namespace hpge
