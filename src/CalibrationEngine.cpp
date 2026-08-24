#include "CalibrationEngine.h"

#include <TAxis.h>
#include <TH1.h>
#include <TSpectrum.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>

namespace hpge {
namespace {

constexpr std::size_t kRadwareParameters = 9;
constexpr double kPi = 3.14159265358979323846;
using RadwareParameters = std::array<double, kRadwareParameters>;

enum RadwareParameter : std::size_t {
    kHeight,
    kCentroid,
    kSigma,
    kTailFraction,
    kBeta,
    kStepFraction,
    kBackground0,
    kBackground1,
    kBackground2
};

double RadwareModel(double charge, const RadwareParameters& parameters,
                    double rangeLow, double rangeHigh) {
    const double sigma = std::max(parameters[kSigma], 1e-9);
    const double beta = std::max(parameters[kBeta], 1e-9);
    const double t = (charge - parameters[kCentroid]) / sigma;
    const double gaussian = parameters[kHeight] * (1.0 - parameters[kTailFraction]) *
                            std::exp(-0.5 * t * t);

    // GF3/RadWare low-energy tail: an exponential convolved with the Gaussian.
    const double k = sigma / beta;
    const double exponent = std::clamp(0.5 * k * k +
                                      (charge - parameters[kCentroid]) / beta,
                                      -700.0, 700.0);
    const double tailShape = std::sqrt(2.0 * kPi) * sigma / (2.0 * beta) *
                             std::exp(exponent) *
                             std::erfc((t + k) / std::sqrt(2.0));
    const double tail = parameters[kHeight] * parameters[kTailFraction] * tailShape;
    const double step = parameters[kHeight] * parameters[kStepFraction] * 0.5 *
                        std::erfc(t / std::sqrt(2.0));

    const double halfRange = std::max(0.5 * (rangeHigh - rangeLow), 1e-9);
    const double z = (charge - 0.5 * (rangeLow + rangeHigh)) / halfRange;
    const double background = parameters[kBackground0] + parameters[kBackground1] * z +
                              parameters[kBackground2] * z * z;
    return gaussian + tail + step + background;
}

void ClampRadwareParameters(RadwareParameters& parameters, double rangeLow,
                            double rangeHigh, double binWidth, double maximumCount) {
    const double width = rangeHigh - rangeLow;
    parameters[kHeight] = std::clamp(parameters[kHeight], 0.0,
                                     std::max(10.0 * maximumCount, 1.0));
    parameters[kCentroid] = std::clamp(parameters[kCentroid], rangeLow + 0.25 * binWidth,
                                       rangeHigh - 0.25 * binWidth);
    parameters[kSigma] = std::clamp(parameters[kSigma], 0.35 * binWidth,
                                    std::max(width / 3.0, 0.36 * binWidth));
    parameters[kTailFraction] = std::clamp(parameters[kTailFraction], 0.0, 0.4);
    parameters[kBeta] = std::clamp(parameters[kBeta],
                                   std::max(0.25 * parameters[kSigma], 0.25 * binWidth),
                                   std::max(width, 0.26 * binWidth));
    parameters[kStepFraction] = std::clamp(parameters[kStepFraction], 0.0, 0.25);
    for (std::size_t index = kBackground0; index <= kBackground2; ++index) {
        parameters[index] = std::clamp(parameters[index], -2.0 * maximumCount,
                                       2.0 * maximumCount);
    }
}

template <std::size_t N>
bool SolveLinearSystem(std::array<std::array<double, N + 1>, N> matrix,
                       std::array<double, N>& solution) {
    for (std::size_t column = 0; column < N; ++column) {
        std::size_t pivot = column;
        for (std::size_t row = column + 1; row < N; ++row) {
            if (std::abs(matrix[row][column]) > std::abs(matrix[pivot][column])) pivot = row;
        }
        if (!std::isfinite(matrix[pivot][column]) ||
            std::abs(matrix[pivot][column]) < 1e-16) return false;
        if (pivot != column) std::swap(matrix[pivot], matrix[column]);
        const double divisor = matrix[column][column];
        for (std::size_t item = column; item <= N; ++item) matrix[column][item] /= divisor;
        for (std::size_t row = 0; row < N; ++row) {
            if (row == column) continue;
            const double factor = matrix[row][column];
            for (std::size_t item = column; item <= N; ++item) {
                matrix[row][item] -= factor * matrix[column][item];
            }
        }
    }
    for (std::size_t index = 0; index < N; ++index) solution[index] = matrix[index][N];
    return true;
}

double RadwareChiSquare(const std::vector<double>& x, const std::vector<double>& y,
                        const RadwareParameters& parameters, double rangeLow,
                        double rangeHigh) {
    double chi2 = 0.0;
    for (std::size_t bin = 0; bin < x.size(); ++bin) {
        const double variance = std::max(y[bin], 1.0);
        const double residual = y[bin] - RadwareModel(x[bin], parameters, rangeLow, rangeHigh);
        chi2 += residual * residual / variance;
    }
    return chi2;
}

void RadwareNormalEquations(
    const std::vector<double>& x, const std::vector<double>& y,
    const RadwareParameters& parameters, double rangeLow, double rangeHigh,
    double binWidth, double maximumCount,
    std::array<std::array<double, kRadwareParameters>, kRadwareParameters>& normal,
    RadwareParameters& rightHandSide) {
    for (auto& row : normal) row.fill(0.0);
    rightHandSide.fill(0.0);
    const std::array<double, kRadwareParameters> naturalScale{
        std::max(maximumCount, 1.0), std::max(rangeHigh - rangeLow, binWidth),
        std::max(rangeHigh - rangeLow, binWidth), 1.0,
        std::max(rangeHigh - rangeLow, binWidth), 1.0,
        std::max(maximumCount, 1.0), std::max(maximumCount, 1.0),
        std::max(maximumCount, 1.0)};

    for (std::size_t bin = 0; bin < x.size(); ++bin) {
        const double variance = std::max(y[bin], 1.0);
        const double weight = 1.0 / std::sqrt(variance);
        const double model = RadwareModel(x[bin], parameters, rangeLow, rangeHigh);
        const double residual = (y[bin] - model) * weight;
        RadwareParameters derivative{};
        for (std::size_t parameter = 0; parameter < kRadwareParameters; ++parameter) {
            const double delta = 1e-4 *
                std::max(std::abs(parameters[parameter]), naturalScale[parameter] * 1e-3);
            auto plus = parameters;
            auto minus = parameters;
            plus[parameter] += delta;
            minus[parameter] -= delta;
            ClampRadwareParameters(plus, rangeLow, rangeHigh, binWidth, maximumCount);
            ClampRadwareParameters(minus, rangeLow, rangeHigh, binWidth, maximumCount);
            const double span = plus[parameter] - minus[parameter];
            derivative[parameter] = std::abs(span) > 1e-15
                ? (RadwareModel(x[bin], plus, rangeLow, rangeHigh) -
                   RadwareModel(x[bin], minus, rangeLow, rangeHigh)) / span * weight
                : 0.0;
        }
        for (std::size_t row = 0; row < kRadwareParameters; ++row) {
            rightHandSide[row] += derivative[row] * residual;
            for (std::size_t column = 0; column < kRadwareParameters; ++column) {
                normal[row][column] += derivative[row] * derivative[column];
            }
        }
    }
}

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

PeakFitResult CalibrationEngine::FitRadwarePeak(const TH1& histogram, double rangeLow,
                                                double rangeHigh) {
    PeakFitResult result;
    if (rangeHigh < rangeLow) std::swap(rangeLow, rangeHigh);
    result.rangeLow = rangeLow;
    result.rangeHigh = rangeHigh;
    const TAxis* axis = histogram.GetXaxis();
    rangeLow = std::max(rangeLow, axis->GetXmin());
    rangeHigh = std::min(rangeHigh, axis->GetXmax());
    result.rangeLow = rangeLow;
    result.rangeHigh = rangeHigh;
    if (!(rangeHigh > rangeLow)) {
        result.status = "invalid peak-fit interval";
        return result;
    }

    const int firstBin = std::max(1, axis->FindBin(rangeLow));
    const int lastBin = std::min(histogram.GetNbinsX(), axis->FindBin(rangeHigh));
    if (lastBin - firstBin + 1 < 12) {
        result.status = "peak-fit interval must contain at least 12 bins";
        return result;
    }
    std::vector<double> x;
    std::vector<double> y;
    x.reserve(lastBin - firstBin + 1);
    y.reserve(lastBin - firstBin + 1);
    double maximumCount = 0.0;
    int maximumIndex = 0;
    for (int bin = firstBin; bin <= lastBin; ++bin) {
        x.push_back(axis->GetBinCenter(bin));
        y.push_back(std::max(histogram.GetBinContent(bin), 0.0));
        if (y.back() > maximumCount) {
            maximumCount = y.back();
            maximumIndex = static_cast<int>(y.size()) - 1;
        }
    }
    if (maximumCount <= 0.0) {
        result.status = "peak-fit interval contains no counts";
        return result;
    }

    const std::size_t edgeBins = std::max<std::size_t>(2, x.size() / 8);
    const double leftBackground = std::accumulate(y.begin(), y.begin() + edgeBins, 0.0) /
                                  static_cast<double>(edgeBins);
    const double rightBackground = std::accumulate(y.end() - edgeBins, y.end(), 0.0) /
                                   static_cast<double>(edgeBins);
    const double binWidth = axis->GetBinWidth(firstBin);
    const double width = rangeHigh - rangeLow;
    double weightedSum = 0.0;
    double weightedSquare = 0.0;
    for (std::size_t bin = 0; bin < x.size(); ++bin) {
        const double z = (x[bin] - rangeLow) / width;
        const double background = leftBackground + z * (rightBackground - leftBackground);
        const double signal = std::max(y[bin] - background, 0.0);
        weightedSum += signal;
        weightedSquare += signal * std::pow(x[bin] - x[maximumIndex], 2);
    }
    const double initialSigma = weightedSum > 0.0
        ? std::sqrt(weightedSquare / weightedSum)
        : width / 10.0;

    RadwareParameters parameters{
        std::max(maximumCount - 0.5 * (leftBackground + rightBackground), 1.0),
        x[maximumIndex],
        std::clamp(initialSigma, binWidth, width / 5.0),
        0.05,
        std::max(3.0 * std::clamp(initialSigma, binWidth, width / 5.0), binWidth),
        0.01,
        0.5 * (leftBackground + rightBackground),
        0.5 * (rightBackground - leftBackground),
        0.0};
    ClampRadwareParameters(parameters, rangeLow, rangeHigh, binWidth, maximumCount);

    double chi2 = RadwareChiSquare(x, y, parameters, rangeLow, rangeHigh);
    double damping = 1e-2;
    bool improved = false;
    for (int iteration = 0; iteration < 120; ++iteration) {
        std::array<std::array<double, kRadwareParameters>, kRadwareParameters> normal{};
        RadwareParameters rightHandSide{};
        RadwareNormalEquations(x, y, parameters, rangeLow, rangeHigh, binWidth,
                               maximumCount, normal, rightHandSide);
        std::array<std::array<double, kRadwareParameters + 1>, kRadwareParameters> system{};
        for (std::size_t row = 0; row < kRadwareParameters; ++row) {
            for (std::size_t column = 0; column < kRadwareParameters; ++column) {
                system[row][column] = normal[row][column];
            }
            system[row][row] += damping * (normal[row][row] + 1e-9);
            system[row][kRadwareParameters] = rightHandSide[row];
        }
        RadwareParameters change{};
        if (!SolveLinearSystem<kRadwareParameters>(system, change)) {
            damping *= 10.0;
            continue;
        }
        auto candidate = parameters;
        for (std::size_t index = 0; index < kRadwareParameters; ++index) {
            candidate[index] += change[index];
        }
        ClampRadwareParameters(candidate, rangeLow, rangeHigh, binWidth, maximumCount);
        const double candidateChi2 = RadwareChiSquare(x, y, candidate, rangeLow, rangeHigh);
        if (std::isfinite(candidateChi2) && candidateChi2 < chi2) {
            const double relativeImprovement = (chi2 - candidateChi2) / std::max(chi2, 1.0);
            parameters = candidate;
            chi2 = candidateChi2;
            damping = std::max(damping / 3.0, 1e-9);
            improved = true;
            if (relativeImprovement < 1e-9 &&
                std::abs(change[kCentroid]) < 1e-5 * binWidth) break;
        } else {
            damping = std::min(damping * 10.0, 1e12);
        }
    }

    result.centroid = parameters[kCentroid];
    result.sigma = parameters[kSigma];
    result.height = parameters[kHeight];
    result.tailFraction = parameters[kTailFraction];
    result.beta = parameters[kBeta];
    result.stepFraction = parameters[kStepFraction];
    result.background0 = parameters[kBackground0];
    result.background1 = parameters[kBackground1];
    result.background2 = parameters[kBackground2];
    result.chi2 = chi2;
    result.ndf = static_cast<int>(x.size()) - static_cast<int>(kRadwareParameters);

    std::array<std::array<double, kRadwareParameters>, kRadwareParameters> normal{};
    RadwareParameters rhs{};
    RadwareNormalEquations(x, y, parameters, rangeLow, rangeHigh, binWidth,
                           maximumCount, normal, rhs);
    std::array<std::array<double, kRadwareParameters + 1>, kRadwareParameters> covarianceSystem{};
    for (std::size_t row = 0; row < kRadwareParameters; ++row) {
        for (std::size_t column = 0; column < kRadwareParameters; ++column) {
            covarianceSystem[row][column] = normal[row][column];
        }
        covarianceSystem[row][kRadwareParameters] = row == kCentroid ? 1.0 : 0.0;
    }
    RadwareParameters covarianceColumn{};
    if (SolveLinearSystem<kRadwareParameters>(covarianceSystem, covarianceColumn) &&
        covarianceColumn[kCentroid] > 0.0) {
        const double reducedChi2 = chi2 / std::max(result.ndf, 1);
        result.centroidError = std::sqrt(covarianceColumn[kCentroid] *
                                         std::max(reducedChi2, 1.0));
    } else {
        result.centroidError = result.sigma /
            std::sqrt(std::max(result.height, 1.0));
    }
    result.success = improved && result.height > 0.0 && result.sigma > 0.0 &&
                     result.centroid > rangeLow && result.centroid < rangeHigh;
    result.status = result.success ? "RadWare peak fit complete" :
                                     "RadWare peak fit did not converge";
    return result;
}

double CalibrationEngine::EvaluateRadwarePeak(double charge, const PeakFitResult& fit) {
    RadwareParameters parameters{
        fit.height, fit.centroid, fit.sigma, fit.tailFraction, fit.beta,
        fit.stepFraction, fit.background0, fit.background1, fit.background2};
    return RadwareModel(charge, parameters, fit.rangeLow, fit.rangeHigh);
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
