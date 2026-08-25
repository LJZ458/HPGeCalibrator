#pragma once

#include "CalibrationTypes.h"

#include <vector>

class TH1;

namespace hpge {

class CalibrationEngine {
public:
    enum class AlignmentModel {
        Auto,
        Affine,
        Quadratic
    };

    struct SearchOptions {
        double sigmaBins = 2.0;
        double threshold = 0.05;
        int maxPeaks = 80;
        double matchToleranceFraction = 0.018;
        double alignmentSensitivity = 0.35;
        bool autoTuneAlignmentSensitivity = true;
        AlignmentModel alignmentModel = AlignmentModel::Auto;
        bool suppressLowEnergyAlignmentCandidates = true;
    };

    static std::vector<double> FindPeakCandidates(const TH1& histogram,
                                                  const SearchOptions& options);

    static double RefinePeak(const TH1& histogram, double approximateCharge,
                             double halfWindow);

    static PeakFitResult FitRadwarePeak(const TH1& histogram, double rangeLow,
                                        double rangeHigh);

    static double EvaluateRadwarePeak(double charge, const PeakFitResult& fit);

    static PeakMatchResult MatchReferencePeaks(const TH1& target,
                                               const std::vector<double>& referenceCharges,
                                               const SearchOptions& options);

    static PeakMatchResult FindCorrespondingPeaks(
        const TH1& reference, const TH1& target,
        const std::vector<double>& referenceCharges,
        const SearchOptions& options);

    static PeakMatchResult AlignSpectrumPatterns(const TH1& reference,
                                                 const TH1& target,
                                                 const SearchOptions& options);

    static double MapReferenceCharge(const PeakMatchResult& alignment,
                                     double referenceCharge);

    static double MapTargetChargeToReference(const PeakMatchResult& alignment,
                                             double targetCharge);

    static CalibrationResult FitSecondOrder(int crystal,
                                            std::vector<CalibrationPoint> points,
                                            double residualRmsLimitKeV);
};

} // namespace hpge
