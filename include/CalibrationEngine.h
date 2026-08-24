#pragma once

#include "CalibrationTypes.h"

#include <vector>

class TH1;

namespace hpge {

class CalibrationEngine {
public:
    struct SearchOptions {
        double sigmaBins = 2.0;
        double threshold = 0.05;
        int maxPeaks = 80;
        double matchToleranceFraction = 0.018;
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

    static CalibrationResult FitSecondOrder(int crystal,
                                            std::vector<CalibrationPoint> points,
                                            double residualRmsLimitKeV);
};

} // namespace hpge
