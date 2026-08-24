#pragma once

#include "CalibrationEngine.h"
#include "CalibrationTypes.h"
#include "CombinedSpectrumAnalyzer.h"
#include "RootDataRepository.h"

#include <QMainWindow>

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

class QComboBox;
class QCheckBox;
class QDoubleSpinBox;
class QGroupBox;
class QLabel;
class QListWidget;
class QSpinBox;
class QTabWidget;
class QWidget;
class TH1D;

namespace hpge {

class SpectrumPlotWidget;

class MainWindow final : public QMainWindow {
public:
    explicit MainWindow(QWidget* parent = nullptr);
    ~MainWindow() override = default;
    bool OpenRootFiles(const std::vector<std::string>& files);

private:
    struct EnergyLine {
        double energy = 0.0;
        std::string label;
        std::string source;
    };
    struct ManualPeak {
        std::string datasetId;
        int crystal = -1;
        double charge = 0.0;
        double energy = 0.0;
        std::string label;
        PeakFitResult peakFit;
    };
    struct CombinedDatasetAnalysis {
        std::string datasetId;
        std::shared_ptr<TH1D> spectrum;
        std::vector<CombinedPeakQuality> peaks;
        int crystalCount = 0;
    };

    void BuildInterface();
    QWidget* BuildDataTab();
    QWidget* BuildPeaksTab();
    QWidget* BuildCalibrationTab();
    void ConnectActions();
    void PopulateEnergyLines();
    void RefreshEnergyList(QListWidget* list, const QComboBox* sourceCombo,
                           std::vector<std::size_t>& indices);
    void AddRootFiles();
    void RefreshDatasetWidgets();
    std::vector<int> SelectedDescriptorIndices() const;
    std::vector<int> SelectedCrystals() const;
    AxisOrientation Orientation() const;
    double AlignmentSensitivity() const;
    CalibrationEngine::AlignmentModel AlignmentModel() const;
    int ReferenceCrystal() const;
    int CurrentResultCrystal() const;
    const HistogramDescriptor* DescriptorForCombo(const QComboBox* combo) const;
    const EnergyLine* SelectedEnergy(const QListWidget* list) const;
    void ShowReferenceSpectrum();
    void ShowCrystalSpectrum(int crystal, const HistogramDescriptor& descriptor);
    void RedrawDisplayedSpectrum(bool preserveView = true);
    void ShowSelectedResultSpectrum();
    void UpdateManualCorrectionForSelection();
    void HandleRangeClick(double charge);
    void AddReferencePeak(const PeakFitResult& fit);
    void AddManualPeak(const PeakFitResult& fit);
    void RefreshReferencePeakList();
    void RefreshManualPeakList();
    std::vector<CalibrationPoint> BuildPointsForCrystal(int crystal);
    CalibrationResult CalibrateCrystal(int crystal);
    void RunCalibration();
    void EvaluateCombinedSpectra();
    void RefreshCombinedQualityList();
    void ShowCombinedSpectrum();
    void ShowSpectrumAlignment();
    void RefitSelectedCrystal();
    void RefreshResults();
    void ShowSelectedCalibration();
    void ExportCsv();
    void SetStatus(const std::string& text);
    void UpdateInteractionMode();
    void SetSecondaryPlotVisible(bool visible);

    RootDataRepository repository_;
    std::vector<HistogramDescriptor> descriptors_;
    std::vector<EnergyLine> energyLines_;
    std::vector<std::string> energySources_;
    std::vector<std::size_t> referenceEnergyIndices_;
    std::vector<std::size_t> manualEnergyIndices_;
    std::vector<ReferencePeak> referencePeaks_;
    std::vector<ManualPeak> manualPeaks_;
    std::map<int, CalibrationResult> results_;
    std::map<std::string, CombinedDatasetAnalysis> combinedAnalyses_;
    std::shared_ptr<TH1D> displayedSpectrum_;
    std::string displayedDatasetId_;
    int displayedCrystal_ = -1;
    std::optional<double> pendingRangeStart_;
    bool updatingWidgets_ = false;

    QTabWidget* tabs_ = nullptr;
    QListWidget* histogramList_ = nullptr;
    QListWidget* crystalList_ = nullptr;
    QComboBox* orientationCombo_ = nullptr;
    QSpinBox* referenceCrystalEntry_ = nullptr;
    QComboBox* referenceHistogramCombo_ = nullptr;
    QComboBox* referenceSourceCombo_ = nullptr;
    QListWidget* energyList_ = nullptr;
    QDoubleSpinBox* customEnergyEntry_ = nullptr;
    QListWidget* referencePeakList_ = nullptr;
    QDoubleSpinBox* sigmaEntry_ = nullptr;
    QDoubleSpinBox* thresholdEntry_ = nullptr;
    QDoubleSpinBox* residualLimitEntry_ = nullptr;
    QDoubleSpinBox* alignmentSensitivityEntry_ = nullptr;
    QCheckBox* autoTuneAlignmentEntry_ = nullptr;
    QComboBox* alignmentModelCombo_ = nullptr;
    QComboBox* alignmentHistogramCombo_ = nullptr;
    QSpinBox* alignmentCrystalEntry_ = nullptr;
    QListWidget* resultList_ = nullptr;
    QComboBox* combinedHistogramCombo_ = nullptr;
    QListWidget* combinedQualityList_ = nullptr;
    QComboBox* manualHistogramCombo_ = nullptr;
    QComboBox* manualSourceCombo_ = nullptr;
    QListWidget* manualEnergyList_ = nullptr;
    QListWidget* manualPeakList_ = nullptr;
    QGroupBox* manualCorrectionGroup_ = nullptr;
    QComboBox* mouseModeCombo_ = nullptr;
    SpectrumPlotWidget* primaryPlot_ = nullptr;
    SpectrumPlotWidget* secondaryPlot_ = nullptr;
    QWidget* secondaryPlotContainer_ = nullptr;
    QLabel* statusLabel_ = nullptr;
};

} // namespace hpge
