#pragma once

#include "CalibrationEngine.h"
#include "CalibrationTypes.h"
#include "RootDataRepository.h"

#include <TGFrame.h>

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

class TGComboBox;
class TGLabel;
class TGListBox;
class TGNumberEntry;
class TGTab;
class TRootEmbeddedCanvas;
class TH1D;
class TTimer;

namespace hpge {

class MainWindow : public TGMainFrame {
public:
    MainWindow(const TGWindow* parent, UInt_t width, UInt_t height);
    ~MainWindow() override;

    Bool_t ProcessMessage(Longptr_t msg, Longptr_t parm1, Longptr_t parm2) override;
    void CloseWindow() override;
    void OnCanvasEvent(Int_t event, Int_t px, Int_t py, TObject* selected);
    void ProcessPendingCanvasClick();
    bool CanvasPeakPickingEnabled() const;

private:
    enum WidgetId {
        kAddFiles = 100,
        kSelectAllCrystals,
        kSelectNoCrystals,
        kPreviewReference,
        kHistogramList,
        kReferenceHistogram,
        kOrientation,
        kEnergyList,
        kRemoveReferencePeak,
        kClearReferencePeaks,
        kAddCustomEnergy,
        kRunCalibration,
        kShowAlignment,
        kAlignmentHistogram,
        kResultList,
        kShowSpectrum,
        kShowCalibration,
        kManualHistogram,
        kManualEnergyList,
        kRemoveManualPeak,
        kRefitCrystal,
        kExportCsv
    };

    struct EnergyLine {
        double energy = 0.0;
        std::string label;
    };

    struct ManualPeak {
        std::string datasetId;
        int crystal = -1;
        double charge = 0.0;
        double energy = 0.0;
        std::string label;
        PeakFitResult peakFit;
    };

    void BuildInterface();
    void BuildDataTab(TGCompositeFrame* parent);
    void BuildPeaksTab(TGCompositeFrame* parent);
    void BuildCalibrationTab(TGCompositeFrame* parent);
    void PopulateEnergyLines();
    void AddRootFiles();
    void RefreshDatasetWidgets();
    std::vector<int> SelectedDescriptorIndices() const;
    std::vector<int> SelectedCrystals() const;
    AxisOrientation Orientation() const;
    int ReferenceCrystal() const;
    int CurrentResultCrystal() const;
    const HistogramDescriptor* DescriptorForCombo(const TGComboBox* combo) const;
    const EnergyLine* SelectedEnergy(const TGListBox* list) const;
    void ShowReferenceSpectrum();
    void ShowCrystalSpectrum(int crystal, const HistogramDescriptor& descriptor);
    void RedrawDisplayedSpectrum();
    void DrawSpectrumOverlays();
    void DrawPeakFitOverlay(const PeakFitResult& fit, int color, int lineStyle,
                            const std::string& label);
    void HandleRangeClick(double charge);
    void AddReferencePeak(const PeakFitResult& fit);
    void AddManualPeak(const PeakFitResult& fit);
    void RefreshReferencePeakList();
    void RefreshManualPeakList();
    std::vector<CalibrationPoint> BuildPointsForCrystal(int crystal);
    CalibrationResult CalibrateCrystal(int crystal);
    void RunCalibration();
    void ShowSpectrumAlignment();
    void RefitSelectedCrystal();
    void RefreshResults();
    void ShowSelectedCalibration();
    void ExportCsv();
    void SetStatus(const std::string& text);
    double ClickCharge(Int_t px) const;

    RootDataRepository repository_;
    std::vector<HistogramDescriptor> descriptors_;
    std::vector<EnergyLine> energyLines_;
    std::vector<ReferencePeak> referencePeaks_;
    std::vector<ManualPeak> manualPeaks_;
    std::map<int, CalibrationResult> results_;

    std::shared_ptr<TH1D> displayedSpectrum_;
    std::string displayedDatasetId_;
    int displayedCrystal_ = -1;
    std::optional<double> pendingRangeStart_;
    int pendingClickPixelX_ = 0;
    bool pendingCanvasClick_ = false;
    bool updatingWidgets_ = false;
    TTimer* peakClickTimer_ = nullptr;

    TGTab* tabs_ = nullptr;
    TGListBox* histogramList_ = nullptr;
    TGListBox* crystalList_ = nullptr;
    TGComboBox* orientationCombo_ = nullptr;
    TGNumberEntry* referenceCrystalEntry_ = nullptr;
    TGComboBox* referenceHistogramCombo_ = nullptr;
    TGListBox* energyList_ = nullptr;
    TGNumberEntry* customEnergyEntry_ = nullptr;
    TGListBox* referencePeakList_ = nullptr;
    TGNumberEntry* sigmaEntry_ = nullptr;
    TGNumberEntry* thresholdEntry_ = nullptr;
    TGNumberEntry* residualLimitEntry_ = nullptr;
    TGComboBox* alignmentHistogramCombo_ = nullptr;
    TGNumberEntry* alignmentCrystalEntry_ = nullptr;
    TGListBox* resultList_ = nullptr;
    TGComboBox* manualHistogramCombo_ = nullptr;
    TGListBox* manualEnergyList_ = nullptr;
    TGListBox* manualPeakList_ = nullptr;
    TRootEmbeddedCanvas* canvas_ = nullptr;
    TGLabel* statusLabel_ = nullptr;

};

} // namespace hpge
