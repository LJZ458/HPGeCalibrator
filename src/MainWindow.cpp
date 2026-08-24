#include "MainWindow.h"

#include <Buttons.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGButton.h>
#include <TGClient.h>
#include <TGComboBox.h>
#include <TGFileDialog.h>
#include <TGLabel.h>
#include <TGListBox.h>
#include <TGMsgBox.h>
#include <TGNumberEntry.h>
#include <TGTab.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TList.h>
#include <TObjString.h>
#include <TPad.h>
#include <TRootEmbeddedCanvas.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace hpge {
namespace {

class PeakCanvas final : public TRootEmbeddedCanvas {
public:
    PeakCanvas(const char* name, const TGWindow* parent, UInt_t width, UInt_t height,
               MainWindow* owner)
        : TRootEmbeddedCanvas(name, parent, width, height), owner_(owner) {}

protected:
    Bool_t HandleContainerButton(Event_t* event) override {
        const Bool_t handled = TRootEmbeddedCanvas::HandleContainerButton(event);
        if (event && event->fType == kButtonPress && event->fCode == kButton1 && owner_) {
            owner_->OnCanvasEvent(kButton1Down, event->fX, event->fY, nullptr);
        }
        return handled;
    }

private:
    MainWindow* owner_ = nullptr;
};

TGLayoutHints* ExpandXY(int left = 3, int right = 3, int top = 3, int bottom = 3) {
    return new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, left, right, top, bottom);
}

TGLayoutHints* ExpandX(int left = 3, int right = 3, int top = 3, int bottom = 3) {
    return new TGLayoutHints(kLHintsExpandX, left, right, top, bottom);
}

TGLayoutHints* Left(int left = 3, int right = 3, int top = 3, int bottom = 3) {
    return new TGLayoutHints(kLHintsLeft | kLHintsCenterY, left, right, top, bottom);
}

TGTextButton* CommandButton(TGCompositeFrame* parent, const char* label, Int_t id,
                            const TGWindow* receiver) {
    auto* button = new TGTextButton(parent, label, id);
    button->Associate(receiver);
    return button;
}

std::string FormatNumber(double value, int precision = 5) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(precision) << value;
    return out.str();
}

} // namespace

MainWindow::MainWindow(const TGWindow* parent, UInt_t width, UInt_t height)
    : TGMainFrame(parent, width, height) {
    SetCleanup(kDeepCleanup);
    SetWindowName("HPGe Crystal Calibrator");
    BuildInterface();
    PopulateEnergyLines();
    MapSubwindows();
    Resize(GetDefaultSize());
    MapWindow();
    SetStatus("Add one or more ROOT files to begin.");
}

MainWindow::~MainWindow() = default;

void MainWindow::CloseWindow() { gApplication->Terminate(0); }

void MainWindow::BuildInterface() {
    auto* main = new TGHorizontalFrame(this);
    AddFrame(main, ExpandXY());

    tabs_ = new TGTab(main, 470, 800);
    auto* dataTab = tabs_->AddTab("1. Data");
    auto* peaksTab = tabs_->AddTab("2. Reference peaks");
    auto* calibrationTab = tabs_->AddTab("3. Calibration & review");
    BuildDataTab(dataTab);
    BuildPeaksTab(peaksTab);
    BuildCalibrationTab(calibrationTab);
    main->AddFrame(tabs_, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 2, 4, 2, 2));

    auto* right = new TGVerticalFrame(main);
    canvas_ = new PeakCanvas("hpgeCanvas", right, 900, 720, this);
    right->AddFrame(canvas_, ExpandXY());
    statusLabel_ = new TGLabel(right, "");
    statusLabel_->SetTextJustify(kTextLeft | kTextCenterY);
    right->AddFrame(statusLabel_, ExpandX(6, 6, 4, 6));
    main->AddFrame(right, ExpandXY());
}

void MainWindow::BuildDataTab(TGCompositeFrame* parent) {
    auto* layout = new TGVerticalFrame(parent);
    parent->AddFrame(layout, ExpandXY());

    auto* add = CommandButton(layout, "Add ROOT file...", kAddFiles, this);
    layout->AddFrame(add, ExpandX());
    layout->AddFrame(new TGLabel(layout, "Discovered TH2 histograms (select one or more):"), Left());
    histogramList_ = new TGListBox(layout, kHistogramList);
    histogramList_->Associate(this);
    histogramList_->SetMultipleSelections(kTRUE);
    layout->AddFrame(histogramList_, new TGLayoutHints(kLHintsExpandX, 3, 3, 2, 3));
    histogramList_->Resize(440, 190);

    auto* axisRow = new TGHorizontalFrame(layout);
    axisRow->AddFrame(new TGLabel(axisRow, "TH2 orientation:"), Left());
    orientationCombo_ = new TGComboBox(axisRow, kOrientation);
    orientationCombo_->Associate(this);
    orientationCombo_->AddEntry("X = charge, Y = crystal", 1);
    orientationCombo_->AddEntry("X = crystal, Y = charge", 2);
    orientationCombo_->Select(1);
    axisRow->AddFrame(orientationCombo_, ExpandX());
    layout->AddFrame(axisRow, ExpandX());

    auto* refRow = new TGHorizontalFrame(layout);
    refRow->AddFrame(new TGLabel(refRow, "Reference crystal (0-63):"), Left());
    referenceCrystalEntry_ = new TGNumberEntry(refRow, 0, 4, -1,
                                               TGNumberFormat::kNESInteger,
                                               TGNumberFormat::kNEANonNegative,
                                               TGNumberFormat::kNELLimitMinMax, 0, 63);
    refRow->AddFrame(referenceCrystalEntry_, Left());
    layout->AddFrame(refRow, ExpandX());

    layout->AddFrame(new TGLabel(layout, "Crystals to calibrate:"), Left());
    crystalList_ = new TGListBox(layout);
    crystalList_->Associate(this);
    crystalList_->SetMultipleSelections(kTRUE);
    for (int crystal = 0; crystal < 64; ++crystal) {
        crystalList_->AddEntry(("Crystal " + std::to_string(crystal)).c_str(), crystal + 1);
        crystalList_->Select(crystal + 1, kTRUE);
    }
    crystalList_->Resize(440, 260);
    layout->AddFrame(crystalList_, new TGLayoutHints(kLHintsExpandX, 3, 3, 2, 3));
    auto* crystalButtons = new TGHorizontalFrame(layout);
    crystalButtons->AddFrame(CommandButton(crystalButtons, "Select all", kSelectAllCrystals, this),
                             ExpandX());
    crystalButtons->AddFrame(CommandButton(crystalButtons, "Select none", kSelectNoCrystals, this),
                             ExpandX());
    layout->AddFrame(crystalButtons, ExpandX());
    layout->AddFrame(CommandButton(layout, "Preview reference spectrum", kPreviewReference, this),
                     ExpandX(3, 3, 8, 3));
}

void MainWindow::BuildPeaksTab(TGCompositeFrame* parent) {
    auto* layout = new TGVerticalFrame(parent);
    parent->AddFrame(layout, ExpandXY());
    layout->AddFrame(new TGLabel(layout,
        "For every source histogram, choose an energy and click its peak\non the reference-crystal spectrum."),
        Left());

    layout->AddFrame(new TGLabel(layout, "Source histogram:"), Left());
    referenceHistogramCombo_ = new TGComboBox(layout, kReferenceHistogram);
    referenceHistogramCombo_->Associate(this);
    layout->AddFrame(referenceHistogramCombo_, ExpandX());
    layout->AddFrame(CommandButton(layout, "Show reference spectrum", kPreviewReference, this), ExpandX());

    layout->AddFrame(new TGLabel(layout, "Known line for the next click:"), Left());
    energyList_ = new TGListBox(layout, kEnergyList);
    energyList_->Associate(this);
    energyList_->Resize(440, 180);
    layout->AddFrame(energyList_, new TGLayoutHints(kLHintsExpandX, 3, 3, 2, 3));

    auto* custom = new TGHorizontalFrame(layout);
    custom->AddFrame(new TGLabel(custom, "Custom energy (keV):"), Left());
    customEnergyEntry_ = new TGNumberEntry(custom, 0.0, 9, -1,
                                           TGNumberFormat::kNESRealThree,
                                           TGNumberFormat::kNEAPositive);
    custom->AddFrame(customEnergyEntry_, ExpandX());
    custom->AddFrame(CommandButton(custom, "Add", kAddCustomEnergy, this), Left());
    layout->AddFrame(custom, ExpandX());

    auto* pickRow = new TGHorizontalFrame(layout);
    pickRow->AddFrame(new TGLabel(pickRow, "Peak picking active; snap window:"), Left());
    peakWindowEntry_ = new TGNumberEntry(pickRow, 12.0, 6, -1,
                                         TGNumberFormat::kNESRealOne,
                                         TGNumberFormat::kNEAPositive);
    pickRow->AddFrame(peakWindowEntry_, Left());
    layout->AddFrame(pickRow, ExpandX());

    layout->AddFrame(new TGLabel(layout, "Selected reference peaks:"), Left());
    referencePeakList_ = new TGListBox(layout);
    referencePeakList_->Associate(this);
    referencePeakList_->Resize(440, 190);
    layout->AddFrame(referencePeakList_, new TGLayoutHints(kLHintsExpandX, 3, 3, 2, 3));
    auto* buttons = new TGHorizontalFrame(layout);
    buttons->AddFrame(CommandButton(buttons, "Remove selected", kRemoveReferencePeak, this), ExpandX());
    buttons->AddFrame(CommandButton(buttons, "Clear all", kClearReferencePeaks, this), ExpandX());
    layout->AddFrame(buttons, ExpandX());
}

void MainWindow::BuildCalibrationTab(TGCompositeFrame* parent) {
    auto* layout = new TGVerticalFrame(parent);
    parent->AddFrame(layout, ExpandXY());

    auto* settings = new TGGroupFrame(layout, "Automatic peak matching");
    auto* sigmaRow = new TGHorizontalFrame(settings);
    sigmaRow->AddFrame(new TGLabel(sigmaRow, "TSpectrum sigma (bins):"), Left());
    sigmaEntry_ = new TGNumberEntry(sigmaRow, 2.0, 6, -1, TGNumberFormat::kNESRealOne,
                                    TGNumberFormat::kNEAPositive);
    sigmaRow->AddFrame(sigmaEntry_, Left());
    settings->AddFrame(sigmaRow, ExpandX());
    auto* thresholdRow = new TGHorizontalFrame(settings);
    thresholdRow->AddFrame(new TGLabel(thresholdRow, "Peak threshold (fraction):"), Left());
    thresholdEntry_ = new TGNumberEntry(thresholdRow, 0.05, 7, -1,
                                        TGNumberFormat::kNESRealThree,
                                        TGNumberFormat::kNEAPositive);
    thresholdRow->AddFrame(thresholdEntry_, Left());
    settings->AddFrame(thresholdRow, ExpandX());
    auto* residualRow = new TGHorizontalFrame(settings);
    residualRow->AddFrame(new TGLabel(residualRow, "Review if residual RMS > keV:"), Left());
    residualLimitEntry_ = new TGNumberEntry(residualRow, 1.0, 7, -1,
                                            TGNumberFormat::kNESRealTwo,
                                            TGNumberFormat::kNEAPositive);
    residualRow->AddFrame(residualLimitEntry_, Left());
    settings->AddFrame(residualRow, ExpandX());
    layout->AddFrame(settings, ExpandX());
    layout->AddFrame(CommandButton(layout, "Calibrate selected crystals", kRunCalibration, this),
                     ExpandX(3, 3, 5, 5));

    layout->AddFrame(new TGLabel(layout,
        "Results: [OK/REVIEW/FAIL] crystal | peaks | RMS | p0, p1, p2"), Left());
    resultList_ = new TGListBox(layout, kResultList);
    resultList_->Associate(this);
    resultList_->Resize(440, 180);
    layout->AddFrame(resultList_, new TGLayoutHints(kLHintsExpandX, 3, 3, 2, 3));
    auto* viewButtons = new TGHorizontalFrame(layout);
    viewButtons->AddFrame(CommandButton(viewButtons, "Show spectrum", kShowSpectrum, this), ExpandX());
    viewButtons->AddFrame(CommandButton(viewButtons, "Fit + residuals", kShowCalibration, this), ExpandX());
    layout->AddFrame(viewButtons, ExpandX());

    auto* manual = new TGGroupFrame(layout, "Manual correction for selected result");
    manual->AddFrame(new TGLabel(manual, "Histogram and energy for next spectrum click:"), Left());
    manualHistogramCombo_ = new TGComboBox(manual, kManualHistogram);
    manualHistogramCombo_->Associate(this);
    manual->AddFrame(manualHistogramCombo_, ExpandX());
    manualEnergyList_ = new TGListBox(manual, kManualEnergyList);
    manualEnergyList_->Associate(this);
    manualEnergyList_->Resize(420, 100);
    manual->AddFrame(manualEnergyList_, new TGLayoutHints(kLHintsExpandX, 3, 3, 2, 3));
    manual->AddFrame(new TGLabel(manual, "Peak picking is active while this tab is open."), Left());
    manualPeakList_ = new TGListBox(manual);
    manualPeakList_->Associate(this);
    manualPeakList_->Resize(420, 90);
    manual->AddFrame(manualPeakList_, new TGLayoutHints(kLHintsExpandX, 3, 3, 2, 3));
    auto* manualButtons = new TGHorizontalFrame(manual);
    manualButtons->AddFrame(CommandButton(manualButtons, "Remove point", kRemoveManualPeak, this), ExpandX());
    manualButtons->AddFrame(CommandButton(manualButtons, "Refit crystal", kRefitCrystal, this), ExpandX());
    manual->AddFrame(manualButtons, ExpandX());
    layout->AddFrame(manual, ExpandX());
    layout->AddFrame(CommandButton(layout, "Export results and residuals to CSV", kExportCsv, this),
                     ExpandX(3, 3, 5, 3));
}

void MainWindow::PopulateEnergyLines() {
    if (energyLines_.empty()) {
        energyLines_ = {
            {511.000, "annihilation"}, {661.657, "Cs-137"},
            {846.771, "Co-56"}, {1037.840, "Co-56"},
            {1173.228, "Co-60"}, {1238.282, "Co-56"},
            {1274.537, "Na-22"}, {1332.492, "Co-60"},
            {1460.822, "K-40 background"}, {1771.351, "Co-56"},
            {2034.755, "Co-56"}, {2598.459, "Co-56"},
            {2614.511, "Tl-208 background"}, {3201.962, "Co-56"},
            {3253.416, "Co-56"}
        };
    }
    auto fill = [&](TGListBox* list) {
        list->RemoveAll();
        for (std::size_t i = 0; i < energyLines_.size(); ++i) {
            const auto& line = energyLines_[i];
            const std::string text = FormatNumber(line.energy, 3) + " keV — " + line.label;
            list->AddEntry(text.c_str(), static_cast<int>(i) + 1);
        }
        if (!energyLines_.empty()) list->Select(1);
        list->Layout();
    };
    fill(energyList_);
    fill(manualEnergyList_);
}

Bool_t MainWindow::ProcessMessage(Longptr_t msg, Longptr_t parm1, Longptr_t parm2) {
    if (GET_MSG(msg) != kC_COMMAND) return TGMainFrame::ProcessMessage(msg, parm1, parm2);
    if (GET_SUBMSG(msg) == kCM_BUTTON) {
        switch (parm1) {
        case kAddFiles: AddRootFiles(); break;
        case kSelectAllCrystals:
            for (int i = 1; i <= 64; ++i) crystalList_->Select(i, kTRUE);
            break;
        case kSelectNoCrystals:
            for (int i = 1; i <= 64; ++i) crystalList_->Select(i, kFALSE);
            break;
        case kPreviewReference: ShowReferenceSpectrum(); break;
        case kAddCustomEnergy: {
            const double energy = customEnergyEntry_->GetNumber();
            if (energy > 0.0) {
                energyLines_.push_back({energy, "custom"});
                PopulateEnergyLines();
                energyList_->Select(static_cast<int>(energyLines_.size()));
                manualEnergyList_->Select(static_cast<int>(energyLines_.size()));
                SetStatus("Added custom energy line " + FormatNumber(energy, 3) + " keV.");
            }
            break;
        }
        case kRemoveReferencePeak: {
            const int selected = referencePeakList_->GetSelected();
            if (selected > 0 && selected <= static_cast<int>(referencePeaks_.size())) {
                referencePeaks_.erase(referencePeaks_.begin() + selected - 1);
                RefreshReferencePeakList();
            }
            break;
        }
        case kClearReferencePeaks:
            referencePeaks_.clear(); RefreshReferencePeakList(); break;
        case kRunCalibration: RunCalibration(); break;
        case kShowSpectrum: {
            const int crystal = CurrentResultCrystal();
            const auto* descriptor = DescriptorForCombo(manualHistogramCombo_);
            if (crystal >= 0 && descriptor) ShowCrystalSpectrum(crystal, *descriptor);
            break;
        }
        case kShowCalibration: ShowSelectedCalibration(); break;
        case kRemoveManualPeak: {
            const int selected = manualPeakList_->GetSelected();
            const int crystal = CurrentResultCrystal();
            std::vector<std::size_t> indices;
            for (std::size_t i = 0; i < manualPeaks_.size(); ++i) {
                if (manualPeaks_[i].crystal == crystal) indices.push_back(i);
            }
            if (selected > 0 && selected <= static_cast<int>(indices.size())) {
                manualPeaks_.erase(manualPeaks_.begin() + indices[selected - 1]);
                RefreshManualPeakList();
            }
            break;
        }
        case kRefitCrystal: RefitSelectedCrystal(); break;
        case kExportCsv: ExportCsv(); break;
        default: break;
        }
    } else if (GET_SUBMSG(msg) == kCM_COMBOBOX) {
        if (parm1 == kReferenceHistogram) ShowReferenceSpectrum();
        if (parm1 == kManualHistogram) {
            const int crystal = CurrentResultCrystal();
            const auto* descriptor = DescriptorForCombo(manualHistogramCombo_);
            if (crystal >= 0 && descriptor) ShowCrystalSpectrum(crystal, *descriptor);
        }
    } else if (GET_SUBMSG(msg) == kCM_LISTBOX && parm1 == kResultList) {
        RefreshManualPeakList();
    }
    return kTRUE;
}

void MainWindow::AddRootFiles() {
    static const char* fileTypes[] = {"ROOT files", "*.root", "All files", "*", nullptr, nullptr};
    TGFileInfo info;
    info.fFileTypes = fileTypes;
    new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &info);
    if (!info.fFilename) return;
    std::string error;
    auto found = repository_.Discover(info.fFilename, error);
    int added = 0;
    for (auto& descriptor : found) {
        const auto duplicate = std::find_if(descriptors_.begin(), descriptors_.end(),
            [&](const HistogramDescriptor& existing) { return existing.id == descriptor.id; });
        if (duplicate == descriptors_.end()) {
            descriptors_.push_back(std::move(descriptor));
            ++added;
        }
    }
    RefreshDatasetWidgets();
    if (added > 0) {
        const int firstNew = static_cast<int>(descriptors_.size()) - added + 1;
        for (int id = firstNew; id <= static_cast<int>(descriptors_.size()); ++id) {
            histogramList_->Select(id, kTRUE);
        }
        SetStatus("Discovered " + std::to_string(added) + " TH2 histogram(s). Select active datasets.");
    } else {
        SetStatus(error.empty() ? "This file contains no new TH2 histograms." : error);
    }
}

void MainWindow::RefreshDatasetWidgets() {
    histogramList_->RemoveAll();
    referenceHistogramCombo_->RemoveEntries(0, 999999);
    manualHistogramCombo_->RemoveEntries(0, 999999);
    for (std::size_t i = 0; i < descriptors_.size(); ++i) {
        const auto& descriptor = descriptors_[i];
        const int id = static_cast<int>(i) + 1;
        const std::string item = descriptor.displayName + " [" +
                                 std::to_string(descriptor.xBins) + "x" +
                                 std::to_string(descriptor.yBins) + "]";
        histogramList_->AddEntry(item.c_str(), id);
        referenceHistogramCombo_->AddEntry(descriptor.displayName.c_str(), id);
        manualHistogramCombo_->AddEntry(descriptor.displayName.c_str(), id);
    }
    if (!descriptors_.empty()) {
        referenceHistogramCombo_->Select(1);
        manualHistogramCombo_->Select(1);
    }
    histogramList_->Layout();
    referenceHistogramCombo_->Layout();
    manualHistogramCombo_->Layout();
}

std::vector<int> MainWindow::SelectedDescriptorIndices() const {
    std::vector<int> indices;
    TList selected;
    histogramList_->GetSelectedEntries(&selected);
    TIter next(&selected);
    while (auto* entry = dynamic_cast<TGTextLBEntry*>(next())) {
        const int index = entry->EntryId() - 1;
        if (index >= 0 && index < static_cast<int>(descriptors_.size())) indices.push_back(index);
    }
    std::sort(indices.begin(), indices.end());
    return indices;
}

std::vector<int> MainWindow::SelectedCrystals() const {
    std::vector<int> crystals;
    TList selected;
    crystalList_->GetSelectedEntries(&selected);
    TIter next(&selected);
    while (auto* entry = dynamic_cast<TGTextLBEntry*>(next())) {
        if (entry->EntryId() >= 1 && entry->EntryId() <= 64) {
            crystals.push_back(entry->EntryId() - 1);
        }
    }
    std::sort(crystals.begin(), crystals.end());
    return crystals;
}

AxisOrientation MainWindow::Orientation() const {
    return orientationCombo_->GetSelected() == 2 ? AxisOrientation::ChargeOnY
                                                  : AxisOrientation::ChargeOnX;
}

int MainWindow::ReferenceCrystal() const {
    return std::clamp(static_cast<int>(referenceCrystalEntry_->GetIntNumber()), 0, 63);
}

int MainWindow::CurrentResultCrystal() const {
    const int selected = resultList_->GetSelected();
    return selected > 0 ? selected - 1 : -1;
}

const HistogramDescriptor* MainWindow::DescriptorForCombo(const TGComboBox* combo) const {
    const int index = combo->GetSelected() - 1;
    return index >= 0 && index < static_cast<int>(descriptors_.size()) ? &descriptors_[index]
                                                                       : nullptr;
}

const MainWindow::EnergyLine* MainWindow::SelectedEnergy(const TGListBox* list) const {
    const int index = list->GetSelected() - 1;
    return index >= 0 && index < static_cast<int>(energyLines_.size()) ? &energyLines_[index]
                                                                      : nullptr;
}

void MainWindow::ShowReferenceSpectrum() {
    const auto* descriptor = DescriptorForCombo(referenceHistogramCombo_);
    if (!descriptor) {
        SetStatus("Choose a source histogram first.");
        return;
    }
    ShowCrystalSpectrum(ReferenceCrystal(), *descriptor);
}

void MainWindow::ShowCrystalSpectrum(int crystal, const HistogramDescriptor& descriptor) {
    std::string error;
    auto spectrum = repository_.ProjectCrystal(descriptor, crystal, Orientation(), error);
    if (!spectrum) {
        SetStatus(error);
        return;
    }
    auto* rootCanvas = canvas_->GetCanvas();
    rootCanvas->Clear();
    displayedSpectrum_ = std::move(spectrum);
    displayedDatasetId_ = descriptor.id;
    displayedCrystal_ = crystal;
    displayedSpectrum_->SetLineColor(kBlue + 1);
    displayedSpectrum_->Draw("hist");
    rootCanvas->Update();
    SetStatus("Showing crystal " + std::to_string(crystal) + ". Click peaks only when picking is enabled.");
}

double MainWindow::ClickCharge(Int_t px) const {
    TVirtualPad* pad = canvas_->GetCanvas()->GetSelectedPad();
    if (!pad) pad = canvas_->GetCanvas();
    return pad->PadtoX(pad->AbsPixeltoX(px));
}

void MainWindow::OnCanvasEvent(Int_t event, Int_t px, Int_t, TObject*) {
    if (event != kButton1Down || !displayedSpectrum_) return;
    const double clicked = ClickCharge(px);
    if (clicked < displayedSpectrum_->GetXaxis()->GetXmin() ||
        clicked > displayedSpectrum_->GetXaxis()->GetXmax()) return;
    const double refined = CalibrationEngine::RefinePeak(
        *displayedSpectrum_, clicked, std::max(0.1, peakWindowEntry_->GetNumber()));
    if (tabs_->GetCurrent() == 1) {
        AddReferencePeak(refined);
    } else if (tabs_->GetCurrent() == 2) {
        AddManualPeak(refined);
    }
}

void MainWindow::AddReferencePeak(double charge) {
    const auto* energy = SelectedEnergy(energyList_);
    const auto* descriptor = DescriptorForCombo(referenceHistogramCombo_);
    if (!energy || !descriptor || displayedDatasetId_ != descriptor->id ||
        displayedCrystal_ != ReferenceCrystal()) {
        SetStatus("Show the selected histogram's reference spectrum and choose an energy first.");
        return;
    }
    auto duplicate = std::find_if(referencePeaks_.begin(), referencePeaks_.end(),
        [&](const ReferencePeak& peak) {
            return peak.datasetId == descriptor->id && std::abs(peak.energy - energy->energy) < 1e-6;
        });
    ReferencePeak peak{descriptor->id, charge, energy->energy, energy->label};
    if (duplicate == referencePeaks_.end()) referencePeaks_.push_back(std::move(peak));
    else *duplicate = std::move(peak);
    RefreshReferencePeakList();
    SetStatus("Reference peak: " + FormatNumber(charge, 3) + " charge -> " +
              FormatNumber(energy->energy, 3) + " keV.");
}

void MainWindow::AddManualPeak(double charge) {
    const auto* energy = SelectedEnergy(manualEnergyList_);
    const auto* descriptor = DescriptorForCombo(manualHistogramCombo_);
    const int crystal = CurrentResultCrystal();
    if (!energy || !descriptor || crystal < 0 || displayedDatasetId_ != descriptor->id ||
        displayedCrystal_ != crystal) {
        SetStatus("Select a result, histogram and energy, then show that crystal's spectrum.");
        return;
    }
    auto duplicate = std::find_if(manualPeaks_.begin(), manualPeaks_.end(),
        [&](const ManualPeak& peak) {
            return peak.datasetId == descriptor->id && peak.crystal == crystal &&
                   std::abs(peak.energy - energy->energy) < 1e-6;
        });
    ManualPeak peak{descriptor->id, crystal, charge, energy->energy, energy->label};
    if (duplicate == manualPeaks_.end()) manualPeaks_.push_back(std::move(peak));
    else *duplicate = std::move(peak);
    RefreshManualPeakList();
    SetStatus("Manual point added. Click Refit crystal to apply it.");
}

void MainWindow::RefreshReferencePeakList() {
    referencePeakList_->RemoveAll();
    for (std::size_t i = 0; i < referencePeaks_.size(); ++i) {
        const auto& peak = referencePeaks_[i];
        auto descriptor = std::find_if(descriptors_.begin(), descriptors_.end(),
            [&](const HistogramDescriptor& item) { return item.id == peak.datasetId; });
        const std::string dataset = descriptor == descriptors_.end() ? "missing" : descriptor->displayName;
        const std::string text = dataset + " | q=" + FormatNumber(peak.charge, 3) +
                                 " -> " + FormatNumber(peak.energy, 3) + " keV (" + peak.label + ")";
        referencePeakList_->AddEntry(text.c_str(), static_cast<int>(i) + 1);
    }
    referencePeakList_->Layout();
}

void MainWindow::RefreshManualPeakList() {
    manualPeakList_->RemoveAll();
    const int crystal = CurrentResultCrystal();
    int id = 1;
    for (const auto& peak : manualPeaks_) {
        if (peak.crystal != crystal) continue;
        const std::string text = "q=" + FormatNumber(peak.charge, 3) + " -> " +
                                 FormatNumber(peak.energy, 3) + " keV (" + peak.label + ")";
        manualPeakList_->AddEntry(text.c_str(), id++);
    }
    manualPeakList_->Layout();
}

std::vector<CalibrationPoint> MainWindow::BuildPointsForCrystal(int crystal) {
    std::vector<CalibrationPoint> points;
    CalibrationEngine::SearchOptions options;
    options.sigmaBins = sigmaEntry_->GetNumber();
    options.threshold = thresholdEntry_->GetNumber();

    for (int descriptorIndex : SelectedDescriptorIndices()) {
        const auto& descriptor = descriptors_[descriptorIndex];
        std::vector<ReferencePeak> references;
        for (const auto& peak : referencePeaks_) {
            if (peak.datasetId == descriptor.id) references.push_back(peak);
        }
        if (references.empty()) continue;
        std::sort(references.begin(), references.end(),
                  [](const ReferencePeak& a, const ReferencePeak& b) { return a.charge < b.charge; });
        if (crystal == ReferenceCrystal()) {
            for (const auto& ref : references) {
                points.push_back({descriptor.id, ref.charge, ref.energy, 0.0, false, 0.0});
            }
            continue;
        }

        std::string error;
        auto spectrum = repository_.ProjectCrystal(descriptor, crystal, Orientation(), error);
        if (!spectrum) continue;
        std::vector<double> referenceCharges;
        for (const auto& ref : references) referenceCharges.push_back(ref.charge);
        const auto matches = CalibrationEngine::MatchReferencePeaks(*spectrum, referenceCharges, options);
        for (std::size_t i = 0; i < references.size(); ++i) {
            if (i < matches.matched.size() && matches.matched[i]) {
                const double refined = CalibrationEngine::RefinePeak(
                    *spectrum, matches.charges[i], std::max(0.1, peakWindowEntry_->GetNumber()));
                points.push_back({descriptor.id, refined, references[i].energy, 0.0, false, 0.0});
            }
        }
    }

    // A manual point replaces the automatic point with the same dataset and energy.
    for (const auto& manual : manualPeaks_) {
        if (manual.crystal != crystal) continue;
        auto existing = std::find_if(points.begin(), points.end(), [&](const CalibrationPoint& point) {
            return point.datasetId == manual.datasetId &&
                   std::abs(point.energy - manual.energy) < 1e-6;
        });
        CalibrationPoint replacement{manual.datasetId, manual.charge, manual.energy, 0.0, true, 0.0};
        if (existing == points.end()) points.push_back(replacement);
        else *existing = replacement;
    }
    return points;
}

CalibrationResult MainWindow::CalibrateCrystal(int crystal) {
    auto points = BuildPointsForCrystal(crystal);
    return CalibrationEngine::FitSecondOrder(crystal, std::move(points),
                                              residualLimitEntry_->GetNumber());
}

void MainWindow::RunCalibration() {
    const auto datasets = SelectedDescriptorIndices();
    const auto crystals = SelectedCrystals();
    if (datasets.empty() || crystals.empty()) {
        SetStatus("Select at least one TH2 histogram and one crystal.");
        return;
    }
    int referenceCount = 0;
    for (int index : datasets) {
        referenceCount += static_cast<int>(std::count_if(referencePeaks_.begin(), referencePeaks_.end(),
            [&](const ReferencePeak& peak) { return peak.datasetId == descriptors_[index].id; }));
    }
    if (referenceCount < 3) {
        SetStatus("At least three reference peaks across the selected histograms are required.");
        return;
    }
    SetStatus("Calibrating " + std::to_string(crystals.size()) + " crystals...");
    gSystem->ProcessEvents();
    results_.clear();
    for (int crystal : crystals) results_[crystal] = CalibrateCrystal(crystal);
    RefreshResults();
    int ok = 0, review = 0, failed = 0;
    for (const auto& [crystal, result] : results_) {
        (void)crystal;
        if (!result.success) ++failed;
        else if (result.needsReview) ++review;
        else ++ok;
    }
    SetStatus("Calibration complete: " + std::to_string(ok) + " OK, " +
              std::to_string(review) + " review, " + std::to_string(failed) + " failed.");
}

void MainWindow::RefitSelectedCrystal() {
    const int crystal = CurrentResultCrystal();
    if (crystal < 0) {
        SetStatus("Select a calibration result to refit.");
        return;
    }
    results_[crystal] = CalibrateCrystal(crystal);
    RefreshResults();
    resultList_->Select(crystal + 1);
    ShowSelectedCalibration();
    SetStatus("Crystal " + std::to_string(crystal) + " refitted with manual overrides.");
}

void MainWindow::RefreshResults() {
    const int previous = resultList_->GetSelected();
    resultList_->RemoveAll();
    for (const auto& [crystal, result] : results_) {
        const std::string state = !result.success ? "FAIL" : result.needsReview ? "REVIEW" : "OK";
        std::ostringstream row;
        row << '[' << state << "] C" << std::setw(2) << std::setfill('0') << crystal
            << " | n=" << result.points.size() << " | RMS=" << std::fixed << std::setprecision(3)
            << result.residualRms << " | " << std::scientific << std::setprecision(5)
            << result.p0 << ", " << result.p1 << ", " << result.p2;
        resultList_->AddEntry(row.str().c_str(), crystal + 1);
    }
    resultList_->Layout();
    if (previous > 0 && results_.count(previous - 1)) resultList_->Select(previous);
    else if (!results_.empty()) resultList_->Select(results_.begin()->first + 1);
    RefreshManualPeakList();
}

void MainWindow::ShowSelectedCalibration() {
    const int crystal = CurrentResultCrystal();
    const auto found = results_.find(crystal);
    if (found == results_.end() || !found->second.success) {
        SetStatus(found == results_.end() ? "Select a result first." : found->second.status);
        return;
    }
    displayedSpectrum_.reset();
    const auto& result = found->second;
    std::vector<double> x, y, residuals;
    for (const auto& point : result.points) {
        x.push_back(point.charge);
        y.push_back(point.energy);
        residuals.push_back(point.residual);
    }
    auto* rootCanvas = canvas_->GetCanvas();
    rootCanvas->Clear();
    rootCanvas->Divide(1, 2);
    rootCanvas->cd(1);
    TGraph graph(static_cast<int>(x.size()), x.data(), y.data());
    graph.SetTitle(("Crystal " + std::to_string(crystal) +
                    " calibration;Charge;Energy (keV)").c_str());
    graph.SetMarkerStyle(20);
    graph.SetMarkerColor(kBlue + 1);
    graph.DrawClone("AP");
    const double xMinimum = *std::min_element(x.begin(), x.end());
    const double xMaximum = *std::max_element(x.begin(), x.end());
    std::vector<double> fitX(200), fitY(200);
    for (std::size_t i = 0; i < fitX.size(); ++i) {
        fitX[i] = xMinimum + (xMaximum - xMinimum) * static_cast<double>(i) /
                               static_cast<double>(fitX.size() - 1);
        fitY[i] = result.p0 + result.p1 * fitX[i] + result.p2 * fitX[i] * fitX[i];
    }
    TGraph fitGraph(static_cast<int>(fitX.size()), fitX.data(), fitY.data());
    fitGraph.SetLineColor(kRed + 1);
    fitGraph.SetLineWidth(2);
    fitGraph.DrawClone("L same");
    rootCanvas->cd(2);
    TGraph residualGraph(static_cast<int>(x.size()), x.data(), residuals.data());
    residualGraph.SetTitle(("Residuals (RMS " + FormatNumber(result.residualRms, 3) +
                            " keV);Charge;Energy - fit (keV)").c_str());
    residualGraph.SetMarkerStyle(21);
    residualGraph.SetMarkerColor(result.needsReview ? kOrange + 7 : kGreen + 2);
    residualGraph.DrawClone("AP");
    TLine zero(*std::min_element(x.begin(), x.end()), 0.0,
               *std::max_element(x.begin(), x.end()), 0.0);
    zero.SetLineStyle(2);
    zero.DrawClone();
    rootCanvas->cd(0);
    rootCanvas->Update();
    SetStatus("Crystal " + std::to_string(crystal) + ": " + result.status);
}

void MainWindow::ExportCsv() {
    if (results_.empty()) {
        SetStatus("There are no calibration results to export.");
        return;
    }
    static const char* fileTypes[] = {"CSV files", "*.csv", "All files", "*", nullptr, nullptr};
    TGFileInfo info;
    info.fFileTypes = fileTypes;
    info.fFilename = StrDup("hpge_calibration.csv");
    new TGFileDialog(gClient->GetRoot(), this, kFDSave, &info);
    if (!info.fFilename) return;
    std::ofstream output(info.fFilename);
    if (!output) {
        SetStatus("Could not write CSV file.");
        return;
    }
    output << "record,crystal,status,needs_review,p0,p1,p2,chi2,ndf,residual_rms_keV,"
              "dataset,charge,energy_keV,residual_keV,manual\n";
    output << std::setprecision(12);
    for (const auto& [crystal, result] : results_) {
        output << "fit," << crystal << ",\"" << result.status << "\"," << result.needsReview
               << ',' << result.p0 << ',' << result.p1 << ',' << result.p2 << ','
               << result.chi2 << ',' << result.ndf << ',' << result.residualRms
               << ",,,,,\n";
        for (const auto& point : result.points) {
            output << "point," << crystal << ",,," << result.p0 << ',' << result.p1 << ','
                   << result.p2 << ",,,," << '"' << point.datasetId << '"' << ','
                   << point.charge << ',' << point.energy << ',' << point.residual << ','
                   << point.manual << '\n';
        }
    }
    SetStatus("Exported coefficients and per-peak residuals to " + std::string(info.fFilename));
}

void MainWindow::SetStatus(const std::string& text) {
    statusLabel_->SetText(text.c_str());
    statusLabel_->Layout();
}

} // namespace hpge
