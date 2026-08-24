#include "MainWindow.h"

#include "SpectrumPlotWidget.h"

#include <QApplication>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QFrame>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QListWidget>
#include <QMenuBar>
#include <QPushButton>
#include <QScrollArea>
#include <QSpinBox>
#include <QSplitter>
#include <QStatusBar>
#include <QTabWidget>
#include <QVBoxLayout>

#include <TH1D.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace hpge {
namespace {

QString Text(const std::string& value) { return QString::fromStdString(value); }

std::string FormatNumber(double value, int precision = 5) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(precision) << value;
    return out.str();
}

QLabel* Hint(const QString& text) {
    auto* label = new QLabel(text);
    label->setWordWrap(true);
    label->setProperty("class", "hint");
    return label;
}

QWidget* Row(std::initializer_list<QWidget*> widgets) {
    auto* widget = new QWidget;
    auto* layout = new QHBoxLayout(widget);
    layout->setContentsMargins(0, 0, 0, 0);
    for (auto* child : widgets) layout->addWidget(child);
    return widget;
}

QScrollArea* ScrollPage(QWidget* content) {
    auto* scroll = new QScrollArea;
    scroll->setWidgetResizable(true);
    scroll->setFrameShape(QFrame::NoFrame);
    scroll->setWidget(content);
    return scroll;
}

QDoubleSpinBox* RealEntry(double value, int decimals, double minimum, double maximum) {
    auto* entry = new QDoubleSpinBox;
    entry->setDecimals(decimals);
    entry->setRange(minimum, maximum);
    entry->setValue(value);
    entry->setKeyboardTracking(false);
    return entry;
}

PlotSeries HistogramSeries(const TH1D& histogram, const QColor& color,
                           const std::string& name = {}) {
    PlotSeries series;
    series.name = name;
    series.color = color;
    series.width = 1;
    series.x.reserve(static_cast<std::size_t>(histogram.GetNbinsX()));
    series.y.reserve(static_cast<std::size_t>(histogram.GetNbinsX()));
    for (int bin = 1; bin <= histogram.GetNbinsX(); ++bin) {
        series.x.push_back(histogram.GetXaxis()->GetBinCenter(bin));
        series.y.push_back(histogram.GetBinContent(bin));
    }
    return series;
}

PlotSeries FitSeries(const PeakFitResult& fit, const std::string& name, bool dashed = false) {
    PlotSeries series;
    series.name = name;
    series.color = QColor("#dc2626");
    series.width = 2;
    series.dashed = dashed;
    constexpr std::size_t samples = 180;
    series.x.resize(samples);
    series.y.resize(samples);
    for (std::size_t i = 0; i < samples; ++i) {
        const double fraction = static_cast<double>(i) / static_cast<double>(samples - 1);
        series.x[i] = fit.rangeLow + fraction * (fit.rangeHigh - fit.rangeLow);
        series.y[i] = CalibrationEngine::EvaluateRadwarePeak(series.x[i], fit);
    }
    return series;
}

} // namespace

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent) {
    setWindowTitle("HPGe Crystal Calibrator");
    resize(1480, 920);
    setMinimumSize(1120, 720);
    BuildInterface();
    PopulateEnergyLines();
    ConnectActions();
    SetStatus("Add one or more ROOT files to begin.");
}

void MainWindow::BuildInterface() {
    auto* fileMenu = menuBar()->addMenu("&File");
    auto* openAction = fileMenu->addAction("&Add ROOT files...");
    openAction->setShortcut(QKeySequence::Open);
    auto* exportAction = fileMenu->addAction("Export calibration CSV...");
    fileMenu->addSeparator();
    auto* quitAction = fileMenu->addAction("Quit");
    connect(openAction, &QAction::triggered, this, [this] { AddRootFiles(); });
    connect(exportAction, &QAction::triggered, this, [this] { ExportCsv(); });
    connect(quitAction, &QAction::triggered, this, &QWidget::close);

    auto* central = new QWidget;
    auto* outer = new QHBoxLayout(central);
    outer->setContentsMargins(12, 12, 12, 12);
    outer->setSpacing(12);

    tabs_ = new QTabWidget;
    tabs_->setMinimumWidth(430);
    tabs_->setMaximumWidth(525);
    tabs_->addTab(BuildDataTab(), "1  Data");
    tabs_->addTab(BuildPeaksTab(), "2  Reference peaks");
    tabs_->addTab(BuildCalibrationTab(), "3  Calibrate + review");
    outer->addWidget(tabs_);

    auto* plotCard = new QFrame;
    plotCard->setObjectName("plotCard");
    auto* plotLayout = new QVBoxLayout(plotCard);
    plotLayout->setContentsMargins(12, 12, 12, 12);
    auto* modeLabel = new QLabel("Mouse interaction");
    mouseModeCombo_ = new QComboBox;
    mouseModeCombo_->addItems({"Select peak-fit range", "Zoom / pan"});
    mouseModeCombo_->setToolTip("Selection mode records two fit boundaries. Zoom mode uses the wheel, left-drag zoom, right-drag pan, and double-click reset.");
    auto* previousView = new QPushButton("Back");
    auto* zoomOut = new QPushButton("Zoom −");
    auto* zoomIn = new QPushButton("Zoom +");
    auto* resetView = new QPushButton("Reset");
    previousView->setToolTip("Return to the previous axis range");
    zoomOut->setToolTip("Zoom out around the center of the current view");
    zoomIn->setToolTip("Zoom in around the center of the current view");
    resetView->setToolTip("Show the complete spectrum range");
    plotLayout->addWidget(Row({modeLabel, mouseModeCombo_, previousView, zoomOut, zoomIn, resetView}));

    auto* plotSplitter = new QSplitter(Qt::Vertical);
    primaryPlot_ = new SpectrumPlotWidget;
    secondaryPlotContainer_ = new QWidget;
    auto* secondaryLayout = new QVBoxLayout(secondaryPlotContainer_);
    secondaryLayout->setContentsMargins(0, 0, 0, 0);
    secondaryPlot_ = new SpectrumPlotWidget;
    secondaryPlot_->setMinimumHeight(220);
    secondaryLayout->addWidget(secondaryPlot_);
    plotSplitter->addWidget(primaryPlot_);
    plotSplitter->addWidget(secondaryPlotContainer_);
    plotSplitter->setStretchFactor(0, 3);
    plotSplitter->setStretchFactor(1, 2);
    plotLayout->addWidget(plotSplitter, 1);
    outer->addWidget(plotCard, 1);
    setCentralWidget(central);

    statusLabel_ = new QLabel;
    statusLabel_->setWordWrap(true);
    statusBar()->addWidget(statusLabel_, 1);
    SetSecondaryPlotVisible(false);

    connect(resetView, &QPushButton::clicked, primaryPlot_, &SpectrumPlotWidget::ResetView);
    connect(resetView, &QPushButton::clicked, secondaryPlot_, &SpectrumPlotWidget::ResetView);
    connect(previousView, &QPushButton::clicked, primaryPlot_, &SpectrumPlotWidget::PreviousView);
    connect(zoomOut, &QPushButton::clicked, primaryPlot_, &SpectrumPlotWidget::ZoomOut);
    connect(zoomIn, &QPushButton::clicked, primaryPlot_, &SpectrumPlotWidget::ZoomIn);

    setStyleSheet(R"(
        QMainWindow { background: #eef1f5; }
        QFrame#plotCard { background: palette(base); border: 1px solid #c9ced6; border-radius: 6px; }
        QTabWidget::pane { border: 1px solid #c9ced6; background: palette(base); }
        QTabBar::tab { padding: 9px 12px; }
        QGroupBox { font-weight: 600; margin-top: 12px; padding-top: 12px; }
        QGroupBox::title { subcontrol-origin: margin; left: 8px; padding: 0 4px; }
        QPushButton { min-height: 26px; padding: 3px 10px; }
        QListWidget { border: 1px solid #b8bec7; border-radius: 3px; }
        QLabel[class="hint"] { color: #59636e; }
        QStatusBar { background: #e2e7ed; border-top: 1px solid #c9ced6; }
    )");
}

QWidget* MainWindow::BuildDataTab() {
    auto* page = new QWidget;
    auto* layout = new QVBoxLayout(page);
    layout->setContentsMargins(12, 12, 12, 12);
    auto* add = new QPushButton("Add ROOT files...");
    add->setIcon(style()->standardIcon(QStyle::SP_DialogOpenButton));
    connect(add, &QPushButton::clicked, this, [this] { AddRootFiles(); });
    layout->addWidget(add);
    layout->addWidget(new QLabel("Discovered TH2 histograms"));
    layout->addWidget(Hint("Select one or more charge-versus-crystal datasets."));
    histogramList_ = new QListWidget;
    histogramList_->setSelectionMode(QAbstractItemView::ExtendedSelection);
    histogramList_->setMinimumHeight(150);
    layout->addWidget(histogramList_);

    orientationCombo_ = new QComboBox;
    orientationCombo_->addItems({"X = charge, Y = crystal", "X = crystal, Y = charge"});
    layout->addWidget(new QLabel("TH2 orientation"));
    layout->addWidget(orientationCombo_);
    referenceCrystalEntry_ = new QSpinBox;
    referenceCrystalEntry_->setRange(0, 63);
    layout->addWidget(Row({new QLabel("Reference crystal (0-63)"), referenceCrystalEntry_}));

    layout->addWidget(new QLabel("Crystals to calibrate"));
    crystalList_ = new QListWidget;
    crystalList_->setSelectionMode(QAbstractItemView::ExtendedSelection);
    crystalList_->setMinimumHeight(250);
    for (int crystal = 0; crystal < 64; ++crystal) {
        auto* item = new QListWidgetItem(QString("Crystal %1").arg(crystal, 2, 10, QChar('0')),
                                         crystalList_);
        item->setData(Qt::UserRole, crystal);
        item->setSelected(true);
    }
    layout->addWidget(crystalList_);
    auto* all = new QPushButton("Select all");
    auto* none = new QPushButton("Select none");
    connect(all, &QPushButton::clicked, crystalList_, &QListWidget::selectAll);
    connect(none, &QPushButton::clicked, crystalList_, &QListWidget::clearSelection);
    layout->addWidget(Row({all, none}));
    auto* preview = new QPushButton("Preview reference spectrum");
    connect(preview, &QPushButton::clicked, this, [this] { ShowReferenceSpectrum(); });
    layout->addWidget(preview);
    layout->addStretch();
    return ScrollPage(page);
}

QWidget* MainWindow::BuildPeaksTab() {
    auto* page = new QWidget;
    auto* layout = new QVBoxLayout(page);
    layout->setContentsMargins(12, 12, 12, 12);
    layout->addWidget(Hint("Choose an energy, then click the lower and upper limits around its peak. The RadWare centroid—not the clicked limits—is used."));
    layout->addWidget(new QLabel("Source histogram"));
    referenceHistogramCombo_ = new QComboBox;
    layout->addWidget(referenceHistogramCombo_);
    auto* show = new QPushButton("Show reference spectrum");
    connect(show, &QPushButton::clicked, this, [this] { ShowReferenceSpectrum(); });
    layout->addWidget(show);
    layout->addWidget(new QLabel("Radioactive source"));
    referenceSourceCombo_ = new QComboBox;
    layout->addWidget(referenceSourceCombo_);
    layout->addWidget(new QLabel("Known line for the next fit"));
    energyList_ = new QListWidget;
    energyList_->setMinimumHeight(150);
    layout->addWidget(energyList_);

    customEnergyEntry_ = RealEntry(0.0, 3, 0.0, 100000.0);
    customEnergyEntry_->setSuffix(" keV");
    auto* addCustom = new QPushButton("Add custom line");
    connect(addCustom, &QPushButton::clicked, this, [this] {
        const double energy = customEnergyEntry_->value();
        if (energy <= 0.0) return;
        energyLines_.push_back({energy, "user supplied", "Custom"});
        PopulateEnergyLines();
        referenceSourceCombo_->setCurrentText("Custom");
        manualSourceCombo_->setCurrentText("Custom");
        RefreshEnergyList(energyList_, referenceSourceCombo_, referenceEnergyIndices_);
        RefreshEnergyList(manualEnergyList_, manualSourceCombo_, manualEnergyIndices_);
        energyList_->setCurrentRow(energyList_->count() - 1);
        manualEnergyList_->setCurrentRow(manualEnergyList_->count() - 1);
        SetStatus("Added custom energy line " + FormatNumber(energy, 3) + " keV.");
    });
    layout->addWidget(Row({customEnergyEntry_, addCustom}));
    layout->addWidget(new QLabel("Selected reference peaks"));
    referencePeakList_ = new QListWidget;
    referencePeakList_->setMinimumHeight(190);
    layout->addWidget(referencePeakList_);
    auto* remove = new QPushButton("Remove selected");
    auto* clear = new QPushButton("Clear all");
    connect(remove, &QPushButton::clicked, this, [this] {
        const int row = referencePeakList_->currentRow();
        if (row < 0 || row >= static_cast<int>(referencePeaks_.size())) return;
        referencePeaks_.erase(referencePeaks_.begin() + row);
        RefreshReferencePeakList();
        RedrawDisplayedSpectrum();
    });
    connect(clear, &QPushButton::clicked, this, [this] {
        referencePeaks_.clear();
        RefreshReferencePeakList();
        RedrawDisplayedSpectrum();
    });
    layout->addWidget(Row({remove, clear}));
    layout->addStretch();
    return ScrollPage(page);
}

QWidget* MainWindow::BuildCalibrationTab() {
    auto* page = new QWidget;
    auto* layout = new QVBoxLayout(page);
    layout->setContentsMargins(12, 12, 12, 12);
    auto* settings = new QGroupBox("Automatic peak matching");
    auto* settingsLayout = new QVBoxLayout(settings);
    sigmaEntry_ = RealEntry(2.0, 1, 0.1, 50.0);
    thresholdEntry_ = RealEntry(0.05, 3, 0.001, 1.0);
    residualLimitEntry_ = RealEntry(1.0, 2, 0.0, 1000.0);
    settingsLayout->addWidget(Row({new QLabel("TSpectrum sigma (bins)"), sigmaEntry_}));
    settingsLayout->addWidget(Row({new QLabel("Peak threshold (fraction)"), thresholdEntry_}));
    settingsLayout->addWidget(Row({new QLabel("Review if residual RMS > keV"), residualLimitEntry_}));
    layout->addWidget(settings);
    auto* run = new QPushButton("Calibrate selected crystals");
    run->setProperty("class", "primary");
    connect(run, &QPushButton::clicked, this, [this] { RunCalibration(); });
    layout->addWidget(run);

    auto* alignment = new QGroupBox("Pre-calibration spectrum alignment");
    auto* alignmentLayout = new QVBoxLayout(alignment);
    alignmentLayout->addWidget(Hint("Rough diagnostic overlay after intensity-independent peak-pattern mapping; calibration still refits every mapped peak."));
    alignmentHistogramCombo_ = new QComboBox;
    alignmentCrystalEntry_ = new QSpinBox;
    alignmentCrystalEntry_->setRange(0, 63);
    alignmentCrystalEntry_->setValue(1);
    auto* align = new QPushButton("Show aligned spectra");
    connect(align, &QPushButton::clicked, this, [this] { ShowSpectrumAlignment(); });
    alignmentLayout->addWidget(alignmentHistogramCombo_);
    alignmentLayout->addWidget(Row({new QLabel("Target crystal"), alignmentCrystalEntry_}));
    alignmentLayout->addWidget(align);
    layout->addWidget(alignment);

    layout->addWidget(new QLabel("Calibration results"));
    resultList_ = new QListWidget;
    resultList_->setMinimumHeight(135);
    layout->addWidget(resultList_);
    auto* showSpectrum = new QPushButton("Show spectrum");
    auto* showFit = new QPushButton("Fit + residuals");
    connect(showSpectrum, &QPushButton::clicked, this, [this] { ShowSelectedResultSpectrum(); });
    connect(showFit, &QPushButton::clicked, this, [this] { ShowSelectedCalibration(); });
    layout->addWidget(Row({showSpectrum, showFit}));

    manualCorrectionGroup_ = new QGroupBox("Manual correction — select one calibration result");
    auto* manualLayout = new QVBoxLayout(manualCorrectionGroup_);
    manualLayout->addWidget(Hint("Corrections below apply only to the single crystal selected in Calibration results."));
    manualHistogramCombo_ = new QComboBox;
    manualSourceCombo_ = new QComboBox;
    manualEnergyList_ = new QListWidget;
    manualEnergyList_->setMaximumHeight(90);
    manualPeakList_ = new QListWidget;
    manualPeakList_->setMaximumHeight(90);
    manualLayout->addWidget(manualHistogramCombo_);
    manualLayout->addWidget(manualSourceCombo_);
    manualLayout->addWidget(manualEnergyList_);
    manualLayout->addWidget(new QLabel("Manual peak overrides"));
    manualLayout->addWidget(manualPeakList_);
    auto* remove = new QPushButton("Remove point");
    auto* refit = new QPushButton("Refit crystal");
    connect(remove, &QPushButton::clicked, this, [this] {
        const int selected = manualPeakList_->currentRow();
        const int crystal = CurrentResultCrystal();
        std::vector<std::size_t> indices;
        for (std::size_t i = 0; i < manualPeaks_.size(); ++i) {
            if (manualPeaks_[i].crystal == crystal) indices.push_back(i);
        }
        if (selected < 0 || selected >= static_cast<int>(indices.size())) return;
        manualPeaks_.erase(manualPeaks_.begin() + static_cast<std::ptrdiff_t>(indices[selected]));
        RefreshManualPeakList();
        RedrawDisplayedSpectrum();
    });
    connect(refit, &QPushButton::clicked, this, [this] { RefitSelectedCrystal(); });
    manualLayout->addWidget(Row({remove, refit}));
    manualCorrectionGroup_->setEnabled(false);
    layout->addWidget(manualCorrectionGroup_);
    auto* exportButton = new QPushButton("Export results and residuals to CSV...");
    connect(exportButton, &QPushButton::clicked, this, [this] { ExportCsv(); });
    layout->addWidget(exportButton);
    layout->addStretch();
    return ScrollPage(page);
}

void MainWindow::ConnectActions() {
    connect(mouseModeCombo_, &QComboBox::currentIndexChanged, this,
            [this] { UpdateInteractionMode(); });
    primaryPlot_->SetRangeSelectedCallback([this](double charge) { HandleRangeClick(charge); });
    connect(referenceSourceCombo_, &QComboBox::currentIndexChanged, this, [this] {
        if (!updatingWidgets_) RefreshEnergyList(energyList_, referenceSourceCombo_, referenceEnergyIndices_);
    });
    connect(manualSourceCombo_, &QComboBox::currentIndexChanged, this, [this] {
        if (!updatingWidgets_) RefreshEnergyList(manualEnergyList_, manualSourceCombo_, manualEnergyIndices_);
    });
    connect(referenceHistogramCombo_, &QComboBox::currentIndexChanged, this, [this] {
        if (!updatingWidgets_) ShowReferenceSpectrum();
    });
    connect(manualHistogramCombo_, &QComboBox::currentIndexChanged, this, [this] {
        if (updatingWidgets_) return;
        const auto* descriptor = DescriptorForCombo(manualHistogramCombo_);
        const int crystal = CurrentResultCrystal();
        if (descriptor && crystal >= 0) ShowCrystalSpectrum(crystal, *descriptor);
    });
    connect(resultList_, &QListWidget::currentRowChanged, this, [this] {
        if (updatingWidgets_) return;
        UpdateManualCorrectionForSelection();
        ShowSelectedResultSpectrum();
    });
    UpdateInteractionMode();
}

void MainWindow::PopulateEnergyLines() {
    if (energyLines_.empty()) {
        energyLines_ = {
            {1173.228, "Co-60 line", "Co-60"}, {1332.492, "Co-60 line", "Co-60"},
            {846.771, "Co-56 line", "Co-56"}, {1037.840, "Co-56 line", "Co-56"},
            {1238.282, "Co-56 line", "Co-56"}, {1771.351, "Co-56 line", "Co-56"},
            {2034.755, "Co-56 line", "Co-56"}, {2598.459, "Co-56 line", "Co-56"},
            {3201.962, "Co-56 line", "Co-56"}, {3253.416, "Co-56 line", "Co-56"},
            {661.657, "Cs-137 line", "Cs-137"}, {511.000, "annihilation", "Na-22"},
            {1274.537, "Na-22 line", "Na-22"},
            {1460.822, "K-40 background", "Background / contaminants"},
            {2614.511, "Tl-208 background", "Background / contaminants"}
        };
    }
    energySources_ = {"Co-60", "Co-56", "Cs-137", "Na-22",
                      "Background / contaminants", "Custom"};
    const QString previousReference = referenceSourceCombo_->currentText();
    const QString previousManual = manualSourceCombo_->currentText();
    updatingWidgets_ = true;
    referenceSourceCombo_->clear();
    manualSourceCombo_->clear();
    for (const auto& source : energySources_) {
        referenceSourceCombo_->addItem(Text(source));
        manualSourceCombo_->addItem(Text(source));
    }
    referenceSourceCombo_->setCurrentText(previousReference.isEmpty() ? "Co-60" : previousReference);
    manualSourceCombo_->setCurrentText(previousManual.isEmpty() ? "Co-60" : previousManual);
    updatingWidgets_ = false;
    RefreshEnergyList(energyList_, referenceSourceCombo_, referenceEnergyIndices_);
    RefreshEnergyList(manualEnergyList_, manualSourceCombo_, manualEnergyIndices_);
}

void MainWindow::RefreshEnergyList(QListWidget* list, const QComboBox* sourceCombo,
                                   std::vector<std::size_t>& indices) {
    list->clear();
    indices.clear();
    const std::string source = sourceCombo->currentText().toStdString();
    for (std::size_t index = 0; index < energyLines_.size(); ++index) {
        const auto& line = energyLines_[index];
        if (line.source != source) continue;
        indices.push_back(index);
        list->addItem(Text(FormatNumber(line.energy, 3) + " keV  —  " + line.label));
    }
    if (list->count() > 0) list->setCurrentRow(0);
}

void MainWindow::AddRootFiles() {
    const QStringList files = QFileDialog::getOpenFileNames(
        this, "Add ROOT files", {}, "ROOT files (*.root);;All files (*)");
    if (files.empty()) return;
    std::vector<std::string> paths;
    paths.reserve(static_cast<std::size_t>(files.size()));
    for (const QString& file : files) paths.push_back(file.toStdString());
    OpenRootFiles(paths);
}

bool MainWindow::OpenRootFiles(const std::vector<std::string>& files) {
    if (files.empty()) return false;
    std::vector<std::string> selectedIds;
    for (int index : SelectedDescriptorIndices()) selectedIds.push_back(descriptors_[index].id);
    int added = 0;
    std::string lastError;
    for (const std::string& file : files) {
        std::string error;
        auto found = repository_.Discover(file, error);
        if (!error.empty()) lastError = error;
        for (auto& descriptor : found) {
            const auto duplicate = std::find_if(descriptors_.begin(), descriptors_.end(),
                [&](const HistogramDescriptor& existing) { return existing.id == descriptor.id; });
            if (duplicate == descriptors_.end()) {
                descriptors_.push_back(std::move(descriptor));
                ++added;
            }
        }
    }
    RefreshDatasetWidgets();
    for (int row = 0; row < static_cast<int>(descriptors_.size()); ++row) {
        const bool old = std::find(selectedIds.begin(), selectedIds.end(), descriptors_[row].id) != selectedIds.end();
        const bool fresh = row >= static_cast<int>(descriptors_.size()) - added;
        histogramList_->item(row)->setSelected(old || fresh);
    }
    if (added > 0) {
        SetStatus("Discovered " + std::to_string(added) + " new TH2 histogram(s) from " +
                  std::to_string(files.size()) + " file(s).");
        ShowReferenceSpectrum();
    } else {
        SetStatus(lastError.empty() ? "The selected files contain no new TH2 histograms." : lastError);
    }
    return added > 0;
}

void MainWindow::RefreshDatasetWidgets() {
    updatingWidgets_ = true;
    histogramList_->clear();
    referenceHistogramCombo_->clear();
    manualHistogramCombo_->clear();
    alignmentHistogramCombo_->clear();
    for (std::size_t i = 0; i < descriptors_.size(); ++i) {
        const auto& descriptor = descriptors_[i];
        histogramList_->addItem(Text(descriptor.displayName + "  [" +
                                     std::to_string(descriptor.xBins) + " x " +
                                     std::to_string(descriptor.yBins) + "]"));
        const QString name = Text(descriptor.displayName);
        referenceHistogramCombo_->addItem(name, static_cast<int>(i));
        manualHistogramCombo_->addItem(name, static_cast<int>(i));
        alignmentHistogramCombo_->addItem(name, static_cast<int>(i));
    }
    updatingWidgets_ = false;
}

std::vector<int> MainWindow::SelectedDescriptorIndices() const {
    std::vector<int> indices;
    for (const auto* item : histogramList_->selectedItems()) indices.push_back(histogramList_->row(item));
    std::sort(indices.begin(), indices.end());
    return indices;
}

std::vector<int> MainWindow::SelectedCrystals() const {
    std::vector<int> crystals;
    for (const auto* item : crystalList_->selectedItems()) crystals.push_back(item->data(Qt::UserRole).toInt());
    std::sort(crystals.begin(), crystals.end());
    return crystals;
}

AxisOrientation MainWindow::Orientation() const {
    return orientationCombo_->currentIndex() == 1 ? AxisOrientation::ChargeOnY : AxisOrientation::ChargeOnX;
}

int MainWindow::ReferenceCrystal() const { return referenceCrystalEntry_->value(); }

int MainWindow::CurrentResultCrystal() const {
    const auto* item = resultList_->currentItem();
    return item ? item->data(Qt::UserRole).toInt() : -1;
}

const HistogramDescriptor* MainWindow::DescriptorForCombo(const QComboBox* combo) const {
    const int index = combo->currentData().toInt();
    return index >= 0 && index < static_cast<int>(descriptors_.size()) ? &descriptors_[index] : nullptr;
}

const MainWindow::EnergyLine* MainWindow::SelectedEnergy(const QListWidget* list) const {
    const int row = list->currentRow();
    const auto& indices = list == energyList_ ? referenceEnergyIndices_ : manualEnergyIndices_;
    return row >= 0 && row < static_cast<int>(indices.size()) ? &energyLines_[indices[row]] : nullptr;
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
    displayedSpectrum_ = std::move(spectrum);
    displayedDatasetId_ = descriptor.id;
    displayedCrystal_ = crystal;
    pendingRangeStart_.reset();
    SetSecondaryPlotVisible(false);
    RedrawDisplayedSpectrum(false);
    SetStatus("Showing crystal " + std::to_string(crystal) +
              ". Select range uses two boundary clicks; Zoom / pan uses wheel, drag, or double-click reset.");
}

void MainWindow::ShowSelectedResultSpectrum() {
    const int crystal = CurrentResultCrystal();
    const auto result = results_.find(crystal);
    if (crystal < 0 || result == results_.end()) {
        SetStatus("Select one calibration result first.");
        return;
    }

    const HistogramDescriptor* descriptor = DescriptorForCombo(manualHistogramCombo_);
    const auto fitCountFor = [&](const std::string& datasetId) {
        return static_cast<int>(std::count_if(result->second.points.begin(), result->second.points.end(),
            [&](const CalibrationPoint& point) {
                return point.datasetId == datasetId && point.peakFit.success;
            }));
    };

    // If the currently chosen source has no stored fit for this result, move to
    // the first source that does. This makes Show spectrum reliably restore the
    // red fitted peaks instead of opening an unrelated, overlay-free source.
    if ((!descriptor || fitCountFor(descriptor->id) == 0) && !result->second.points.empty()) {
        const std::string& datasetId = result->second.points.front().datasetId;
        const auto found = std::find_if(descriptors_.begin(), descriptors_.end(),
            [&](const HistogramDescriptor& item) { return item.id == datasetId; });
        if (found != descriptors_.end()) {
            const int index = static_cast<int>(std::distance(descriptors_.begin(), found));
            updatingWidgets_ = true;
            manualHistogramCombo_->setCurrentIndex(index);
            updatingWidgets_ = false;
            descriptor = &descriptors_[index];
        }
    }
    if (!descriptor) {
        SetStatus("No source histogram is available for the selected result.");
        return;
    }
    ShowCrystalSpectrum(crystal, *descriptor);
    const int fittedPeaks = fitCountFor(descriptor->id);
    SetStatus("Showing selected result Crystal " + std::to_string(crystal) + " with " +
              std::to_string(fittedPeaks) + " stored fitted peak(s) in red.");
}

void MainWindow::UpdateManualCorrectionForSelection() {
    const int crystal = CurrentResultCrystal();
    const bool selected = crystal >= 0 && results_.count(crystal) != 0;
    manualCorrectionGroup_->setEnabled(selected);
    manualCorrectionGroup_->setTitle(selected
        ? QString("Manual correction — Crystal %1 only").arg(crystal, 2, 10, QChar('0'))
        : QString("Manual correction — select one calibration result"));
    pendingRangeStart_.reset();
    RefreshManualPeakList();
}

void MainWindow::RedrawDisplayedSpectrum(bool preserveView) {
    if (!displayedSpectrum_) return;
    std::vector<PlotSeries> series{HistogramSeries(*displayedSpectrum_, QColor("#2563eb"), "Spectrum")};
    std::vector<PlotMarker> markers;
    const auto addFit = [&](const PeakFitResult& fit, double energy, bool manual) {
        if (!fit.success) return;
        series.push_back(FitSeries(fit, FormatNumber(energy, 1) + " keV", manual));
        markers.push_back({fit.centroid, FormatNumber(energy, 1) + " keV", QColor("#dc2626"), manual});
    };
    if (displayedCrystal_ == ReferenceCrystal()) {
        for (const auto& peak : referencePeaks_) {
            if (peak.datasetId == displayedDatasetId_) addFit(peak.peakFit, peak.energy, false);
        }
    }
    const auto result = results_.find(displayedCrystal_);
    if (result != results_.end()) {
        for (const auto& point : result->second.points) {
            if (point.datasetId == displayedDatasetId_) addFit(point.peakFit, point.energy, point.manual);
        }
    }
    for (const auto& peak : manualPeaks_) {
        if (peak.datasetId != displayedDatasetId_ || peak.crystal != displayedCrystal_) continue;
        const bool applied = result != results_.end() &&
            std::any_of(result->second.points.begin(), result->second.points.end(), [&](const CalibrationPoint& point) {
                return point.manual && point.datasetId == peak.datasetId &&
                       std::abs(point.energy - peak.energy) < 1e-6 &&
                       std::abs(point.charge - peak.peakFit.centroid) < 1e-6;
            });
        if (!applied) addFit(peak.peakFit, peak.energy, true);
    }
    if (pendingRangeStart_) {
        markers.push_back({*pendingRangeStart_, "first fit limit", QColor("#d97706"), true});
    }
    primaryPlot_->SetPlot("Crystal " + std::to_string(displayedCrystal_) + " spectrum",
                          "Charge", "Counts", std::move(series), std::move(markers),
                          {preserveView, 0.0});
}

void MainWindow::UpdateInteractionMode() {
    const bool selecting = mouseModeCombo_->currentIndex() == 0;
    primaryPlot_->SetInteractionMode(selecting ? SpectrumPlotWidget::InteractionMode::SelectRange
                                               : SpectrumPlotWidget::InteractionMode::ZoomPan);
    secondaryPlot_->SetInteractionMode(SpectrumPlotWidget::InteractionMode::ZoomPan);
    if (!selecting) pendingRangeStart_.reset();
    if (displayedSpectrum_) RedrawDisplayedSpectrum();
    SetStatus(selecting ? "Peak-range mode: click two fit limits; wheel zoom and right-drag pan remain available."
                        : "Zoom mode: wheel or buttons zoom, left-drag a window, right-drag pans, double-click resets.");
}

void MainWindow::HandleRangeClick(double charge) {
    if (mouseModeCombo_->currentIndex() != 0 || !displayedSpectrum_ ||
        (tabs_->currentIndex() != 1 && tabs_->currentIndex() != 2)) return;
    if (!pendingRangeStart_) {
        pendingRangeStart_ = charge;
        RedrawDisplayedSpectrum();
        SetStatus("Lower peak-fit limit selected at " + FormatNumber(charge, 3) + ". Click the upper limit.");
        return;
    }
    const double low = std::min(*pendingRangeStart_, charge);
    const double high = std::max(*pendingRangeStart_, charge);
    pendingRangeStart_.reset();
    const double minimumWidth = 12.0 * displayedSpectrum_->GetXaxis()->GetBinWidth(1);
    if (high - low < minimumWidth) {
        RedrawDisplayedSpectrum();
        SetStatus("Peak-fit interval is too narrow; select at least 12 histogram bins.");
        return;
    }
    const auto fit = CalibrationEngine::FitRadwarePeak(*displayedSpectrum_, low, high);
    if (!fit.success) {
        RedrawDisplayedSpectrum();
        SetStatus("Peak fit failed: " + fit.status);
        return;
    }
    if (tabs_->currentIndex() == 1) AddReferencePeak(fit);
    else AddManualPeak(fit);
}

void MainWindow::AddReferencePeak(const PeakFitResult& fit) {
    const auto* energy = SelectedEnergy(energyList_);
    const auto* descriptor = DescriptorForCombo(referenceHistogramCombo_);
    if (!energy || !descriptor || displayedDatasetId_ != descriptor->id || displayedCrystal_ != ReferenceCrystal()) {
        SetStatus("Show the selected histogram's reference spectrum and choose an energy first.");
        RedrawDisplayedSpectrum();
        return;
    }
    auto duplicate = std::find_if(referencePeaks_.begin(), referencePeaks_.end(), [&](const ReferencePeak& peak) {
        return peak.datasetId == descriptor->id && std::abs(peak.energy - energy->energy) < 1e-6;
    });
    ReferencePeak peak{descriptor->id, fit.centroid, energy->energy,
                       energy->source + " / " + energy->label, fit};
    if (duplicate == referencePeaks_.end()) referencePeaks_.push_back(std::move(peak));
    else *duplicate = std::move(peak);
    RefreshReferencePeakList();
    RedrawDisplayedSpectrum();
    SetStatus("RadWare centroid " + FormatNumber(fit.centroid, 3) + " +/- " +
              FormatNumber(fit.centroidError, 3) + " charge -> " +
              FormatNumber(energy->energy, 3) + " keV.");
}

void MainWindow::AddManualPeak(const PeakFitResult& fit) {
    const auto* energy = SelectedEnergy(manualEnergyList_);
    const auto* descriptor = DescriptorForCombo(manualHistogramCombo_);
    const int crystal = CurrentResultCrystal();
    if (!energy || !descriptor || crystal < 0 || displayedDatasetId_ != descriptor->id || displayedCrystal_ != crystal) {
        SetStatus("Select a result, histogram and energy, then show that crystal's spectrum.");
        RedrawDisplayedSpectrum();
        return;
    }
    auto duplicate = std::find_if(manualPeaks_.begin(), manualPeaks_.end(), [&](const ManualPeak& peak) {
        return peak.datasetId == descriptor->id && peak.crystal == crystal &&
               std::abs(peak.energy - energy->energy) < 1e-6;
    });
    ManualPeak peak{descriptor->id, crystal, fit.centroid, energy->energy,
                    energy->source + " / " + energy->label, fit};
    if (duplicate == manualPeaks_.end()) manualPeaks_.push_back(std::move(peak));
    else *duplicate = std::move(peak);
    RefreshManualPeakList();
    RedrawDisplayedSpectrum();
    SetStatus("Manual RadWare centroid fitted. Click Refit crystal to apply it.");
}

void MainWindow::RefreshReferencePeakList() {
    referencePeakList_->clear();
    for (const auto& peak : referencePeaks_) {
        const auto descriptor = std::find_if(descriptors_.begin(), descriptors_.end(), [&](const HistogramDescriptor& item) {
            return item.id == peak.datasetId;
        });
        const std::string dataset = descriptor == descriptors_.end() ? "missing" : descriptor->displayName;
        referencePeakList_->addItem(Text(dataset + " | [" + FormatNumber(peak.peakFit.rangeLow, 2) +
            ", " + FormatNumber(peak.peakFit.rangeHigh, 2) + "] centroid " +
            FormatNumber(peak.charge, 3) + " +/- " + FormatNumber(peak.peakFit.centroidError, 3) +
            " -> " + FormatNumber(peak.energy, 3) + " keV (" + peak.label + ")"));
    }
}

void MainWindow::RefreshManualPeakList() {
    manualPeakList_->clear();
    const int crystal = CurrentResultCrystal();
    for (const auto& peak : manualPeaks_) {
        if (peak.crystal != crystal) continue;
        manualPeakList_->addItem(Text("[" + FormatNumber(peak.peakFit.rangeLow, 2) + ", " +
            FormatNumber(peak.peakFit.rangeHigh, 2) + "] centroid " + FormatNumber(peak.charge, 3) +
            " -> " + FormatNumber(peak.energy, 3) + " keV (" + peak.label + ")"));
    }
}

std::vector<CalibrationPoint> MainWindow::BuildPointsForCrystal(int crystal) {
    std::vector<CalibrationPoint> points;
    CalibrationEngine::SearchOptions options;
    options.sigmaBins = sigmaEntry_->value();
    options.threshold = thresholdEntry_->value();
    for (int descriptorIndex : SelectedDescriptorIndices()) {
        const auto& descriptor = descriptors_[descriptorIndex];
        std::vector<ReferencePeak> references;
        for (const auto& peak : referencePeaks_) if (peak.datasetId == descriptor.id) references.push_back(peak);
        if (references.empty()) continue;
        std::sort(references.begin(), references.end(), [](const ReferencePeak& a, const ReferencePeak& b) {
            return a.charge < b.charge;
        });
        if (crystal == ReferenceCrystal()) {
            for (const auto& ref : references) {
                points.push_back({descriptor.id, ref.peakFit.centroid, ref.energy,
                                  ref.peakFit.centroidError, false, 0.0, ref.peakFit});
            }
            continue;
        }
        std::string error;
        auto referenceSpectrum = repository_.ProjectCrystal(descriptor, ReferenceCrystal(), Orientation(), error);
        auto spectrum = repository_.ProjectCrystal(descriptor, crystal, Orientation(), error);
        if (!referenceSpectrum || !spectrum) continue;
        const auto alignment = CalibrationEngine::AlignSpectrumPatterns(*referenceSpectrum, *spectrum, options);
        if (!alignment.success) continue;
        for (const auto& reference : references) {
            const double low = CalibrationEngine::MapReferenceCharge(alignment, reference.peakFit.rangeLow);
            const double high = CalibrationEngine::MapReferenceCharge(alignment, reference.peakFit.rangeHigh);
            const auto fit = CalibrationEngine::FitRadwarePeak(*spectrum, low, high);
            if (fit.success) {
                points.push_back({descriptor.id, fit.centroid, reference.energy,
                                  fit.centroidError, false, 0.0, fit});
            }
        }
    }
    for (const auto& manual : manualPeaks_) {
        if (manual.crystal != crystal) continue;
        auto existing = std::find_if(points.begin(), points.end(), [&](const CalibrationPoint& point) {
            return point.datasetId == manual.datasetId && std::abs(point.energy - manual.energy) < 1e-6;
        });
        CalibrationPoint replacement{manual.datasetId, manual.peakFit.centroid, manual.energy,
                                     manual.peakFit.centroidError, true, 0.0, manual.peakFit};
        if (existing == points.end()) points.push_back(replacement);
        else *existing = replacement;
    }
    return points;
}

CalibrationResult MainWindow::CalibrateCrystal(int crystal) {
    return CalibrationEngine::FitSecondOrder(crystal, BuildPointsForCrystal(crystal), residualLimitEntry_->value());
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
    QApplication::processEvents();
    results_.clear();
    for (int crystal : crystals) {
        results_[crystal] = CalibrateCrystal(crystal);
        QApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
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

void MainWindow::ShowSpectrumAlignment() {
    const auto* descriptor = DescriptorForCombo(alignmentHistogramCombo_);
    if (!descriptor) {
        SetStatus("Choose a source histogram for the alignment preview.");
        return;
    }
    const int referenceCrystal = ReferenceCrystal();
    const int targetCrystal = alignmentCrystalEntry_->value();
    std::string error;
    auto reference = repository_.ProjectCrystal(*descriptor, referenceCrystal, Orientation(), error);
    auto target = repository_.ProjectCrystal(*descriptor, targetCrystal, Orientation(), error);
    if (!reference || !target) {
        SetStatus(error);
        return;
    }
    CalibrationEngine::SearchOptions options;
    options.sigmaBins = sigmaEntry_->value();
    options.threshold = thresholdEntry_->value();
    const auto match = CalibrationEngine::AlignSpectrumPatterns(*reference, *target, options);
    if (!match.success || !(match.scale > 0.0) || !std::isfinite(match.scale) || !std::isfinite(match.offset)) {
        SetStatus("Could not align this spectrum: not enough corresponding peaks were found.");
        return;
    }
    const double referenceMaximum = std::max(reference->GetMaximum(), 1.0);
    const double targetMaximum = std::max(target->GetMaximum(), 1.0);
    PlotSeries referenceSeries = HistogramSeries(*reference, QColor("#2563eb"),
                                                  "Reference C" + std::to_string(referenceCrystal));
    for (double& value : referenceSeries.y) value /= referenceMaximum;
    PlotSeries targetSeries;
    targetSeries.name = "Aligned C" + std::to_string(targetCrystal);
    targetSeries.color = QColor("#dc2626");
    targetSeries.width = 1;
    for (int bin = 1; bin <= target->GetNbinsX(); ++bin) {
        const double charge = target->GetXaxis()->GetBinCenter(bin);
        targetSeries.x.push_back(CalibrationEngine::MapTargetChargeToReference(match, charge));
        targetSeries.y.push_back(target->GetBinContent(bin) / targetMaximum);
    }
    std::vector<PlotMarker> markers;
    for (std::size_t i = 0; i < match.matched.size(); ++i) {
        if (match.matched[i]) markers.push_back({match.referenceCharges[i], {}, QColor("#16a34a"), true});
    }
    displayedSpectrum_.reset();
    displayedDatasetId_.clear();
    displayedCrystal_ = -1;
    pendingRangeStart_.reset();
    SetSecondaryPlotVisible(false);
    primaryPlot_->SetPlot("Pre-calibration alignment", "Reference-spectrum charge", "Normalized counts",
                          {std::move(referenceSeries), std::move(targetSeries)}, std::move(markers));
    const int matched = static_cast<int>(std::count(match.matched.begin(), match.matched.end(), true));
    SetStatus("Alignment preview: " + std::to_string(matched) + "/" +
              std::to_string(match.referenceCharges.size()) + " pattern peaks; target charge = " +
              FormatNumber(match.offset, 3) + " + " + FormatNumber(match.scale, 6) + " q + " +
              FormatNumber(match.quadratic, 9) + " q^2. No energy calibration applied.");
}

void MainWindow::RefitSelectedCrystal() {
    const int crystal = CurrentResultCrystal();
    if (crystal < 0) {
        SetStatus("Select a calibration result to refit.");
        return;
    }
    results_[crystal] = CalibrateCrystal(crystal);
    RefreshResults();
    for (int row = 0; row < resultList_->count(); ++row) {
        if (resultList_->item(row)->data(Qt::UserRole).toInt() == crystal) {
            resultList_->setCurrentRow(row);
            break;
        }
    }
    const auto* descriptor = DescriptorForCombo(manualHistogramCombo_);
    if (descriptor) ShowCrystalSpectrum(crystal, *descriptor);
    SetStatus("Crystal " + std::to_string(crystal) + " refitted with manual overrides.");
}

void MainWindow::RefreshResults() {
    const int previous = CurrentResultCrystal();
    updatingWidgets_ = true;
    resultList_->clear();
    int selectedRow = -1;
    for (const auto& [crystal, result] : results_) {
        const std::string state = !result.success ? "FAIL" : result.needsReview ? "REVIEW" : "OK";
        std::ostringstream row;
        row << '[' << state << "] C" << std::setw(2) << std::setfill('0') << crystal
            << " | n=" << result.points.size() << " | RMS=" << std::fixed << std::setprecision(3)
            << result.residualRms << " | " << std::scientific << std::setprecision(4)
            << result.p0 << ", " << result.p1 << ", " << result.p2;
        auto* item = new QListWidgetItem(Text(row.str()), resultList_);
        item->setData(Qt::UserRole, crystal);
        if (!result.success) item->setForeground(QColor("#b91c1c"));
        else if (result.needsReview) item->setForeground(QColor("#b45309"));
        else item->setForeground(QColor("#15803d"));
        if (crystal == previous) selectedRow = resultList_->count() - 1;
    }
    if (selectedRow < 0 && resultList_->count() > 0) selectedRow = 0;
    resultList_->setCurrentRow(selectedRow);
    updatingWidgets_ = false;
    UpdateManualCorrectionForSelection();
}

void MainWindow::ShowSelectedCalibration() {
    const int crystal = CurrentResultCrystal();
    const auto found = results_.find(crystal);
    if (found == results_.end() || !found->second.success) {
        SetStatus(found == results_.end() ? "Select a result first." : found->second.status);
        return;
    }
    const auto& result = found->second;
    PlotSeries points;
    points.name = "Peak centroids";
    points.color = QColor("#2563eb");
    points.points = true;
    PlotSeries residuals;
    residuals.name = "Residual";
    residuals.color = result.needsReview ? QColor("#d97706") : QColor("#16a34a");
    residuals.points = true;
    for (const auto& point : result.points) {
        points.x.push_back(point.charge);
        points.y.push_back(point.energy);
        residuals.x.push_back(point.charge);
        residuals.y.push_back(point.residual);
    }
    const auto [minimum, maximum] = std::minmax_element(points.x.begin(), points.x.end());
    const double dataWidth = std::max(*maximum - *minimum, 1.0);
    const double plotMinimum = *minimum - 0.15 * dataWidth;
    const double plotMaximum = *maximum + 0.15 * dataWidth;
    PlotSeries fit;
    fit.name = "Second-order fit";
    fit.color = QColor("#dc2626");
    fit.width = 2;
    for (int i = 0; i < 240; ++i) {
        const double charge = plotMinimum + (plotMaximum - plotMinimum) * i / 239.0;
        fit.x.push_back(charge);
        fit.y.push_back(result.p0 + result.p1 * charge + result.p2 * charge * charge);
    }
    displayedSpectrum_.reset();
    displayedDatasetId_.clear();
    displayedCrystal_ = -1;
    pendingRangeStart_.reset();
    SetSecondaryPlotVisible(true);
    primaryPlot_->SetPlot("Crystal " + std::to_string(crystal) + " calibration",
                          "Charge", "Energy (keV)", {std::move(fit), std::move(points)});
    secondaryPlot_->SetPlot("Residuals — RMS " + FormatNumber(result.residualRms, 3) + " keV",
                            "Charge", "Energy - fit (keV)", {std::move(residuals)}, {},
                            {false, 0.15});
    SetStatus("Crystal " + std::to_string(crystal) + ": " + result.status);
}

void MainWindow::ExportCsv() {
    if (results_.empty()) {
        SetStatus("There are no calibration results to export.");
        return;
    }
    const QString path = QFileDialog::getSaveFileName(this, "Export calibration",
        "hpge_calibration.csv", "CSV files (*.csv);;All files (*)");
    if (path.isEmpty()) return;
    std::ofstream output(path.toStdString());
    if (!output) {
        SetStatus("Could not write CSV file.");
        return;
    }
    output << "record,crystal,status,needs_review,p0,p1,p2,chi2,ndf,residual_rms_keV,"
              "dataset,charge,energy_keV,residual_keV,manual,charge_error,range_low,range_high,"
              "peak_sigma,peak_chi2,peak_ndf,tail_fraction,beta,step_fraction\n";
    output << std::setprecision(12);
    for (const auto& [crystal, result] : results_) {
        output << "fit," << crystal << ",\"" << result.status << "\"," << result.needsReview
               << ',' << result.p0 << ',' << result.p1 << ',' << result.p2 << ','
               << result.chi2 << ',' << result.ndf << ',' << result.residualRms << ",,,,,,,,,,,,,,\n";
        for (const auto& point : result.points) {
            output << "point," << crystal << ",,," << result.p0 << ',' << result.p1 << ','
                   << result.p2 << ",,,,\"" << point.datasetId << "\"," << point.charge << ','
                   << point.energy << ',' << point.residual << ',' << point.manual << ','
                   << point.chargeError << ',' << point.peakFit.rangeLow << ','
                   << point.peakFit.rangeHigh << ',' << point.peakFit.sigma << ','
                   << point.peakFit.chi2 << ',' << point.peakFit.ndf << ','
                   << point.peakFit.tailFraction << ',' << point.peakFit.beta << ','
                   << point.peakFit.stepFraction << '\n';
        }
    }
    SetStatus("Exported coefficients and per-peak residuals to " + path.toStdString());
}

void MainWindow::SetStatus(const std::string& text) {
    if (statusLabel_) statusLabel_->setText(Text(text));
}

void MainWindow::SetSecondaryPlotVisible(bool visible) {
    secondaryPlotContainer_->setVisible(visible);
}

} // namespace hpge
