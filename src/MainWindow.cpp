#include "MainWindow.h"

#include "CalibrationExporter.h"
#include "SpectrumPlotWidget.h"

#include <QApplication>
#include <QCheckBox>
#include <QClipboard>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QFile>
#include <QFileInfo>
#include <QFrame>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QMenuBar>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QPushButton>
#include <QSaveFile>
#include <QScrollArea>
#include <QSettings>
#include <QSpinBox>
#include <QSplitter>
#include <QStatusBar>
#include <QTabWidget>
#include <QVBoxLayout>
#include <QDir>

#include <TH1D.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <set>
#include <sstream>

namespace hpge {
namespace {

QString Text(const std::string& value) { return QString::fromStdString(value); }

std::string FormatNumber(double value, int precision = 5) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(precision) << value;
    return out.str();
}

std::string FormatPrecise(double value) {
    std::ostringstream out;
    out << std::setprecision(12) << value;
    return out.str();
}

std::string CsvField(const std::string& value) {
    if (value.find_first_of(",\"\r\n") == std::string::npos) return value;
    std::string escaped = "\"";
    for (char character : value) {
        if (character == '"') escaped += '"';
        escaped += character;
    }
    escaped += '"';
    return escaped;
}

void WriteCsvRow(std::ostream& output, const std::vector<std::string>& fields) {
    for (std::size_t index = 0; index < fields.size(); ++index) {
        if (index != 0) output << ',';
        output << CsvField(fields[index]);
    }
    output << '\n';
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

QJsonArray DoubleArray(const std::vector<double>& values) {
    QJsonArray array;
    for (double value : values) array.append(value);
    return array;
}

QJsonArray BoolArray(const std::vector<bool>& values) {
    QJsonArray array;
    for (bool value : values) array.append(value);
    return array;
}

std::vector<double> ReadDoubleArray(const QJsonValue& value) {
    std::vector<double> output;
    for (const auto item : value.toArray()) output.push_back(item.toDouble());
    return output;
}

std::vector<bool> ReadBoolArray(const QJsonValue& value) {
    std::vector<bool> output;
    for (const auto item : value.toArray()) output.push_back(item.toBool());
    return output;
}

QJsonObject PeakFitJson(const PeakFitResult& fit) {
    return {{"success", fit.success}, {"status", Text(fit.status)},
            {"rangeLow", fit.rangeLow}, {"rangeHigh", fit.rangeHigh},
            {"centroid", fit.centroid}, {"centroidError", fit.centroidError},
            {"sigma", fit.sigma}, {"height", fit.height},
            {"tailFraction", fit.tailFraction}, {"beta", fit.beta},
            {"stepFraction", fit.stepFraction}, {"background0", fit.background0},
            {"background1", fit.background1}, {"background2", fit.background2},
            {"chi2", fit.chi2}, {"ndf", fit.ndf}};
}

PeakFitResult ReadPeakFit(const QJsonValue& value) {
    const auto object = value.toObject();
    PeakFitResult fit;
    fit.success = object["success"].toBool();
    fit.status = object["status"].toString().toStdString();
    fit.rangeLow = object["rangeLow"].toDouble();
    fit.rangeHigh = object["rangeHigh"].toDouble();
    fit.centroid = object["centroid"].toDouble();
    fit.centroidError = object["centroidError"].toDouble();
    fit.sigma = object["sigma"].toDouble();
    fit.height = object["height"].toDouble();
    fit.tailFraction = object["tailFraction"].toDouble();
    fit.beta = object["beta"].toDouble();
    fit.stepFraction = object["stepFraction"].toDouble();
    fit.background0 = object["background0"].toDouble();
    fit.background1 = object["background1"].toDouble();
    fit.background2 = object["background2"].toDouble();
    fit.chi2 = object["chi2"].toDouble();
    fit.ndf = object["ndf"].toInt();
    return fit;
}

QJsonObject CalibrationPointJson(const CalibrationPoint& point) {
    return {{"datasetId", Text(point.datasetId)}, {"charge", point.charge},
            {"energy", point.energy}, {"chargeError", point.chargeError},
            {"manual", point.manual}, {"residual", point.residual},
            {"peakFit", PeakFitJson(point.peakFit)}};
}

CalibrationPoint ReadCalibrationPoint(const QJsonValue& value) {
    const auto object = value.toObject();
    CalibrationPoint point;
    point.datasetId = object["datasetId"].toString().toStdString();
    point.charge = object["charge"].toDouble();
    point.energy = object["energy"].toDouble();
    point.chargeError = object["chargeError"].toDouble();
    point.manual = object["manual"].toBool();
    point.residual = object["residual"].toDouble();
    point.peakFit = ReadPeakFit(object["peakFit"]);
    return point;
}

QJsonObject CalibrationResultJson(const CalibrationResult& result) {
    QJsonArray points;
    for (const auto& point : result.points) points.append(CalibrationPointJson(point));
    return {{"crystal", result.crystal}, {"success", result.success},
            {"needsReview", result.needsReview}, {"status", Text(result.status)},
            {"p0", result.p0}, {"p1", result.p1}, {"p2", result.p2},
            {"chi2", result.chi2}, {"ndf", result.ndf},
            {"residualRms", result.residualRms}, {"points", points}};
}

CalibrationResult ReadCalibrationResult(const QJsonValue& value) {
    const auto object = value.toObject();
    CalibrationResult result;
    result.crystal = object["crystal"].toInt(-1);
    result.success = object["success"].toBool();
    result.needsReview = object["needsReview"].toBool(true);
    result.status = object["status"].toString().toStdString();
    result.p0 = object["p0"].toDouble();
    result.p1 = object["p1"].toDouble();
    result.p2 = object["p2"].toDouble();
    result.chi2 = object["chi2"].toDouble();
    result.ndf = object["ndf"].toInt();
    result.residualRms = object["residualRms"].toDouble();
    for (const auto point : object["points"].toArray()) {
        result.points.push_back(ReadCalibrationPoint(point));
    }
    return result;
}

QJsonObject PeakMatchJson(const PeakMatchResult& match) {
    return {{"success", match.success},
            {"referenceCharges", DoubleArray(match.referenceCharges)},
            {"charges", DoubleArray(match.charges)},
            {"matched", BoolArray(match.matched)},
            {"scale", match.scale}, {"offset", match.offset},
            {"quadratic", match.quadratic}, {"score", match.score},
            {"alignmentCost", match.alignmentCost},
            {"referenceSensitivity", match.referenceSensitivity},
            {"targetSensitivity", match.targetSensitivity},
            {"quadraticModel", match.quadraticModel}};
}

PeakMatchResult ReadPeakMatch(const QJsonValue& value) {
    const auto object = value.toObject();
    PeakMatchResult match;
    match.success = object["success"].toBool();
    match.referenceCharges = ReadDoubleArray(object["referenceCharges"]);
    match.charges = ReadDoubleArray(object["charges"]);
    match.matched = ReadBoolArray(object["matched"]);
    match.scale = object["scale"].toDouble(1.0);
    match.offset = object["offset"].toDouble();
    match.quadratic = object["quadratic"].toDouble();
    match.score = object["score"].toDouble();
    match.alignmentCost = object["alignmentCost"].toDouble();
    match.referenceSensitivity = object["referenceSensitivity"].toDouble();
    match.targetSensitivity = object["targetSensitivity"].toDouble();
    match.quadraticModel = object["quadraticModel"].toBool();
    return match;
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
    auto* openProjectAction = fileMenu->addAction("&Open project...");
    openProjectAction->setShortcut(QKeySequence("Ctrl+Shift+O"));
    auto* saveProjectAction = fileMenu->addAction("&Save project...");
    saveProjectAction->setShortcut(QKeySequence("Ctrl+Shift+S"));
    fileMenu->addSeparator();
    auto* openAction = fileMenu->addAction("&Add ROOT files...");
    openAction->setShortcut(QKeySequence::Open);
    auto* exportAction = fileMenu->addAction("Export calibration CSV...");
    fileMenu->addSeparator();
    auto* quitAction = fileMenu->addAction("Quit");
    connect(openProjectAction, &QAction::triggered, this, [this] { OpenProjectDialog(); });
    connect(saveProjectAction, &QAction::triggered, this, [this] { SaveProjectDialog(); });
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
    primaryPlot_->setObjectName("primaryPlot");
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
    statusLabel_->setObjectName("statusLabel");
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
    referenceSourceCombo_->setObjectName("referenceSourceCombo");
    layout->addWidget(referenceSourceCombo_);
    layout->addWidget(new QLabel("Known line for the next fit"));
    energyList_ = new QListWidget;
    energyList_->setObjectName("referenceEnergyList");
    energyList_->setMinimumHeight(150);
    layout->addWidget(energyList_);

    layout->addWidget(new QLabel("Custom line"));
    customEnergyEntry_ = RealEntry(0.0, 3, 0.0, 100000.0);
    customEnergyEntry_->setObjectName("customEnergyEntry");
    customEnergyEntry_->setSuffix(" keV");
    customLabelEntry_ = new QLineEdit;
    customLabelEntry_->setObjectName("customLineLabelEntry");
    customLabelEntry_->setPlaceholderText("Optional label");
    auto* addCustom = new QPushButton("Add for session");
    addCustom->setObjectName("addCustomLineButton");
    auto* saveCustom = new QPushButton("Save for future");
    saveCustom->setObjectName("saveCustomLineButton");
    auto* removeCustom = new QPushButton("Remove selected custom line");
    removeCustom->setObjectName("removeCustomLineButton");
    connect(addCustom, &QPushButton::clicked, this,
            [this] { AddCustomLine(false); });
    connect(saveCustom, &QPushButton::clicked, this,
            [this] { AddCustomLine(true); });
    connect(removeCustom, &QPushButton::clicked, this,
            [this] { RemoveSelectedCustomLine(); });
    layout->addWidget(Row({customEnergyEntry_, customLabelEntry_}));
    layout->addWidget(Row({addCustom, saveCustom}));
    layout->addWidget(removeCustom);
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
    auto* settings = new QGroupBox("Automatic peak fitting");
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
    alignmentLayout->addWidget(Hint("Alignment searches each projected crystal independently. Auto-tuning adjusts the requested sensitivity for that spectrum's count level and spike content; intense low-energy X-rays are compressed during peak discovery so they cannot hide weaker gamma lines."));
    alignmentSensitivityEntry_ = RealEntry(35.0, 0, 0.0, 100.0);
    alignmentSensitivityEntry_->setSuffix(" %");
    alignmentSensitivityEntry_->setToolTip(
        "Continuous starting sensitivity: lower values reject more narrow/noisy candidates; "
        "higher values retain weaker peaks.");
    alignedFitHalfRangeEntry_ = RealEntry(35.0, 1, 1.0, 1000.0);
    alignedFitHalfRangeEntry_->setObjectName("alignedFitHalfRangeEntry");
    alignedFitHalfRangeEntry_->setSuffix(" charge");
    alignedFitHalfRangeEntry_->setToolTip(
        "Half-width in aligned reference-spectrum charge. The exact interval is mapped "
        "into the target spectrum before the RadWare centroid fit.");
    autoTuneAlignmentEntry_ = new QCheckBox("Auto-tune for each crystal spectrum");
    autoTuneAlignmentEntry_->setChecked(true);
    alignmentModelCombo_ = new QComboBox;
    alignmentModelCombo_->addItem("Auto — affine or 2nd order");
    alignmentModelCombo_->addItem("Affine: a0 + a1 q");
    alignmentModelCombo_->addItem("2nd-order polynomial: a0 + a1 q + a2 q²");
    alignmentHistogramCombo_ = new QComboBox;
    alignmentHistogramCombo_->setObjectName("alignmentHistogramCombo");
    alignmentCrystalEntry_ = new QSpinBox;
    alignmentCrystalEntry_->setRange(0, 63);
    alignmentCrystalEntry_->setValue(1);
    auto* align = new QPushButton("Show aligned spectra");
    align->setObjectName("showAlignmentButton");
    connect(align, &QPushButton::clicked, this, [this] { ShowSpectrumAlignment(); });
    alignmentLayout->addWidget(Row({new QLabel("Peak sensitivity"),
                                    alignmentSensitivityEntry_}));
    alignmentLayout->addWidget(Row({new QLabel("Aligned peak-fit half-range"),
                                    alignedFitHalfRangeEntry_}));
    alignmentLayout->addWidget(autoTuneAlignmentEntry_);
    alignmentLayout->addWidget(Row({new QLabel("Charge mapping"), alignmentModelCombo_}));
    alignmentLayout->addWidget(alignmentHistogramCombo_);
    alignmentLayout->addWidget(Row({new QLabel("Target crystal"), alignmentCrystalEntry_}));
    alignmentLayout->addWidget(align);
    alignmentPreviewParameters_ = new QLabel("No alignment has been previewed.");
    alignmentPreviewParameters_->setObjectName("alignmentPreviewParameters");
    alignmentPreviewParameters_->setWordWrap(true);
    alignmentPreviewParameters_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    alignmentLayout->addWidget(alignmentPreviewParameters_);
    layout->addWidget(alignment);

    layout->addWidget(new QLabel("Calibration results"));
    resultList_ = new QListWidget;
    resultList_->setMinimumHeight(135);
    layout->addWidget(resultList_);
    layout->addWidget(new QLabel("Alignment parameters for selected crystal"));
    alignmentParameterList_ = new QListWidget;
    alignmentParameterList_->setObjectName("alignmentParameterList");
    alignmentParameterList_->setMaximumHeight(105);
    layout->addWidget(alignmentParameterList_);
    auto* copyAlignment = new QPushButton("Copy selected alignment parameters");
    copyAlignment->setObjectName("copyAlignmentParametersButton");
    connect(copyAlignment, &QPushButton::clicked, this, [this] {
        const auto* item = alignmentParameterList_->currentItem();
        if (!item) {
            SetStatus("Select an alignment-parameter row to copy.");
            return;
        }
        QApplication::clipboard()->setText(item->text());
        SetStatus("Copied the selected alignment parameters to the clipboard.");
    });
    layout->addWidget(copyAlignment);
    layout->addWidget(new QLabel("Fitted source spectra used by selected result"));
    resultSpectrumCombo_ = new QComboBox;
    resultSpectrumCombo_->setObjectName("resultSpectrumCombo");
    layout->addWidget(resultSpectrumCombo_);
    auto* showSpectrum = new QPushButton("Show selected source");
    showSpectrum->setObjectName("showResultSpectrumButton");
    auto* showAllSpectra = new QPushButton("Show all fitted sources");
    showAllSpectra->setObjectName("showAllResultSpectraButton");
    auto* showFit = new QPushButton("Fit + residuals");
    connect(showSpectrum, &QPushButton::clicked, this, [this] { ShowSelectedResultSpectrum(); });
    connect(showAllSpectra, &QPushButton::clicked, this, [this] { ShowAllResultSpectra(); });
    connect(showFit, &QPushButton::clicked, this, [this] { ShowSelectedCalibration(); });
    layout->addWidget(Row({showSpectrum, showAllSpectra, showFit}));

    auto* residualOverview = new QGroupBox("Residual overview across detector crystals");
    auto* residualOverviewLayout = new QVBoxLayout(residualOverview);
    residualOverviewLayout->addWidget(Hint(
        "Plots every fitted energy residual versus crystal 0-63 and summarizes each "
        "crystal's RMS and maximum absolute residual."));
    residualDatasetCombo_ = new QComboBox;
    residualDatasetCombo_->setObjectName("residualDatasetCombo");
    auto* showResidualOverview = new QPushButton("Show residuals vs crystals");
    showResidualOverview->setObjectName("showResidualsByCrystalButton");
    connect(showResidualOverview, &QPushButton::clicked, this,
            [this] { ShowResidualsByCrystal(); });
    residualOverviewLayout->addWidget(residualDatasetCombo_);
    residualOverviewLayout->addWidget(showResidualOverview);
    layout->addWidget(residualOverview);

    auto* combined = new QGroupBox("Combined calibrated spectrum quality");
    auto* combinedLayout = new QVBoxLayout(combined);
    combinedLayout->addWidget(Hint("After calibration, sum all successful crystals in energy space and refit each assigned line."));
    combinedHistogramCombo_ = new QComboBox;
    combinedQualityList_ = new QListWidget;
    combinedQualityList_->setMinimumHeight(105);
    auto* showCombined = new QPushButton("Show combined spectrum + energy residuals");
    connect(showCombined, &QPushButton::clicked, this, [this] { ShowCombinedSpectrum(); });
    combinedLayout->addWidget(combinedHistogramCombo_);
    combinedLayout->addWidget(combinedQualityList_);
    combinedLayout->addWidget(showCombined);
    layout->addWidget(combined);

    manualCorrectionGroup_ = new QGroupBox("Manual correction — select one calibration result");
    auto* manualLayout = new QVBoxLayout(manualCorrectionGroup_);
    manualLayout->addWidget(Hint("Corrections below apply only to the single crystal selected in Calibration results."));
    manualHistogramCombo_ = new QComboBox;
    manualSourceCombo_ = new QComboBox;
    manualEnergyList_ = new QListWidget;
    manualEnergyList_->setMaximumHeight(90);
    fittedPointList_ = new QListWidget;
    fittedPointList_->setObjectName("fittedPointList");
    fittedPointList_->setMaximumHeight(125);
    manualPeakList_ = new QListWidget;
    manualPeakList_->setMaximumHeight(90);
    manualLayout->addWidget(new QLabel(
        "Existing fitted peaks — select one to replace, or choose another energy to add"));
    manualLayout->addWidget(fittedPointList_);
    manualLayout->addWidget(manualHistogramCombo_);
    manualLayout->addWidget(manualSourceCombo_);
    manualLayout->addWidget(manualEnergyList_);
    manualLayout->addWidget(new QLabel("Pending peak replacements / additions"));
    manualLayout->addWidget(manualPeakList_);
    auto* remove = new QPushButton("Remove point");
    auto* refit = new QPushButton("Apply pending peaks (keep all others)");
    refit->setObjectName("applySelectiveRefitButton");
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
    auto* exportButton = new QPushButton("Export CSV + three C++ coefficient lists...");
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
    connect(fittedPointList_, &QListWidget::currentRowChanged, this, [this](int row) {
        if (updatingWidgets_ || row < 0) return;
        const int crystal = CurrentResultCrystal();
        const auto result = results_.find(crystal);
        const auto* item = fittedPointList_->item(row);
        const int pointIndex = item ? item->data(Qt::UserRole).toInt() : -1;
        if (result == results_.end() || pointIndex < 0 ||
            pointIndex >= static_cast<int>(result->second.points.size())) return;
        const auto& point = result->second.points[pointIndex];
        const auto descriptor = std::find_if(
            descriptors_.begin(), descriptors_.end(), [&](const HistogramDescriptor& value) {
                return value.id == point.datasetId;
            });
        if (descriptor == descriptors_.end()) return;
        const int descriptorIndex = static_cast<int>(std::distance(descriptors_.begin(), descriptor));
        const auto line = std::find_if(energyLines_.begin(), energyLines_.end(),
            [&](const EnergyLine& value) { return std::abs(value.energy - point.energy) < 1e-6; });
        updatingWidgets_ = true;
        manualHistogramCombo_->setCurrentIndex(descriptorIndex);
        if (line != energyLines_.end()) manualSourceCombo_->setCurrentText(Text(line->source));
        updatingWidgets_ = false;
        RefreshEnergyList(manualEnergyList_, manualSourceCombo_, manualEnergyIndices_);
        for (int energyRow = 0; energyRow < static_cast<int>(manualEnergyIndices_.size()); ++energyRow) {
            if (std::abs(energyLines_[manualEnergyIndices_[energyRow]].energy - point.energy) < 1e-6) {
                manualEnergyList_->setCurrentRow(energyRow);
                break;
            }
        }
        ShowCrystalSpectrum(crystal, *descriptor);
        SetStatus("Selected fitted peak " + FormatNumber(point.energy, 3) +
                  " keV for replacement. Select a new fit range, or choose another energy to add.");
    });
    connect(combinedHistogramCombo_, &QComboBox::currentIndexChanged, this, [this] {
        if (!updatingWidgets_) RefreshCombinedQualityList();
    });
    connect(alignmentHistogramCombo_, &QComboBox::currentIndexChanged, this, [this] {
        if (!updatingWidgets_ && alignmentHistogramCombo_->currentIndex() >= 0) {
            ShowSpectrumAlignment();
        }
    });
    connect(resultSpectrumCombo_, &QComboBox::currentIndexChanged, this, [this] {
        if (!updatingWidgets_ && resultSpectrumCombo_->currentIndex() >= 0) {
            ShowSelectedResultSpectrum();
        }
    });
    connect(resultList_, &QListWidget::currentRowChanged, this, [this] {
        if (updatingWidgets_) return;
        UpdateManualCorrectionForSelection();
        RefreshAlignmentParameterList();
        RefreshResultSpectrumChoices();
        ShowAllResultSpectra();
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
        LoadSavedCustomLines();
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
        list->addItem(Text(FormatNumber(line.energy, 3) + " keV  —  " + line.label +
                           (line.saved ? "  [saved]" : "")));
    }
    if (list->count() > 0) list->setCurrentRow(0);
}

void MainWindow::LoadSavedCustomLines() {
    QSettings settings;
    const int count = settings.beginReadArray("customPeaks");
    for (int index = 0; index < count; ++index) {
        settings.setArrayIndex(index);
        const double energy = settings.value("energy").toDouble();
        const std::string label = settings.value("label", "user supplied").toString()
                                      .trimmed().toStdString();
        if (!(energy > 0.0) || !std::isfinite(energy)) continue;
        const bool duplicate = std::any_of(
            energyLines_.begin(), energyLines_.end(), [&](const EnergyLine& line) {
                return line.source == "Custom" && std::abs(line.energy - energy) < 1e-6;
            });
        if (!duplicate) {
            energyLines_.push_back({energy, label.empty() ? "user supplied" : label,
                                    "Custom", true});
        }
    }
    settings.endArray();
}

void MainWindow::SaveCustomLines() const {
    QSettings settings;
    int savedCount = 0;
    for (const auto& line : energyLines_) {
        if (line.source == "Custom" && line.saved) ++savedCount;
    }
    settings.beginWriteArray("customPeaks", savedCount);
    int savedIndex = 0;
    for (const auto& line : energyLines_) {
        if (line.source != "Custom" || !line.saved) continue;
        settings.setArrayIndex(savedIndex++);
        settings.setValue("energy", line.energy);
        settings.setValue("label", Text(line.label));
    }
    settings.endArray();
    settings.sync();
}

void MainWindow::AddCustomLine(bool saveForFuture) {
    const double energy = customEnergyEntry_->value();
    if (!(energy > 0.0) || !std::isfinite(energy)) {
        SetStatus("Enter a positive custom peak energy first.");
        return;
    }
    const QString enteredLabel = customLabelEntry_->text().trimmed();
    const std::string label = enteredLabel.isEmpty()
        ? "user supplied" : enteredLabel.toStdString();
    auto existing = std::find_if(
        energyLines_.begin(), energyLines_.end(), [&](const EnergyLine& line) {
            return line.source == "Custom" && std::abs(line.energy - energy) < 1e-6;
        });
    if (existing == energyLines_.end()) {
        energyLines_.push_back({energy, label, "Custom", saveForFuture});
    } else if (saveForFuture || !existing->saved) {
        existing->label = label;
        existing->saved = existing->saved || saveForFuture;
    }
    if (saveForFuture) SaveCustomLines();
    PopulateEnergyLines();
    referenceSourceCombo_->setCurrentText("Custom");
    manualSourceCombo_->setCurrentText("Custom");
    RefreshEnergyList(energyList_, referenceSourceCombo_, referenceEnergyIndices_);
    RefreshEnergyList(manualEnergyList_, manualSourceCombo_, manualEnergyIndices_);
    for (int row = 0; row < static_cast<int>(referenceEnergyIndices_.size()); ++row) {
        if (std::abs(energyLines_[referenceEnergyIndices_[row]].energy - energy) < 1e-6) {
            energyList_->setCurrentRow(row);
            manualEnergyList_->setCurrentRow(row);
            break;
        }
    }
    SetStatus(std::string(saveForFuture ? "Saved" : "Added") +
              " custom energy line " + FormatNumber(energy, 3) + " keV" +
              (saveForFuture ? " for future sessions." : " for this session."));
}

void MainWindow::RemoveSelectedCustomLine() {
    if (referenceSourceCombo_->currentText() != "Custom") {
        SetStatus("Select the Custom source and a custom line to remove.");
        return;
    }
    const int row = energyList_->currentRow();
    if (row < 0 || row >= static_cast<int>(referenceEnergyIndices_.size())) {
        SetStatus("Select a custom line to remove.");
        return;
    }
    const std::size_t energyIndex = referenceEnergyIndices_[row];
    if (energyIndex >= energyLines_.size() || energyLines_[energyIndex].source != "Custom") return;
    const double energy = energyLines_[energyIndex].energy;
    const bool wasSaved = energyLines_[energyIndex].saved;
    energyLines_.erase(energyLines_.begin() + static_cast<std::ptrdiff_t>(energyIndex));
    if (wasSaved) SaveCustomLines();
    PopulateEnergyLines();
    referenceSourceCombo_->setCurrentText("Custom");
    manualSourceCombo_->setCurrentText("Custom");
    RefreshEnergyList(energyList_, referenceSourceCombo_, referenceEnergyIndices_);
    RefreshEnergyList(manualEnergyList_, manualSourceCombo_, manualEnergyIndices_);
    SetStatus("Removed custom energy line " + FormatNumber(energy, 3) + " keV" +
              (wasSaved ? " from future sessions." : " from this session."));
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

void MainWindow::SaveProjectDialog() {
    QString path = QFileDialog::getSaveFileName(
        this, "Save HPGe calibration project", {},
        "HPGe calibration projects (*.hpgecal.json);;JSON files (*.json)");
    if (path.isEmpty()) return;
    if (!path.endsWith(".json", Qt::CaseInsensitive)) path += ".hpgecal.json";
    std::string error;
    if (!SaveProject(path.toStdString(), error)) {
        SetStatus(error);
        return;
    }
    SetStatus("Saved the complete calibration project to " + path.toStdString() + ".");
}

void MainWindow::OpenProjectDialog() {
    const QString path = QFileDialog::getOpenFileName(
        this, "Open HPGe calibration project", {},
        "HPGe calibration projects (*.hpgecal.json *.json);;All files (*)");
    if (path.isEmpty()) return;
    std::string error;
    if (!OpenProject(path.toStdString(), error)) {
        SetStatus(error);
        return;
    }
    SetStatus("Restored the complete calibration project from " + path.toStdString() + ".");
}

bool MainWindow::SaveProject(const std::string& path, std::string& error) const {
    error.clear();
    if (path.empty()) {
        error = "A project output path is required.";
        return false;
    }
    if (descriptors_.empty()) {
        error = "Add at least one ROOT file before saving a calibration project.";
        return false;
    }
    QJsonObject root;
    root["format"] = "HPGeCalibratorProject";
    root["version"] = 1;
    const QDir projectDirectory = QFileInfo(Text(path)).absoluteDir();
    std::map<std::string, int> fileIndices;
    QJsonArray files;
    for (const auto& descriptor : descriptors_) {
        if (fileIndices.find(descriptor.filePath) != fileIndices.end()) continue;
        const int index = fileIndices.size();
        fileIndices[descriptor.filePath] = index;
        const QString absolutePath = QFileInfo(Text(descriptor.filePath)).absoluteFilePath();
        files.append(QJsonObject{{"path", absolutePath},
                                 {"relativePath", projectDirectory.relativeFilePath(absolutePath)}});
    }
    root["rootFiles"] = files;

    const auto selectedDescriptorIndices = SelectedDescriptorIndices();
    const std::set<int> selectedDescriptors(selectedDescriptorIndices.begin(),
                                             selectedDescriptorIndices.end());
    QJsonArray datasets;
    for (std::size_t index = 0; index < descriptors_.size(); ++index) {
        const auto& descriptor = descriptors_[index];
        datasets.append(QJsonObject{
            {"id", Text(descriptor.id)},
            {"rootFile", fileIndices.at(descriptor.filePath)},
            {"objectPath", Text(descriptor.objectPath)},
            {"selected", selectedDescriptors.count(static_cast<int>(index)) != 0}});
    }
    root["datasets"] = datasets;

    QJsonArray selectedCrystals;
    for (int crystal : SelectedCrystals()) selectedCrystals.append(crystal);
    const auto descriptorId = [this](const QComboBox* combo) {
        const auto* descriptor = DescriptorForCombo(combo);
        return descriptor ? Text(descriptor->id) : QString();
    };
    root["ui"] = QJsonObject{
        {"orientation", orientationCombo_->currentIndex()},
        {"referenceCrystal", ReferenceCrystal()},
        {"selectedCrystals", selectedCrystals},
        {"sigmaBins", sigmaEntry_->value()},
        {"threshold", thresholdEntry_->value()},
        {"residualLimit", residualLimitEntry_->value()},
        {"alignmentSensitivityPercent", alignmentSensitivityEntry_->value()},
        {"autoTuneAlignment", autoTuneAlignmentEntry_->isChecked()},
        {"alignmentModel", alignmentModelCombo_->currentIndex()},
        {"alignedFitHalfRange", alignedFitHalfRangeEntry_->value()},
        {"referenceDataset", descriptorId(referenceHistogramCombo_)},
        {"manualDataset", descriptorId(manualHistogramCombo_)},
        {"alignmentDataset", descriptorId(alignmentHistogramCombo_)},
        {"residualDataset", residualDatasetCombo_->currentData().toString()},
        {"alignmentCrystal", alignmentCrystalEntry_->value()},
        {"selectedResultCrystal", CurrentResultCrystal()}};

    QJsonArray customLines;
    for (const auto& line : energyLines_) {
        if (line.source != "Custom") continue;
        customLines.append(QJsonObject{{"energy", line.energy},
                                       {"label", Text(line.label)}});
    }
    root["customLines"] = customLines;

    QJsonArray referencePeaks;
    for (const auto& peak : referencePeaks_) {
        referencePeaks.append(QJsonObject{
            {"datasetId", Text(peak.datasetId)}, {"charge", peak.charge},
            {"energy", peak.energy}, {"label", Text(peak.label)},
            {"peakFit", PeakFitJson(peak.peakFit)}});
    }
    root["referencePeaks"] = referencePeaks;

    QJsonArray manualPeaks;
    for (const auto& peak : manualPeaks_) {
        manualPeaks.append(QJsonObject{
            {"datasetId", Text(peak.datasetId)}, {"crystal", peak.crystal},
            {"charge", peak.charge}, {"energy", peak.energy},
            {"label", Text(peak.label)}, {"peakFit", PeakFitJson(peak.peakFit)}});
    }
    root["pendingManualPeaks"] = manualPeaks;

    QJsonArray results;
    for (const auto& [crystal, result] : results_) {
        (void)crystal;
        results.append(CalibrationResultJson(result));
    }
    root["calibrationResults"] = results;

    QJsonArray alignments;
    for (const auto& [key, alignment] : alignmentResults_) {
        alignments.append(QJsonObject{{"crystal", key.first},
                                      {"datasetId", Text(key.second)},
                                      {"alignment", PeakMatchJson(alignment)}});
    }
    root["alignments"] = alignments;

    QSaveFile output(Text(path));
    if (!output.open(QIODevice::WriteOnly)) {
        error = "Could not create project file: " + path;
        return false;
    }
    const QByteArray projectData = QJsonDocument(root).toJson(QJsonDocument::Indented);
    if (output.write(projectData) != projectData.size() || !output.commit()) {
        error = "Could not finish writing project file: " + path;
        return false;
    }
    return true;
}

bool MainWindow::OpenProject(const std::string& path, std::string& error) {
    error.clear();
    QFile input(Text(path));
    if (!input.open(QIODevice::ReadOnly)) {
        error = "Could not open project file: " + path;
        return false;
    }
    QJsonParseError parseError;
    const QJsonDocument document = QJsonDocument::fromJson(input.readAll(), &parseError);
    if (parseError.error != QJsonParseError::NoError || !document.isObject()) {
        error = "Invalid project JSON: " + parseError.errorString().toStdString();
        return false;
    }
    const QJsonObject root = document.object();
    if (root["format"].toString() != "HPGeCalibratorProject" ||
        root["version"].toInt() != 1) {
        error = "Unsupported HPGe calibration project format or version.";
        return false;
    }

    const QDir projectDirectory = QFileInfo(Text(path)).absoluteDir();
    std::vector<std::string> resolvedFiles;
    for (const auto value : root["rootFiles"].toArray()) {
        const auto object = value.toObject();
        QString resolved = object["path"].toString();
        if (!QFileInfo::exists(resolved)) {
            resolved = projectDirectory.absoluteFilePath(object["relativePath"].toString());
        }
        if (!QFileInfo::exists(resolved)) {
            error = "A project ROOT file is missing: " + object["path"].toString().toStdString();
            return false;
        }
        resolvedFiles.push_back(QFileInfo(resolved).absoluteFilePath().toStdString());
    }
    if (resolvedFiles.empty()) {
        error = "The project does not contain any ROOT files.";
        return false;
    }
    const QJsonArray datasetObjects = root["datasets"].toArray();
    RootDataRepository validator;
    std::vector<std::vector<HistogramDescriptor>> discoveredFiles;
    for (const auto& file : resolvedFiles) {
        std::string validationError;
        auto discovered = validator.Discover(file, validationError);
        if (discovered.empty()) {
            error = validationError.empty() ? "No TH2 histograms found in " + file
                                            : validationError;
            return false;
        }
        discoveredFiles.push_back(std::move(discovered));
    }
    for (const auto value : datasetObjects) {
        const auto object = value.toObject();
        const int rootIndex = object["rootFile"].toInt(-1);
        const std::string objectPath = object["objectPath"].toString().toStdString();
        if (rootIndex < 0 || rootIndex >= static_cast<int>(discoveredFiles.size()) ||
            std::none_of(discoveredFiles[static_cast<std::size_t>(rootIndex)].begin(),
                         discoveredFiles[static_cast<std::size_t>(rootIndex)].end(),
                         [&](const HistogramDescriptor& item) {
                             return item.objectPath == objectPath;
                         })) {
            error = "A project TH2 histogram no longer exists: " + objectPath;
            return false;
        }
    }

    repository_.ClearCache();
    descriptors_.clear();
    referencePeaks_.clear();
    manualPeaks_.clear();
    results_.clear();
    alignmentResults_.clear();
    combinedAnalyses_.clear();
    displayedSpectrum_.reset();
    displayedDatasetId_.clear();
    displayedCrystal_ = -1;
    pendingRangeStart_.reset();
    RefreshDatasetWidgets();
    if (!OpenRootFiles(resolvedFiles)) {
        error = "The project ROOT files could not be rediscovered.";
        return false;
    }

    std::map<std::string, std::string> datasetMap;
    for (const auto value : datasetObjects) {
        const auto object = value.toObject();
        const int rootIndex = object["rootFile"].toInt(-1);
        if (rootIndex < 0 || rootIndex >= static_cast<int>(resolvedFiles.size())) continue;
        const std::string objectPath = object["objectPath"].toString().toStdString();
        const auto descriptor = std::find_if(
            descriptors_.begin(), descriptors_.end(), [&](const HistogramDescriptor& item) {
                return item.filePath == resolvedFiles[static_cast<std::size_t>(rootIndex)] &&
                       item.objectPath == objectPath;
            });
        if (descriptor != descriptors_.end()) {
            datasetMap[object["id"].toString().toStdString()] = descriptor->id;
        }
    }
    if (datasetMap.size() != static_cast<std::size_t>(datasetObjects.size())) {
        error = "One or more project TH2 histograms no longer exist in the ROOT files.";
        return false;
    }
    const auto remapDataset = [&](const std::string& oldId) -> std::string {
        const auto mapped = datasetMap.find(oldId);
        return mapped == datasetMap.end() ? std::string() : mapped->second;
    };

    energyLines_.clear();
    PopulateEnergyLines();
    for (const auto value : root["customLines"].toArray()) {
        const auto object = value.toObject();
        const double energy = object["energy"].toDouble();
        const std::string label = object["label"].toString("user supplied").toStdString();
        if (!(energy > 0.0) || !std::isfinite(energy)) continue;
        const auto existing = std::find_if(
            energyLines_.begin(), energyLines_.end(), [&](const EnergyLine& line) {
                return line.source == "Custom" && std::abs(line.energy - energy) < 1e-6;
            });
        if (existing == energyLines_.end()) {
            energyLines_.push_back({energy, label, "Custom", false});
        }
    }
    PopulateEnergyLines();

    const QJsonObject ui = root["ui"].toObject();
    updatingWidgets_ = true;
    orientationCombo_->setCurrentIndex(std::clamp(ui["orientation"].toInt(), 0, 1));
    referenceCrystalEntry_->setValue(ui["referenceCrystal"].toInt());
    sigmaEntry_->setValue(ui["sigmaBins"].toDouble(2.0));
    thresholdEntry_->setValue(ui["threshold"].toDouble(0.05));
    residualLimitEntry_->setValue(ui["residualLimit"].toDouble(1.0));
    alignmentSensitivityEntry_->setValue(ui["alignmentSensitivityPercent"].toDouble(35.0));
    autoTuneAlignmentEntry_->setChecked(ui["autoTuneAlignment"].toBool(true));
    alignmentModelCombo_->setCurrentIndex(std::clamp(ui["alignmentModel"].toInt(), 0, 2));
    alignedFitHalfRangeEntry_->setValue(ui["alignedFitHalfRange"].toDouble(35.0));
    alignmentCrystalEntry_->setValue(ui["alignmentCrystal"].toInt(1));

    histogramList_->clearSelection();
    for (const auto value : datasetObjects) {
        const auto object = value.toObject();
        if (!object["selected"].toBool()) continue;
        const std::string mappedId = remapDataset(object["id"].toString().toStdString());
        for (std::size_t index = 0; index < descriptors_.size(); ++index) {
            if (descriptors_[index].id == mappedId) {
                histogramList_->item(static_cast<int>(index))->setSelected(true);
                break;
            }
        }
    }
    crystalList_->clearSelection();
    for (const auto value : ui["selectedCrystals"].toArray()) {
        const int crystal = value.toInt(-1);
        if (crystal >= 0 && crystal < crystalList_->count()) {
            crystalList_->item(crystal)->setSelected(true);
        }
    }
    const auto setDescriptorCombo = [&](QComboBox* combo, const QJsonValue& oldIdValue) {
        const std::string mappedId = remapDataset(oldIdValue.toString().toStdString());
        for (std::size_t index = 0; index < descriptors_.size(); ++index) {
            if (descriptors_[index].id == mappedId) {
                combo->setCurrentIndex(static_cast<int>(index));
                return;
            }
        }
    };
    setDescriptorCombo(referenceHistogramCombo_, ui["referenceDataset"]);
    setDescriptorCombo(manualHistogramCombo_, ui["manualDataset"]);
    setDescriptorCombo(alignmentHistogramCombo_, ui["alignmentDataset"]);
    const std::string residualDataset = remapDataset(
        ui["residualDataset"].toString().toStdString());
    for (int index = 0; index < residualDatasetCombo_->count(); ++index) {
        if (residualDatasetCombo_->itemData(index).toString().toStdString() == residualDataset) {
            residualDatasetCombo_->setCurrentIndex(index);
            break;
        }
    }
    updatingWidgets_ = false;

    for (const auto value : root["referencePeaks"].toArray()) {
        const auto object = value.toObject();
        const std::string datasetId = remapDataset(
            object["datasetId"].toString().toStdString());
        if (datasetId.empty()) continue;
        referencePeaks_.push_back({datasetId, object["charge"].toDouble(),
            object["energy"].toDouble(), object["label"].toString().toStdString(),
            ReadPeakFit(object["peakFit"])});
    }
    for (const auto value : root["pendingManualPeaks"].toArray()) {
        const auto object = value.toObject();
        const std::string datasetId = remapDataset(
            object["datasetId"].toString().toStdString());
        if (datasetId.empty()) continue;
        manualPeaks_.push_back({datasetId, object["crystal"].toInt(-1),
            object["charge"].toDouble(), object["energy"].toDouble(),
            object["label"].toString().toStdString(), ReadPeakFit(object["peakFit"])});
    }
    for (const auto value : root["calibrationResults"].toArray()) {
        CalibrationResult result = ReadCalibrationResult(value);
        for (auto& point : result.points) point.datasetId = remapDataset(point.datasetId);
        result.points.erase(std::remove_if(result.points.begin(), result.points.end(),
            [](const CalibrationPoint& point) { return point.datasetId.empty(); }),
            result.points.end());
        if (result.crystal >= 0 && result.crystal < 64) {
            results_[result.crystal] = std::move(result);
        }
    }
    for (const auto value : root["alignments"].toArray()) {
        const auto object = value.toObject();
        const int crystal = object["crystal"].toInt(-1);
        const std::string datasetId = remapDataset(
            object["datasetId"].toString().toStdString());
        if (crystal >= 0 && crystal < 64 && !datasetId.empty()) {
            alignmentResults_[{crystal, datasetId}] = ReadPeakMatch(object["alignment"]);
        }
    }

    RefreshReferencePeakList();
    RefreshResults();
    const int selectedResult = ui["selectedResultCrystal"].toInt(-1);
    for (int row = 0; row < resultList_->count(); ++row) {
        if (resultList_->item(row)->data(Qt::UserRole).toInt() == selectedResult) {
            resultList_->setCurrentRow(row);
            break;
        }
    }
    RefreshAlignmentParameterList();
    RefreshFittedPointList();
    RefreshManualPeakList();
    EvaluateCombinedSpectra();
    ShowReferenceSpectrum();
    return true;
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
    residualDatasetCombo_->clear();
    residualDatasetCombo_->addItem("All source histograms", QString());
    for (std::size_t i = 0; i < descriptors_.size(); ++i) {
        const auto& descriptor = descriptors_[i];
        histogramList_->addItem(Text(descriptor.displayName + "  [" +
                                     std::to_string(descriptor.xBins) + " x " +
                                     std::to_string(descriptor.yBins) + "]"));
        const QString name = Text(descriptor.displayName);
        referenceHistogramCombo_->addItem(name, static_cast<int>(i));
        manualHistogramCombo_->addItem(name, static_cast<int>(i));
        alignmentHistogramCombo_->addItem(name, static_cast<int>(i));
        residualDatasetCombo_->addItem(name, Text(descriptor.id));
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

double MainWindow::AlignmentSensitivity() const {
    return alignmentSensitivityEntry_->value() / 100.0;
}

CalibrationEngine::AlignmentModel MainWindow::AlignmentModel() const {
    if (alignmentModelCombo_->currentIndex() == 1) {
        return CalibrationEngine::AlignmentModel::Affine;
    }
    if (alignmentModelCombo_->currentIndex() == 2) {
        return CalibrationEngine::AlignmentModel::Quadratic;
    }
    return CalibrationEngine::AlignmentModel::Auto;
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

void MainWindow::RefreshResultSpectrumChoices() {
    const int crystal = CurrentResultCrystal();
    const auto result = results_.find(crystal);
    const int previousDescriptor = resultSpectrumCombo_->currentData().toInt();
    const bool wasUpdating = updatingWidgets_;
    updatingWidgets_ = true;
    resultSpectrumCombo_->clear();
    if (crystal >= 0 && result != results_.end()) {
        for (std::size_t index = 0; index < descriptors_.size(); ++index) {
            const auto& descriptor = descriptors_[index];
            const bool used = std::any_of(
                result->second.points.begin(), result->second.points.end(),
                [&](const CalibrationPoint& point) {
                    return point.datasetId == descriptor.id && point.peakFit.success;
                });
            if (used) {
                resultSpectrumCombo_->addItem(Text(descriptor.displayName),
                                              static_cast<int>(index));
            }
        }
    }
    const int previousRow = resultSpectrumCombo_->findData(previousDescriptor);
    if (previousRow >= 0) resultSpectrumCombo_->setCurrentIndex(previousRow);
    else if (resultSpectrumCombo_->count() > 0) resultSpectrumCombo_->setCurrentIndex(0);
    updatingWidgets_ = wasUpdating;
}

void MainWindow::ShowSelectedResultSpectrum() {
    const int crystal = CurrentResultCrystal();
    const auto result = results_.find(crystal);
    if (crystal < 0 || result == results_.end()) {
        SetStatus("Select one calibration result first.");
        return;
    }

    if (resultSpectrumCombo_->count() == 0) RefreshResultSpectrumChoices();
    const HistogramDescriptor* descriptor = DescriptorForCombo(resultSpectrumCombo_);
    const auto fitCountFor = [&](const std::string& datasetId) {
        return static_cast<int>(std::count_if(result->second.points.begin(), result->second.points.end(),
            [&](const CalibrationPoint& point) {
                return point.datasetId == datasetId && point.peakFit.success;
            }));
    };
    if (!descriptor) {
        SetStatus("No fitted source spectrum is available for the selected result.");
        return;
    }
    const int descriptorIndex = static_cast<int>(descriptor - descriptors_.data());
    const bool wasUpdating = updatingWidgets_;
    updatingWidgets_ = true;
    manualHistogramCombo_->setCurrentIndex(descriptorIndex);
    updatingWidgets_ = wasUpdating;
    ShowCrystalSpectrum(crystal, *descriptor);
    const int fittedPeaks = fitCountFor(descriptor->id);
    SetStatus("Showing " + descriptor->displayName + " for Crystal " +
              std::to_string(crystal) + " with " + std::to_string(fittedPeaks) +
              " stored fitted peak(s) in red (source " +
              std::to_string(resultSpectrumCombo_->currentIndex() + 1) + "/" +
              std::to_string(resultSpectrumCombo_->count()) + ").");
}

void MainWindow::ShowAllResultSpectra() {
    const int crystal = CurrentResultCrystal();
    const auto result = results_.find(crystal);
    if (crystal < 0 || result == results_.end()) {
        SetStatus("Select one calibration result first.");
        return;
    }
    if (resultSpectrumCombo_->count() == 0) RefreshResultSpectrumChoices();
    if (resultSpectrumCombo_->count() == 0) {
        SetStatus("The selected result has no stored fitted source spectra.");
        return;
    }

    const std::array<QColor, 8> colors{
        QColor("#2563eb"), QColor("#16a34a"), QColor("#7c3aed"), QColor("#d97706"),
        QColor("#0891b2"), QColor("#db2777"), QColor("#4f46e5"), QColor("#65a30d")};
    std::vector<PlotSeries> series;
    std::vector<PlotMarker> markers;
    int sourceCount = 0;
    int fitCount = 0;
    for (int row = 0; row < resultSpectrumCombo_->count(); ++row) {
        const int descriptorIndex = resultSpectrumCombo_->itemData(row).toInt();
        if (descriptorIndex < 0 || descriptorIndex >= static_cast<int>(descriptors_.size())) continue;
        const auto& descriptor = descriptors_[descriptorIndex];
        std::string error;
        auto spectrum = repository_.ProjectCrystal(descriptor, crystal, Orientation(), error);
        if (!spectrum) continue;
        const double maximum = std::max(spectrum->GetMaximum(), 1.0);
        const double offset = 1.15 * sourceCount;
        PlotSeries spectrumSeries = HistogramSeries(
            *spectrum, colors[static_cast<std::size_t>(sourceCount) % colors.size()],
            descriptor.displayName + "  (+" + FormatNumber(offset, 2) + ")");
        for (double& value : spectrumSeries.y) value = value / maximum + offset;
        series.push_back(std::move(spectrumSeries));
        for (const auto& point : result->second.points) {
            if (point.datasetId != descriptor.id || !point.peakFit.success) continue;
            PlotSeries fitted = FitSeries(
                point.peakFit, FormatNumber(point.energy, 1) + " keV (" +
                                   descriptor.objectPath + ")", point.manual);
            for (double& value : fitted.y) value = value / maximum + offset;
            series.push_back(std::move(fitted));
            markers.push_back({point.peakFit.centroid,
                               FormatNumber(point.energy, 1) + " keV",
                               QColor("#dc2626"), point.manual});
            ++fitCount;
        }
        ++sourceCount;
    }
    if (sourceCount == 0) {
        primaryPlot_->Clear("No fitted source spectra could be projected for this result.");
        SetStatus("No fitted source spectra could be projected for Crystal " +
                  std::to_string(crystal) + ".");
        return;
    }
    displayedSpectrum_.reset();
    displayedDatasetId_.clear();
    displayedCrystal_ = -1;
    pendingRangeStart_.reset();
    SetSecondaryPlotVisible(false);
    primaryPlot_->SetPlot("All fitted source spectra — Crystal " + std::to_string(crystal),
                          "Charge", "Normalized counts + source offset",
                          std::move(series), std::move(markers));
    SetStatus("Showing all " + std::to_string(sourceCount) +
              " fitted source spectra used by Crystal " + std::to_string(crystal) +
              " with " + std::to_string(fitCount) + " fitted peak curve(s) in red.");
}

void MainWindow::UpdateManualCorrectionForSelection() {
    const int crystal = CurrentResultCrystal();
    const bool selected = crystal >= 0 && results_.count(crystal) != 0;
    manualCorrectionGroup_->setEnabled(selected);
    manualCorrectionGroup_->setTitle(selected
        ? QString("Manual correction — Crystal %1 only").arg(crystal, 2, 10, QChar('0'))
        : QString("Manual correction — select one calibration result"));
    pendingRangeStart_.reset();
    RefreshFittedPointList();
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
    const auto descriptor = std::find_if(
        descriptors_.begin(), descriptors_.end(), [&](const HistogramDescriptor& item) {
            return item.id == displayedDatasetId_;
        });
    const std::string datasetName = descriptor == descriptors_.end()
        ? displayedDatasetId_ : descriptor->displayName;
    primaryPlot_->SetPlot("Crystal " + std::to_string(displayedCrystal_) +
                          " spectrum — " + datasetName,
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
    const auto result = results_.find(crystal);
    const bool replaces = result != results_.end() && std::any_of(
        result->second.points.begin(), result->second.points.end(),
        [&](const CalibrationPoint& point) {
            return point.datasetId == descriptor->id &&
                   std::abs(point.energy - energy->energy) < 1e-6;
        });
    SetStatus(std::string("Manual RadWare centroid fitted as a pending ") +
              (replaces ? "replacement" : "new peak") +
              ". Apply pending peaks to recompute only the calibration polynomial.");
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

void MainWindow::RefreshAlignmentParameterList() {
    alignmentParameterList_->clear();
    const int crystal = CurrentResultCrystal();
    if (crystal < 0) return;
    for (const auto& descriptor : descriptors_) {
        const auto stored = alignmentResults_.find({crystal, descriptor.id});
        if (stored == alignmentResults_.end()) continue;
        const auto& alignment = stored->second;
        const int matched = static_cast<int>(
            std::count(alignment.matched.begin(), alignment.matched.end(), true));
        const std::string state = alignment.success ? "OK" : "FAILED";
        const std::string model = alignment.quadraticModel ? "quadratic" : "affine";
        alignmentParameterList_->addItem(Text(
            "[" + state + "] " + descriptor.displayName + " | C" +
            std::to_string(crystal) + " | " + model + " | a0=" +
            FormatPrecise(alignment.offset) + " | a1=" +
            FormatPrecise(alignment.scale) + " | a2=" +
            FormatPrecise(alignment.quadratic) + " | cost=" +
            FormatNumber(alignment.alignmentCost, 4) + " | matched=" +
            std::to_string(matched) + "/" +
            std::to_string(alignment.referenceCharges.size()) + " | sensitivity=" +
            FormatNumber(100.0 * alignment.referenceSensitivity, 1) + "%/" +
            FormatNumber(100.0 * alignment.targetSensitivity, 1) + "%"));
    }
    if (alignmentParameterList_->count() > 0) alignmentParameterList_->setCurrentRow(0);
}

void MainWindow::RefreshFittedPointList() {
    fittedPointList_->clear();
    const int crystal = CurrentResultCrystal();
    const auto result = results_.find(crystal);
    if (result == results_.end()) return;
    for (std::size_t index = 0; index < result->second.points.size(); ++index) {
        const auto& point = result->second.points[index];
        const auto descriptor = std::find_if(
            descriptors_.begin(), descriptors_.end(), [&](const HistogramDescriptor& value) {
                return value.id == point.datasetId;
            });
        const std::string dataset = descriptor == descriptors_.end()
            ? point.datasetId : descriptor->displayName;
        const std::string fitQuality = point.peakFit.ndf > 0
            ? " chi2/ndf=" + FormatNumber(point.peakFit.chi2 /
                                           static_cast<double>(point.peakFit.ndf), 2)
            : "";
        auto* item = new QListWidgetItem(Text(
            (point.manual ? "[manual] " : "[auto] ") + dataset + " | " +
            FormatNumber(point.energy, 3) + " keV | centroid " +
            FormatNumber(point.charge, 3) + " +/- " +
            FormatNumber(point.chargeError, 3) + fitQuality), fittedPointList_);
        item->setData(Qt::UserRole, static_cast<int>(index));
        item->setForeground(point.manual ? QColor("#b45309") : QColor("#15803d"));
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
    options.alignmentSensitivity = AlignmentSensitivity();
    options.autoTuneAlignmentSensitivity = autoTuneAlignmentEntry_->isChecked();
    options.alignmentModel = AlignmentModel();
    for (int descriptorIndex : SelectedDescriptorIndices()) {
        const auto& descriptor = descriptors_[descriptorIndex];
        std::vector<ReferencePeak> references;
        for (const auto& peak : referencePeaks_) if (peak.datasetId == descriptor.id) references.push_back(peak);
        if (references.empty()) continue;
        std::sort(references.begin(), references.end(), [](const ReferencePeak& a, const ReferencePeak& b) {
            return a.charge < b.charge;
        });
        if (crystal == ReferenceCrystal()) {
            PeakMatchResult identity;
            identity.success = true;
            identity.referenceCharges.reserve(references.size());
            identity.charges.reserve(references.size());
            identity.matched.assign(references.size(), true);
            for (const auto& ref : references) {
                identity.referenceCharges.push_back(ref.peakFit.centroid);
                identity.charges.push_back(ref.peakFit.centroid);
                points.push_back({descriptor.id, ref.peakFit.centroid, ref.energy,
                                  ref.peakFit.centroidError, false, 0.0, ref.peakFit});
            }
            alignmentResults_[{crystal, descriptor.id}] = std::move(identity);
            continue;
        }
        std::string error;
        auto referenceSpectrum = repository_.ProjectCrystal(descriptor, ReferenceCrystal(), Orientation(), error);
        auto spectrum = repository_.ProjectCrystal(descriptor, crystal, Orientation(), error);
        if (!referenceSpectrum || !spectrum) continue;
        const auto alignment = CalibrationEngine::AlignSpectrumPatterns(
            *referenceSpectrum, *spectrum, options);
        alignmentResults_[{crystal, descriptor.id}] = alignment;
        if (!alignment.success) continue;
        const double halfRange = alignedFitHalfRangeEntry_->value();
        for (const auto& reference : references) {
            // The interval is defined exactly in the aligned reference-charge
            // coordinate and mapped at both boundaries. No second peak finder
            // is allowed to move it to an unrelated target candidate.
            const double low = CalibrationEngine::MapReferenceCharge(
                alignment, reference.peakFit.centroid - halfRange);
            const double high = CalibrationEngine::MapReferenceCharge(
                alignment, reference.peakFit.centroid + halfRange);
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
    alignmentResults_.clear();
    for (int crystal : crystals) {
        results_[crystal] = CalibrateCrystal(crystal);
        QApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    RefreshResults();
    EvaluateCombinedSpectra();
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

void MainWindow::EvaluateCombinedSpectra() {
    combinedAnalyses_.clear();
    updatingWidgets_ = true;
    combinedHistogramCombo_->clear();
    combinedQualityList_->clear();
    updatingWidgets_ = false;

    for (int descriptorIndex : SelectedDescriptorIndices()) {
        const auto& descriptor = descriptors_[descriptorIndex];
        std::vector<std::shared_ptr<TH1D>> ownedSpectra;
        std::vector<CalibratedSpectrumInput> inputs;
        for (const auto& [crystal, result] : results_) {
            if (!result.success) continue;
            std::string error;
            auto spectrum = repository_.ProjectCrystal(descriptor, crystal, Orientation(), error);
            if (!spectrum) continue;
            inputs.push_back({spectrum.get(), result.p0, result.p1, result.p2});
            ownedSpectra.push_back(std::move(spectrum));
        }
        std::string error;
        auto combined = CombinedSpectrumAnalyzer::Combine(descriptor.id, inputs, error);
        if (!combined) continue;

        CombinedDatasetAnalysis analysis;
        analysis.datasetId = descriptor.id;
        analysis.spectrum = std::move(combined);
        analysis.crystalCount = static_cast<int>(inputs.size());
        std::vector<const ReferencePeak*> references;
        for (const auto& reference : referencePeaks_) {
            if (reference.datasetId == descriptor.id) references.push_back(&reference);
        }
        std::sort(references.begin(), references.end(), [](const ReferencePeak* first,
                                                            const ReferencePeak* second) {
            return first->energy < second->energy;
        });
        for (const ReferencePeak* reference : references) {
            double halfWindow = 10.0;
            const auto referenceResult = results_.find(ReferenceCrystal());
            if (referenceResult != results_.end() && referenceResult->second.success) {
                const auto& calibration = referenceResult->second;
                const auto energyAt = [&](double charge) {
                    return calibration.p0 + calibration.p1 * charge +
                           calibration.p2 * charge * charge;
                };
                halfWindow = 1.25 * std::max(
                    std::abs(energyAt(reference->peakFit.rangeLow) - reference->energy),
                    std::abs(energyAt(reference->peakFit.rangeHigh) - reference->energy));
                halfWindow = std::clamp(halfWindow, 4.0, 100.0);
            }
            analysis.peaks.push_back(CombinedSpectrumAnalyzer::EvaluatePeak(
                *analysis.spectrum, descriptor.id, reference->energy, halfWindow));
        }
        combinedAnalyses_[descriptor.id] = std::move(analysis);
        combinedHistogramCombo_->addItem(Text(descriptor.displayName), Text(descriptor.id));
    }
    if (combinedHistogramCombo_->count() > 0) combinedHistogramCombo_->setCurrentIndex(0);
    RefreshCombinedQualityList();
}

void MainWindow::RefreshCombinedQualityList() {
    combinedQualityList_->clear();
    const std::string datasetId = combinedHistogramCombo_->currentData().toString().toStdString();
    const auto found = combinedAnalyses_.find(datasetId);
    if (found == combinedAnalyses_.end()) return;
    for (const auto& peak : found->second.peaks) {
        if (!peak.success) {
            combinedQualityList_->addItem(Text("FAIL | " + FormatNumber(peak.expectedEnergy, 3) +
                " keV | " + peak.status));
            continue;
        }
        combinedQualityList_->addItem(Text(
            FormatNumber(peak.expectedEnergy, 3) + " keV | centroid " +
            FormatNumber(peak.fittedEnergy, 3) + " | residual " +
            FormatNumber(peak.residualKeV, 3) + " keV | FWHM " +
            FormatNumber(peak.fwhmKeV, 3) + " keV | resolution " +
            FormatNumber(peak.resolutionPercent, 3) + "%"));
    }
}

void MainWindow::ShowCombinedSpectrum() {
    const std::string datasetId = combinedHistogramCombo_->currentData().toString().toStdString();
    const auto found = combinedAnalyses_.find(datasetId);
    if (found == combinedAnalyses_.end() || !found->second.spectrum) {
        SetStatus("Run calibration before viewing the combined spectrum.");
        return;
    }
    const auto descriptor = std::find_if(descriptors_.begin(), descriptors_.end(),
        [&](const HistogramDescriptor& item) { return item.id == datasetId; });
    const std::string datasetName = descriptor == descriptors_.end()
        ? datasetId : descriptor->displayName;
    std::vector<PlotSeries> spectrumSeries{
        HistogramSeries(*found->second.spectrum, QColor("#2563eb"), "Combined spectrum")};
    std::vector<PlotMarker> markers;
    PlotSeries residuals;
    residuals.name = "Combined centroid residual";
    residuals.color = QColor("#16a34a");
    residuals.points = true;
    int successful = 0;
    for (const auto& peak : found->second.peaks) {
        if (!peak.success) continue;
        ++successful;
        spectrumSeries.push_back(FitSeries(
            peak.peakFit, FormatNumber(peak.expectedEnergy, 1) + " keV / FWHM " +
                          FormatNumber(peak.fwhmKeV, 2)));
        markers.push_back({peak.fittedEnergy,
            FormatNumber(peak.expectedEnergy, 1) + " keV", QColor("#dc2626"), false});
        residuals.x.push_back(peak.expectedEnergy);
        residuals.y.push_back(peak.residualKeV);
    }
    displayedSpectrum_.reset();
    displayedDatasetId_.clear();
    displayedCrystal_ = -1;
    pendingRangeStart_.reset();
    SetSecondaryPlotVisible(true);
    primaryPlot_->SetPlot("Combined calibrated spectrum — " + datasetName,
                          "Energy (keV)", "Counts", std::move(spectrumSeries),
                          std::move(markers));
    if (!residuals.x.empty()) {
        secondaryPlot_->SetPlot("Combined-spectrum residuals versus energy",
                                "Peak energy (keV)", "Fitted - expected (keV)",
                                {std::move(residuals)}, {}, {false, 0.15});
    } else {
        secondaryPlot_->Clear("No combined peaks were fitted");
    }
    SetStatus("Combined " + std::to_string(found->second.crystalCount) +
              " calibrated crystals and refitted " + std::to_string(successful) + " peak(s).");
}

void MainWindow::ShowSpectrumAlignment() {
    const auto* descriptor = DescriptorForCombo(alignmentHistogramCombo_);
    if (!descriptor) {
        alignmentPreviewParameters_->setText("No alignment parameters are available.");
        SetStatus("Choose a source histogram for the alignment preview.");
        return;
    }
    const std::string previewName = descriptor->displayName;
    displayedSpectrum_.reset();
    displayedDatasetId_.clear();
    displayedCrystal_ = -1;
    pendingRangeStart_.reset();
    SetSecondaryPlotVisible(false);
    primaryPlot_->Clear("Loading alignment preview for " + previewName + "...");
    const int referenceCrystal = ReferenceCrystal();
    const int targetCrystal = alignmentCrystalEntry_->value();
    std::string error;
    auto reference = repository_.ProjectCrystal(*descriptor, referenceCrystal, Orientation(), error);
    auto target = repository_.ProjectCrystal(*descriptor, targetCrystal, Orientation(), error);
    if (!reference || !target) {
        alignmentPreviewParameters_->setText("Alignment failed: " + Text(error));
        primaryPlot_->Clear("Alignment preview unavailable for " + previewName + ": " + error);
        SetStatus(error);
        return;
    }
    CalibrationEngine::SearchOptions options;
    options.sigmaBins = sigmaEntry_->value();
    options.threshold = thresholdEntry_->value();
    options.alignmentSensitivity = AlignmentSensitivity();
    options.autoTuneAlignmentSensitivity = autoTuneAlignmentEntry_->isChecked();
    options.alignmentModel = AlignmentModel();
    const auto match = CalibrationEngine::AlignSpectrumPatterns(*reference, *target, options);
    if (!match.success || !(match.scale > 0.0) || !std::isfinite(match.scale) || !std::isfinite(match.offset)) {
        alignmentPreviewParameters_->setText(
            "Alignment failed: not enough corresponding pattern peaks.");
        primaryPlot_->Clear("Alignment preview unavailable for " + previewName +
                            ": not enough corresponding peaks were found.");
        SetStatus("Could not align " + previewName +
                  ": not enough corresponding peaks were found.");
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
    primaryPlot_->SetPlot("Pre-calibration alignment — " + previewName,
                          "Reference-spectrum charge", "Normalized counts",
                          {std::move(referenceSeries), std::move(targetSeries)}, std::move(markers));
    const int matched = static_cast<int>(std::count(match.matched.begin(), match.matched.end(), true));
    alignmentPreviewParameters_->setText(Text(
        std::string(match.quadraticModel ? "Quadratic" : "Affine") +
        " parameters (select to copy): a0=" + FormatPrecise(match.offset) +
        ", a1=" + FormatPrecise(match.scale) +
        ", a2=" + FormatPrecise(match.quadratic) +
        "; cost=" + FormatNumber(match.alignmentCost, 4) +
        "; matched=" + std::to_string(matched) + "/" +
        std::to_string(match.referenceCharges.size()) +
        "; sensitivity R/T=" +
        FormatNumber(100.0 * match.referenceSensitivity, 1) + "%/" +
        FormatNumber(100.0 * match.targetSensitivity, 1) + "%"));
    SetStatus("Alignment preview for " + previewName + ": " +
              std::to_string(matched) + "/" +
              std::to_string(match.referenceCharges.size()) + " pattern peaks; " +
              (match.quadraticModel ? "2nd-order polynomial" : "affine") +
              " mapping; tuned sensitivity R/T " +
              FormatNumber(100.0 * match.referenceSensitivity, 0) + "%/" +
              FormatNumber(100.0 * match.targetSensitivity, 0) +
              "%; minimized pattern cost " + FormatNumber(match.alignmentCost, 3) +
              "; target charge = " +
              FormatNumber(match.offset, 3) + " + " + FormatNumber(match.scale, 6) + " q + " +
              FormatNumber(match.quadratic, 9) + " q^2. No energy calibration applied.");
}

void MainWindow::RefitSelectedCrystal() {
    const int crystal = CurrentResultCrystal();
    const auto existingResult = results_.find(crystal);
    if (crystal < 0 || existingResult == results_.end()) {
        SetStatus("Select a calibration result to refit.");
        return;
    }
    std::vector<CalibrationPoint> points = existingResult->second.points;
    int replaced = 0;
    int added = 0;
    for (const auto& manual : manualPeaks_) {
        if (manual.crystal != crystal) continue;
        const auto point = std::find_if(points.begin(), points.end(),
            [&](const CalibrationPoint& value) {
                return value.datasetId == manual.datasetId &&
                       std::abs(value.energy - manual.energy) < 1e-6;
            });
        CalibrationPoint replacement{manual.datasetId, manual.peakFit.centroid,
            manual.energy, manual.peakFit.centroidError, true, 0.0, manual.peakFit};
        if (point == points.end()) {
            points.push_back(std::move(replacement));
            ++added;
        } else {
            *point = std::move(replacement);
            ++replaced;
        }
    }
    if (replaced + added == 0) {
        SetStatus("No pending fitted peaks were selected. Existing peak fits were left unchanged.");
        return;
    }
    results_[crystal] = CalibrationEngine::FitSecondOrder(
        crystal, std::move(points), residualLimitEntry_->value());
    manualPeaks_.erase(std::remove_if(manualPeaks_.begin(), manualPeaks_.end(),
        [crystal](const ManualPeak& peak) { return peak.crystal == crystal; }),
        manualPeaks_.end());
    RefreshResults();
    EvaluateCombinedSpectra();
    for (int row = 0; row < resultList_->count(); ++row) {
        if (resultList_->item(row)->data(Qt::UserRole).toInt() == crystal) {
            resultList_->setCurrentRow(row);
            break;
        }
    }
    RefreshFittedPointList();
    RefreshManualPeakList();
    const auto* descriptor = DescriptorForCombo(manualHistogramCombo_);
    if (descriptor) ShowCrystalSpectrum(crystal, *descriptor);
    SetStatus("Crystal " + std::to_string(crystal) + " calibration updated: " +
              std::to_string(replaced) + " peak(s) replaced, " +
              std::to_string(added) +
              " added; every other fitted centroid and fit curve was preserved.");
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
    RefreshAlignmentParameterList();
    RefreshResultSpectrumChoices();
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
        residuals.x.push_back(point.energy);
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
                            "Peak energy (keV)", "Energy - fit (keV)", {std::move(residuals)}, {},
                            {false, 0.15});
    SetStatus("Crystal " + std::to_string(crystal) + ": " + result.status);
}

void MainWindow::ShowResidualsByCrystal() {
    if (results_.empty()) {
        SetStatus("Run or restore a calibration before plotting residuals by crystal.");
        return;
    }
    const std::string selectedDataset =
        residualDatasetCombo_->currentData().toString().toStdString();
    using ResidualKey = std::pair<std::string, double>;
    std::map<ResidualKey, PlotSeries> byLine;
    std::map<int, std::vector<double>> byCrystal;
    constexpr std::array<const char*, 10> colors{
        "#2563eb", "#dc2626", "#16a34a", "#d97706", "#7c3aed",
        "#0891b2", "#be185d", "#4d7c0f", "#9333ea", "#0f766e"};
    int colorIndex = 0;
    for (const auto& [crystal, result] : results_) {
        if (!result.success) continue;
        for (const auto& point : result.points) {
            if (!selectedDataset.empty() && point.datasetId != selectedDataset) continue;
            const ResidualKey key{point.datasetId, point.energy};
            auto [series, inserted] = byLine.try_emplace(key);
            if (inserted) {
                const auto descriptor = std::find_if(
                    descriptors_.begin(), descriptors_.end(), [&](const HistogramDescriptor& item) {
                        return item.id == point.datasetId;
                    });
                const std::string dataset = descriptor == descriptors_.end()
                    ? point.datasetId : descriptor->objectPath;
                series->second.name = dataset + " / " + FormatNumber(point.energy, 3) + " keV";
                series->second.color = QColor(colors[static_cast<std::size_t>(colorIndex) %
                                                     colors.size()]);
                series->second.points = true;
                ++colorIndex;
            }
            series->second.x.push_back(crystal);
            series->second.y.push_back(point.residual);
            byCrystal[crystal].push_back(point.residual);
        }
    }
    if (byLine.empty()) {
        SetStatus("No successful fitted residuals are available for the selected source.");
        return;
    }
    PlotSeries bounds;
    bounds.x = {0.0, 63.0};
    bounds.y = {0.0, 0.0};
    bounds.color = QColor(0, 0, 0, 0);
    std::vector<PlotSeries> residualSeries{std::move(bounds)};
    for (auto& [key, series] : byLine) {
        (void)key;
        residualSeries.push_back(std::move(series));
    }
    PlotSeries rms;
    rms.name = "RMS residual";
    rms.color = QColor("#2563eb");
    rms.points = true;
    PlotSeries maximum;
    maximum.name = "Maximum |residual|";
    maximum.color = QColor("#dc2626");
    maximum.points = true;
    for (const auto& [crystal, residuals] : byCrystal) {
        double sumSquares = 0.0;
        double maximumAbsolute = 0.0;
        for (double residual : residuals) {
            sumSquares += residual * residual;
            maximumAbsolute = std::max(maximumAbsolute, std::abs(residual));
        }
        rms.x.push_back(crystal);
        rms.y.push_back(std::sqrt(sumSquares / static_cast<double>(residuals.size())));
        maximum.x.push_back(crystal);
        maximum.y.push_back(maximumAbsolute);
    }
    PlotSeries qualityBounds;
    qualityBounds.x = {0.0, 63.0};
    qualityBounds.y = {0.0, 0.0};
    qualityBounds.color = QColor(0, 0, 0, 0);
    displayedSpectrum_.reset();
    displayedDatasetId_.clear();
    displayedCrystal_ = -1;
    pendingRangeStart_.reset();
    SetSecondaryPlotVisible(true);
    primaryPlot_->SetPlot("Energy residuals versus detector crystal",
                          "Detector crystal (0-63)", "Energy residual (keV)",
                          std::move(residualSeries));
    secondaryPlot_->SetPlot("Per-crystal residual quality",
                            "Detector crystal (0-63)", "Residual magnitude (keV)",
                            {std::move(qualityBounds), std::move(rms), std::move(maximum)});
    SetStatus("Showing " + std::to_string(byLine.size()) +
              " energy-residual series across " + std::to_string(byCrystal.size()) +
              " calibrated crystal(s) on the 0-63 detector axis.");
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
    WriteCsvRow(output, {"record", "crystal", "status", "needs_review", "p0", "p1", "p2",
        "chi2", "ndf", "residual_rms_keV", "dataset", "charge", "energy_keV",
        "residual_keV", "manual", "charge_error", "range_low", "range_high",
        "peak_sigma", "peak_chi2", "peak_ndf", "tail_fraction", "beta", "step_fraction",
        "fwhm_keV", "resolution_percent", "combined_crystals", "combined_fitted_energy_keV"});
    for (const auto& [crystal, result] : results_) {
        std::vector<std::string> fitRow(28);
        fitRow[0] = "fit";
        fitRow[1] = std::to_string(crystal);
        fitRow[2] = result.status;
        fitRow[3] = result.needsReview ? "1" : "0";
        fitRow[4] = FormatPrecise(result.p0);
        fitRow[5] = FormatPrecise(result.p1);
        fitRow[6] = FormatPrecise(result.p2);
        fitRow[7] = FormatPrecise(result.chi2);
        fitRow[8] = std::to_string(result.ndf);
        fitRow[9] = FormatPrecise(result.residualRms);
        WriteCsvRow(output, fitRow);
        for (const auto& point : result.points) {
            std::vector<std::string> pointRow(28);
            pointRow[0] = "point";
            pointRow[1] = std::to_string(crystal);
            pointRow[4] = FormatPrecise(result.p0);
            pointRow[5] = FormatPrecise(result.p1);
            pointRow[6] = FormatPrecise(result.p2);
            pointRow[10] = point.datasetId;
            pointRow[11] = FormatPrecise(point.charge);
            pointRow[12] = FormatPrecise(point.energy);
            pointRow[13] = FormatPrecise(point.residual);
            pointRow[14] = point.manual ? "1" : "0";
            pointRow[15] = FormatPrecise(point.chargeError);
            pointRow[16] = FormatPrecise(point.peakFit.rangeLow);
            pointRow[17] = FormatPrecise(point.peakFit.rangeHigh);
            pointRow[18] = FormatPrecise(point.peakFit.sigma);
            pointRow[19] = FormatPrecise(point.peakFit.chi2);
            pointRow[20] = std::to_string(point.peakFit.ndf);
            pointRow[21] = FormatPrecise(point.peakFit.tailFraction);
            pointRow[22] = FormatPrecise(point.peakFit.beta);
            pointRow[23] = FormatPrecise(point.peakFit.stepFraction);
            WriteCsvRow(output, pointRow);
        }
    }
    for (const auto& [datasetId, analysis] : combinedAnalyses_) {
        for (const auto& peak : analysis.peaks) {
            std::vector<std::string> row(28);
            row[0] = "combined_peak";
            row[2] = peak.status;
            row[10] = datasetId;
            row[12] = FormatPrecise(peak.expectedEnergy);
            if (peak.success) {
                row[13] = FormatPrecise(peak.residualKeV);
                row[16] = FormatPrecise(peak.peakFit.rangeLow);
                row[17] = FormatPrecise(peak.peakFit.rangeHigh);
                row[18] = FormatPrecise(peak.peakFit.sigma);
                row[19] = FormatPrecise(peak.peakFit.chi2);
                row[20] = std::to_string(peak.peakFit.ndf);
                row[21] = FormatPrecise(peak.peakFit.tailFraction);
                row[22] = FormatPrecise(peak.peakFit.beta);
                row[23] = FormatPrecise(peak.peakFit.stepFraction);
                row[24] = FormatPrecise(peak.fwhmKeV);
                row[25] = FormatPrecise(peak.resolutionPercent);
                row[27] = FormatPrecise(peak.fittedEnergy);
            }
            row[26] = std::to_string(analysis.crystalCount);
            WriteCsvRow(output, row);
        }
    }
    output.close();
    const std::filesystem::path csvPath(path.toStdString());
    const std::filesystem::path cppPath = csvPath.parent_path() /
        (csvPath.stem().string() + "_coefficients.hpp");
    if (!CalibrationExporter::WriteCppCoefficientLists(cppPath.string(), results_)) {
        SetStatus("Exported CSV, but could not write C++ coefficient lists to " + cppPath.string());
        return;
    }
    SetStatus("Exported CSV to " + path.toStdString() + " and p0/p1/p2 C++ lists to " +
              cppPath.string());
}

void MainWindow::SetStatus(const std::string& text) {
    if (statusLabel_) statusLabel_->setText(Text(text));
}

void MainWindow::SetSecondaryPlotVisible(bool visible) {
    secondaryPlotContainer_->setVisible(visible);
}

} // namespace hpge
