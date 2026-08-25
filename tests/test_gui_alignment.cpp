#include "MainWindow.h"
#include "SpectrumPlotWidget.h"

#include <QApplication>
#include <QCheckBox>
#include <QClipboard>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QFile>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QPushButton>
#include <QSettings>
#include <QTemporaryDir>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

namespace hpge {

class MainWindowTestAccess {
public:
    static bool InstallTwoSourceResult(MainWindow& window) {
        if (window.descriptors_.size() < 2) return false;
        CalibrationResult result;
        result.crystal = 3;
        result.success = true;
        result.needsReview = false;
        result.status = "test result";
        for (std::size_t source = 0; source < 2; ++source) {
            for (int peak = 0; peak < 2; ++peak) {
                CalibrationPoint point;
                point.datasetId = window.descriptors_[source].id;
                point.energy = 800.0 + 500.0 * peak + 100.0 * source;
                point.charge = 1050.0 + 650.0 * peak + 90.0 * source;
                point.peakFit.success = true;
                point.peakFit.rangeLow = point.charge - 24.0;
                point.peakFit.rangeHigh = point.charge + 24.0;
                point.peakFit.centroid = point.charge;
                point.peakFit.sigma = 3.5;
                point.peakFit.height = 400.0;
                point.peakFit.background0 = 3.0;
                result.points.push_back(point);
            }
        }
        window.results_.clear();
        window.results_[result.crystal] = std::move(result);
        window.RefreshResults();
        return true;
    }

    static bool ApplySelectiveReplacementAndAddition(MainWindow& window) {
        const int crystal = window.CurrentResultCrystal();
        const auto existing = window.results_.find(crystal);
        if (existing == window.results_.end() || existing->second.points.size() != 4) return false;
        const CalibrationPoint replacedPoint = existing->second.points.front();
        const std::vector<CalibrationPoint> untouched(
            std::next(existing->second.points.begin()), existing->second.points.end());
        MainWindow::ManualPeak replacement;
        replacement.datasetId = replacedPoint.datasetId;
        replacement.crystal = crystal;
        replacement.energy = replacedPoint.energy;
        replacement.charge = replacedPoint.charge + 4.25;
        replacement.label = "replacement";
        replacement.peakFit = replacedPoint.peakFit;
        replacement.peakFit.centroid = replacement.charge;
        replacement.peakFit.rangeLow += 4.25;
        replacement.peakFit.rangeHigh += 4.25;

        MainWindow::ManualPeak addition;
        addition.datasetId = window.descriptors_[1].id;
        addition.crystal = crystal;
        addition.energy = 2000.0;
        addition.charge = 2550.0;
        addition.label = "addition";
        addition.peakFit.success = true;
        addition.peakFit.rangeLow = 2525.0;
        addition.peakFit.rangeHigh = 2575.0;
        addition.peakFit.centroid = addition.charge;
        addition.peakFit.sigma = 4.0;
        addition.peakFit.height = 300.0;
        addition.peakFit.background0 = 2.0;

        window.manualPeaks_.push_back(std::move(replacement));
        window.manualPeaks_.push_back(std::move(addition));
        window.RefreshManualPeakList();
        window.RefitSelectedCrystal();
        const auto updated = window.results_.find(crystal);
        if (updated == window.results_.end() || updated->second.points.size() != 5) return false;
        for (const auto& original : untouched) {
            const auto preserved = std::find_if(
                updated->second.points.begin(), updated->second.points.end(),
                [&](const CalibrationPoint& point) {
                    return point.datasetId == original.datasetId &&
                           std::abs(point.energy - original.energy) < 1e-12;
                });
            if (preserved == updated->second.points.end() ||
                preserved->charge != original.charge ||
                preserved->chargeError != original.chargeError ||
                preserved->peakFit.centroid != original.peakFit.centroid ||
                preserved->peakFit.rangeLow != original.peakFit.rangeLow ||
                preserved->peakFit.rangeHigh != original.peakFit.rangeHigh ||
                preserved->peakFit.sigma != original.peakFit.sigma ||
                preserved->peakFit.height != original.peakFit.height) {
                return false;
            }
        }
        int manualCount = 0;
        for (const auto& point : updated->second.points) manualCount += point.manual ? 1 : 0;
        return manualCount == 2 && std::none_of(
            window.manualPeaks_.begin(), window.manualPeaks_.end(),
            [crystal](const MainWindow::ManualPeak& peak) { return peak.crystal == crystal; });
    }

    static bool RunDirectAlignedRangeCalibration(MainWindow& window, int crystal,
                                                 double halfRange) {
        const auto referenceCharge = [](double energy) {
            constexpr double p0 = -1.5;
            constexpr double p1 = 0.68;
            constexpr double p2 = 1.8e-5;
            return (-p1 + std::sqrt(p1 * p1 - 4.0 * p2 * (p0 - energy))) /
                   (2.0 * p2);
        };
        window.referenceCrystalEntry_->setValue(0);
        window.histogramList_->selectAll();
        window.alignedFitHalfRangeEntry_->setValue(halfRange);
        window.referencePeaks_.clear();
        for (const auto& descriptor : window.descriptors_) {
            const bool co60 = descriptor.objectPath.find("co60") != std::string::npos;
            const std::vector<double> energies = co60
                ? std::vector<double>{1173.228, 1332.492}
                : std::vector<double>{846.771, 1238.282, 1771.351, 2598.459};
            for (double energy : energies) {
                const double charge = referenceCharge(energy);
                PeakFitResult fit;
                fit.success = true;
                fit.rangeLow = charge - 25.0;
                fit.rangeHigh = charge + 25.0;
                fit.centroid = charge;
                fit.centroidError = 0.2;
                fit.sigma = 4.0;
                fit.height = 1000.0;
                window.referencePeaks_.push_back(
                    {descriptor.id, charge, energy, "test reference", fit});
            }
        }
        window.results_.clear();
        window.alignmentResults_.clear();
        window.results_[crystal] = window.CalibrateCrystal(crystal);
        const auto result = window.results_.find(crystal);
        if (result == window.results_.end() ||
            result->second.points.size() != window.referencePeaks_.size()) return false;
        for (const auto& point : result->second.points) {
            const auto alignment = window.alignmentResults_.find({crystal, point.datasetId});
            const auto reference = std::find_if(
                window.referencePeaks_.begin(), window.referencePeaks_.end(),
                [&](const ReferencePeak& peak) {
                    return peak.datasetId == point.datasetId &&
                           std::abs(peak.energy - point.energy) < 1e-6;
                });
            if (alignment == window.alignmentResults_.end() ||
                !alignment->second.success || reference == window.referencePeaks_.end()) {
                return false;
            }
            double expectedLow = CalibrationEngine::MapReferenceCharge(
                alignment->second, reference->peakFit.centroid - halfRange);
            double expectedHigh = CalibrationEngine::MapReferenceCharge(
                alignment->second, reference->peakFit.centroid + halfRange);
            if (expectedHigh < expectedLow) std::swap(expectedLow, expectedHigh);
            if (std::abs(point.peakFit.rangeLow - expectedLow) > 1e-9 ||
                std::abs(point.peakFit.rangeHigh - expectedHigh) > 1e-9) {
                return false;
            }
        }
        window.RefreshResults();
        return true;
    }

    static bool PrepareCompleteProject(MainWindow& window) {
        if (!InstallTwoSourceResult(window) ||
            !ApplySelectiveReplacementAndAddition(window)) return false;
        window.referencePeaks_.clear();
        for (std::size_t dataset = 0; dataset < window.descriptors_.size(); ++dataset) {
            const auto result = window.results_.find(3);
            if (result == window.results_.end()) return false;
            const auto point = std::find_if(
                result->second.points.begin(), result->second.points.end(),
                [&](const CalibrationPoint& item) {
                    return item.datasetId == window.descriptors_[dataset].id;
                });
            if (point == result->second.points.end()) return false;
            window.referencePeaks_.push_back({point->datasetId, point->charge,
                point->energy, "project reference", point->peakFit});
            PeakMatchResult alignment;
            alignment.success = true;
            alignment.referenceCharges = {500.0, 1500.0, 2500.0};
            alignment.charges = {525.0, 1575.0, 2625.0};
            alignment.matched = {true, true, true};
            alignment.scale = 1.05 + 0.01 * static_cast<double>(dataset);
            alignment.offset = 12.0 + static_cast<double>(dataset);
            alignment.quadratic = 1.0e-6 * static_cast<double>(dataset);
            alignment.alignmentCost = 2.5 + static_cast<double>(dataset);
            alignment.referenceSensitivity = 0.42;
            alignment.targetSensitivity = 0.37;
            alignment.quadraticModel = dataset != 0;
            window.alignmentResults_[{3, point->datasetId}] = std::move(alignment);
        }
        MainWindow::ManualPeak pending;
        pending.datasetId = window.descriptors_.front().id;
        pending.crystal = 3;
        pending.charge = 2875.25;
        pending.energy = 2222.75;
        pending.label = "pending manual project peak";
        pending.peakFit.success = true;
        pending.peakFit.status = "manual test fit";
        pending.peakFit.rangeLow = 2850.0;
        pending.peakFit.rangeHigh = 2900.0;
        pending.peakFit.centroid = pending.charge;
        pending.peakFit.centroidError = 0.35;
        pending.peakFit.sigma = 4.2;
        pending.peakFit.height = 275.0;
        pending.peakFit.tailFraction = 0.05;
        pending.peakFit.beta = 15.0;
        pending.peakFit.stepFraction = 0.02;
        pending.peakFit.background0 = 3.0;
        pending.peakFit.background1 = 0.2;
        pending.peakFit.background2 = 0.01;
        pending.peakFit.chi2 = 18.5;
        pending.peakFit.ndf = 31;
        window.manualPeaks_.push_back(std::move(pending));
        window.energyLines_.push_back({2222.75, "project-only custom", "Custom", false});
        window.PopulateEnergyLines();
        window.histogramList_->selectAll();
        window.crystalList_->clearSelection();
        window.crystalList_->item(3)->setSelected(true);
        window.crystalList_->item(7)->setSelected(true);
        window.orientationCombo_->setCurrentIndex(0);
        window.referenceCrystalEntry_->setValue(0);
        window.sigmaEntry_->setValue(3.4);
        window.thresholdEntry_->setValue(0.027);
        window.residualLimitEntry_->setValue(2.75);
        window.alignmentSensitivityEntry_->setValue(58.0);
        window.autoTuneAlignmentEntry_->setChecked(false);
        window.alignmentModelCombo_->setCurrentIndex(2);
        window.alignedFitHalfRangeEntry_->setValue(42.5);
        window.alignmentCrystalEntry_->setValue(7);
        window.RefreshReferencePeakList();
        window.RefreshManualPeakList();
        return true;
    }

    static bool VerifyCompleteProject(const MainWindow& window) {
        if (window.descriptors_.size() != 2 ||
            window.SelectedDescriptorIndices().size() != 2 ||
            window.SelectedCrystals() != std::vector<int>({3, 7}) ||
            window.referencePeaks_.size() != 2 || window.manualPeaks_.size() != 1 ||
            window.results_.size() != 1 || window.alignmentResults_.size() != 2) return false;
        const auto result = window.results_.find(3);
        if (result == window.results_.end() || result->second.points.size() != 5) return false;
        const int appliedManual = static_cast<int>(std::count_if(
            result->second.points.begin(), result->second.points.end(),
            [](const CalibrationPoint& point) { return point.manual; }));
        const auto& pending = window.manualPeaks_.front();
        const bool customRestored = std::any_of(
            window.energyLines_.begin(), window.energyLines_.end(),
            [](const MainWindow::EnergyLine& line) {
                return line.source == "Custom" && std::abs(line.energy - 2222.75) < 1e-9 &&
                       line.label == "project-only custom";
            });
        return appliedManual == 2 && customRestored &&
               pending.label == "pending manual project peak" &&
               pending.peakFit.status == "manual test fit" &&
               pending.peakFit.rangeLow == 2850.0 &&
               pending.peakFit.rangeHigh == 2900.0 &&
               pending.peakFit.centroid == 2875.25 &&
               pending.peakFit.sigma == 4.2 && pending.peakFit.ndf == 31 &&
               window.sigmaEntry_->value() == 3.4 &&
               window.thresholdEntry_->value() == 0.027 &&
               window.residualLimitEntry_->value() == 2.75 &&
               window.alignmentSensitivityEntry_->value() == 58.0 &&
               !window.autoTuneAlignmentEntry_->isChecked() &&
               window.alignmentModelCombo_->currentIndex() == 2 &&
               window.alignedFitHalfRangeEntry_->value() == 42.5 &&
               window.alignmentCrystalEntry_->value() == 7;
    }

    static bool ContinueRestoredProject(MainWindow& window,
                                        const std::string& additionalRootFile) {
        const std::size_t previousDescriptors = window.descriptors_.size();
        const std::size_t previousReferences = window.referencePeaks_.size();
        const std::size_t previousResults = window.results_.size();
        ReferencePeak extra = window.referencePeaks_.front();
        extra.energy += 0.5;
        extra.label = "added after restore";
        window.referencePeaks_.push_back(std::move(extra));
        if (!window.OpenRootFiles({additionalRootFile})) return false;
        return window.descriptors_.size() > previousDescriptors &&
               window.referencePeaks_.size() == previousReferences + 1 &&
               window.results_.size() == previousResults;
    }

    static bool VerifyContinuedProject(const MainWindow& window) {
        return window.descriptors_.size() == 4 && window.referencePeaks_.size() == 3 &&
               window.results_.size() == 1 && std::any_of(
                   window.referencePeaks_.begin(), window.referencePeaks_.end(),
                   [](const ReferencePeak& peak) { return peak.label == "added after restore"; });
    }

    static void InstallResidualOverview(MainWindow& window) {
        window.results_.clear();
        for (int crystal = 0; crystal < 64; ++crystal) {
            CalibrationResult result;
            result.crystal = crystal;
            result.success = true;
            result.needsReview = false;
            result.status = "residual overview test";
            for (std::size_t dataset = 0; dataset < window.descriptors_.size(); ++dataset) {
                for (int line = 0; line < 2; ++line) {
                    CalibrationPoint point;
                    point.datasetId = window.descriptors_[dataset].id;
                    point.energy = 800.0 + 500.0 * line + 100.0 * dataset;
                    point.charge = 1000.0 + 600.0 * line + 5.0 * crystal;
                    point.residual = 0.04 * static_cast<double>((crystal % 9) - 4) +
                                     0.12 * line - 0.05 * static_cast<double>(dataset);
                    point.peakFit.success = true;
                    point.peakFit.centroid = point.charge;
                    result.points.push_back(point);
                }
            }
            window.results_[crystal] = std::move(result);
        }
        window.RefreshResults();
    }
};

} // namespace hpge

namespace {

bool TestAlignmentPreview(hpge::MainWindow& window, QApplication& application) {
    auto* histogram = window.findChild<QComboBox*>("alignmentHistogramCombo");
    auto* preview = window.findChild<QPushButton*>("showAlignmentButton");
    auto* plot = dynamic_cast<hpge::SpectrumPlotWidget*>(
        window.findChild<QWidget*>("primaryPlot"));
    auto* status = window.findChild<QLabel*>("statusLabel");
    auto* parameters = window.findChild<QLabel*>("alignmentPreviewParameters");
    if (!histogram || !preview || !plot || !status || !parameters || histogram->count() < 2) {
        std::cerr << "Alignment preview controls or multiple histograms are missing\n";
        return false;
    }

    histogram->setCurrentIndex(0);
    preview->click();
    application.processEvents();
    const std::string firstName = histogram->currentText().toStdString();
    if (plot->Title().find(firstName) == std::string::npos || plot->Series().size() != 2) {
        std::cerr << "First alignment preview does not identify or plot its histogram\n";
        return false;
    }
    const std::vector<double> firstReference = plot->Series().front().y;

    histogram->setCurrentIndex(1);
    application.processEvents();
    const std::string secondName = histogram->currentText().toStdString();
    if (secondName == firstName || plot->Title().find(secondName) == std::string::npos ||
        status->text().toStdString().find(secondName) == std::string::npos ||
        plot->Series().size() != 2) {
        std::cerr << "Changing the alignment histogram did not refresh the preview\n";
        return false;
    }
    const auto& secondReference = plot->Series().front().y;
    if (secondReference.size() != firstReference.size()) {
        return true;
    }
    double difference = 0.0;
    for (std::size_t index = 0; index < firstReference.size(); ++index) {
        difference += std::abs(firstReference[index] - secondReference[index]);
    }
    if (difference <= 1e-6) {
        std::cerr << "The refreshed preview still contains the first histogram's data\n";
        return false;
    }
    return parameters->text().contains("a0=") &&
           parameters->text().contains("a1=") &&
           parameters->text().contains("cost=");
}

bool TestMultipleSourceResultReview(hpge::MainWindow& window,
                                    QApplication& application) {
    if (!hpge::MainWindowTestAccess::InstallTwoSourceResult(window)) {
        std::cerr << "Could not install the two-source calibration result\n";
        return false;
    }
    auto* sources = window.findChild<QComboBox*>("resultSpectrumCombo");
    auto* showAll = window.findChild<QPushButton*>("showAllResultSpectraButton");
    auto* plot = dynamic_cast<hpge::SpectrumPlotWidget*>(
        window.findChild<QWidget*>("primaryPlot"));
    auto* status = window.findChild<QLabel*>("statusLabel");
    if (!sources || !showAll || !plot || !status || sources->count() != 2) {
        std::cerr << "Result review did not list both fitted source spectra\n";
        return false;
    }
    showAll->click();
    application.processEvents();
    int redFits = 0;
    for (const auto& series : plot->Series()) {
        if (series.color == QColor("#dc2626")) ++redFits;
    }
    if (plot->Title().find("All fitted source spectra") == std::string::npos ||
        plot->Series().size() < 6 || redFits < 4 ||
        status->text().toStdString().find("all 2 fitted source spectra") == std::string::npos) {
        std::cerr << "All-source result review did not show both spectra and their red fits\n";
        return false;
    }

    sources->setCurrentIndex(1);
    application.processEvents();
    const std::string selectedName = sources->currentText().toStdString();
    if (plot->Title().find(selectedName) == std::string::npos ||
        plot->Series().size() < 3 ||
        status->text().toStdString().find("source 2/2") == std::string::npos) {
        std::cerr << "Selected-source result review did not switch spectra\n";
        return false;
    }
    return true;
}

bool TestSelectivePeakRefit(hpge::MainWindow& window, QApplication& application) {
    if (!hpge::MainWindowTestAccess::InstallTwoSourceResult(window)) {
        std::cerr << "Could not install the selective-refit source result\n";
        return false;
    }
    auto* fittedPoints = window.findChild<QListWidget*>("fittedPointList");
    auto* status = window.findChild<QLabel*>("statusLabel");
    if (!fittedPoints || !status || fittedPoints->count() != 4) {
        std::cerr << "Existing fitted peaks are not exposed for selective review\n";
        return false;
    }
    if (!hpge::MainWindowTestAccess::ApplySelectiveReplacementAndAddition(window)) {
        std::cerr << "Selective replacement/addition did not preserve the result point set\n";
        return false;
    }
    application.processEvents();
    int manualItems = 0;
    for (int row = 0; row < fittedPoints->count(); ++row) {
        if (fittedPoints->item(row)->text().startsWith("[manual]")) ++manualItems;
    }
    if (fittedPoints->count() != 5 || manualItems != 2 ||
        status->text().toStdString().find("1 peak(s) replaced, 1 added") == std::string::npos) {
        std::cerr << "Selective-refit review did not identify replaced and added peaks\n";
        return false;
    }
    return true;
}

bool TestCustomPeakPersistence(hpge::MainWindow& window,
                               QApplication& application) {
    auto* energy = window.findChild<QDoubleSpinBox*>("customEnergyEntry");
    auto* label = window.findChild<QLineEdit*>("customLineLabelEntry");
    auto* save = window.findChild<QPushButton*>("saveCustomLineButton");
    auto* add = window.findChild<QPushButton*>("addCustomLineButton");
    auto* source = window.findChild<QComboBox*>("referenceSourceCombo");
    auto* lines = window.findChild<QListWidget*>("referenceEnergyList");
    if (!energy || !label || !save || !add || !source || !lines) {
        std::cerr << "Custom peak persistence controls are missing\n";
        return false;
    }
    energy->setValue(1115.432);
    label->setText("laboratory background");
    save->click();
    energy->setValue(1777.111);
    label->setText("session only");
    add->click();
    application.processEvents();
    source->setCurrentText("Custom");
    application.processEvents();
    if (lines->count() != 2 ||
        !lines->item(0)->text().contains("laboratory background") ||
        !lines->item(0)->text().contains("[saved]") ||
        lines->item(1)->text().contains("[saved]")) {
        std::cerr << "Saved and session-only custom peaks are not distinguished\n";
        return false;
    }

    {
        hpge::MainWindow restored;
        auto* restoredSource = restored.findChild<QComboBox*>("referenceSourceCombo");
        auto* restoredLines = restored.findChild<QListWidget*>("referenceEnergyList");
        auto* remove = restored.findChild<QPushButton*>("removeCustomLineButton");
        if (!restoredSource || !restoredLines || !remove) return false;
        restoredSource->setCurrentText("Custom");
        application.processEvents();
        if (restoredLines->count() != 1 ||
            !restoredLines->item(0)->text().contains("1115.432") ||
            !restoredLines->item(0)->text().contains("laboratory background") ||
            !restoredLines->item(0)->text().contains("[saved]")) {
            std::cerr << "Saved custom peak was not restored in a new window\n";
            return false;
        }
        restoredLines->setCurrentRow(0);
        remove->click();
        application.processEvents();
        if (restoredLines->count() != 0) {
            std::cerr << "Saved custom peak was not removed\n";
            return false;
        }
    }

    hpge::MainWindow afterRemoval;
    auto* finalSource = afterRemoval.findChild<QComboBox*>("referenceSourceCombo");
    auto* finalLines = afterRemoval.findChild<QListWidget*>("referenceEnergyList");
    if (!finalSource || !finalLines) return false;
    finalSource->setCurrentText("Custom");
    application.processEvents();
    if (finalLines->count() != 0) {
        std::cerr << "Removed custom peak returned in a later window\n";
        return false;
    }
    return true;
}

bool TestDirectAlignedFitRange(hpge::MainWindow& window,
                               QApplication& application) {
    constexpr double halfRange = 18.0;
    if (!hpge::MainWindowTestAccess::RunDirectAlignedRangeCalibration(
            window, 3, halfRange)) {
        std::cerr << "Automatic fit did not use the exact mapped alignment interval\n";
        return false;
    }
    application.processEvents();
    auto* range = window.findChild<QDoubleSpinBox*>("alignedFitHalfRangeEntry");
    auto* parameters = window.findChild<QListWidget*>("alignmentParameterList");
    auto* copy = window.findChild<QPushButton*>("copyAlignmentParametersButton");
    if (!range || !parameters || !copy || std::abs(range->value() - halfRange) > 1e-12 ||
        parameters->count() != 2) {
        std::cerr << "Aligned fit-range or stored parameter controls are unavailable\n";
        return false;
    }
    for (int row = 0; row < parameters->count(); ++row) {
        const QString text = parameters->item(row)->text();
        if (!text.contains("a0=") || !text.contains("a1=") || !text.contains("a2=") ||
            !text.contains("cost=") || !text.contains("matched=")) {
            std::cerr << "Stored alignment row omits accessible parameters\n";
            return false;
        }
    }
    parameters->setCurrentRow(1);
    copy->click();
    application.processEvents();
    if (QApplication::clipboard()->text() != parameters->currentItem()->text()) {
        std::cerr << "Alignment parameters were not copied to the clipboard\n";
        return false;
    }
    return true;
}

bool TestCompleteProjectPersistence(hpge::MainWindow& window,
                                    QApplication& application,
                                    const std::string& samplePath) {
    if (!hpge::MainWindowTestAccess::PrepareCompleteProject(window)) {
        std::cerr << "Could not prepare complete fitting state for project save\n";
        return false;
    }
    QTemporaryDir projectDirectory;
    if (!projectDirectory.isValid()) return false;
    const std::string projectPath =
        (projectDirectory.path() + "/calibration.hpgecal.json").toStdString();
    std::string error;
    if (!window.SaveProject(projectPath, error)) {
        std::cerr << "Could not save complete project: " << error << '\n';
        return false;
    }
    hpge::MainWindow restored;
    if (!restored.OpenProject(projectPath, error) ||
        !hpge::MainWindowTestAccess::VerifyCompleteProject(restored)) {
        std::cerr << "Complete project was not restored: " << error << '\n';
        return false;
    }
    application.processEvents();

    const QString additionalPath = projectDirectory.path() + "/additional.root";
    if (!QFile::copy(QString::fromStdString(samplePath), additionalPath) ||
        !hpge::MainWindowTestAccess::ContinueRestoredProject(
            restored, additionalPath.toStdString())) {
        std::cerr << "Could not add peaks and ROOT files after project restoration\n";
        return false;
    }
    const std::string resumedPath =
        (projectDirectory.path() + "/resumed.hpgecal.json").toStdString();
    if (!restored.SaveProject(resumedPath, error)) return false;
    hpge::MainWindow resumed;
    if (!resumed.OpenProject(resumedPath, error) ||
        !hpge::MainWindowTestAccess::VerifyContinuedProject(resumed)) {
        std::cerr << "Extended project did not survive a second restore: " << error << '\n';
        return false;
    }
    return true;
}

bool TestResidualsByCrystal(hpge::MainWindow& window,
                            QApplication& application) {
    hpge::MainWindowTestAccess::InstallResidualOverview(window);
    auto* dataset = window.findChild<QComboBox*>("residualDatasetCombo");
    auto* show = window.findChild<QPushButton*>("showResidualsByCrystalButton");
    auto* primary = dynamic_cast<hpge::SpectrumPlotWidget*>(
        window.findChild<QWidget*>("primaryPlot"));
    if (!dataset || !show || !primary || dataset->count() != 3) {
        std::cerr << "Residual-overview controls are unavailable\n";
        return false;
    }
    dataset->setCurrentIndex(0);
    show->click();
    application.processEvents();
    const auto [minimum, maximum] = primary->FullXRange();
    if (primary->Title().find("residuals versus detector crystal") == std::string::npos ||
        primary->Series().size() != 5 || minimum != 0.0 || maximum != 63.0) {
        std::cerr << "All-source residual overview does not span crystals 0-63\n";
        return false;
    }
    for (std::size_t index = 1; index < primary->Series().size(); ++index) {
        if (primary->Series()[index].x.size() != 64 ||
            primary->Series()[index].y.size() != 64) {
            std::cerr << "A fitted energy does not include all 64 crystal residuals\n";
            return false;
        }
    }
    dataset->setCurrentIndex(1);
    show->click();
    application.processEvents();
    if (primary->Series().size() != 3) {
        std::cerr << "Dataset-filtered residual overview has the wrong energy series\n";
        return false;
    }
    return true;
}

} // namespace

int main(int argc, char** argv) {
    QApplication application(argc, argv);
    if (argc < 2 || argc > 3) {
        std::cerr << "Expected a ROOT sample path and optional test mode\n";
        return 2;
    }

    const std::string mode = argc == 3 ? argv[2] : "alignment-preview";
    QTemporaryDir settingsDirectory;
    if (!settingsDirectory.isValid()) {
        std::cerr << "Could not create isolated GUI-test settings directory\n";
        return 1;
    }
    QCoreApplication::setOrganizationName("HPGeCalibratorTests");
    QCoreApplication::setApplicationName("HPGeCalibratorGuiTests");
    QSettings::setDefaultFormat(QSettings::IniFormat);
    QSettings::setPath(QSettings::IniFormat, QSettings::UserScope,
                       settingsDirectory.path());

    hpge::MainWindow window;
    if (!window.OpenRootFiles({argv[1]})) {
        std::cerr << "Could not load the multi-histogram ROOT sample\n";
        return 1;
    }
    const bool passed = mode == "alignment-preview"
        ? TestAlignmentPreview(window, application)
        : mode == "result-review"
            ? TestMultipleSourceResultReview(window, application)
            : mode == "selective-refit"
                ? TestSelectivePeakRefit(window, application)
                : mode == "custom-peak-persistence"
                    ? TestCustomPeakPersistence(window, application)
                    : mode == "direct-aligned-fit"
                        ? TestDirectAlignedFitRange(window, application)
                        : mode == "project-persistence"
                            ? TestCompleteProjectPersistence(
                                  window, application, argv[1])
                            : mode == "residuals-by-crystal"
                                ? TestResidualsByCrystal(window, application)
                                : false;
    if (!passed && mode != "alignment-preview" && mode != "result-review" &&
        mode != "selective-refit" && mode != "custom-peak-persistence" &&
        mode != "direct-aligned-fit" && mode != "project-persistence" &&
        mode != "residuals-by-crystal") {
        std::cerr << "Unknown GUI test mode: " << mode << '\n';
        return 2;
    }
    if (passed) std::cout << "PASS: " << mode << '\n';
    return passed ? 0 : 1;
}
