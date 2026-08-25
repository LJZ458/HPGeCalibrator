#include "MainWindow.h"
#include "SpectrumPlotWidget.h"

#include <QApplication>
#include <QClipboard>
#include <QComboBox>
#include <QDoubleSpinBox>
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
                        : false;
    if (!passed && mode != "alignment-preview" && mode != "result-review" &&
        mode != "selective-refit" && mode != "custom-peak-persistence" &&
        mode != "direct-aligned-fit") {
        std::cerr << "Unknown GUI test mode: " << mode << '\n';
        return 2;
    }
    if (passed) std::cout << "PASS: " << mode << '\n';
    return passed ? 0 : 1;
}
