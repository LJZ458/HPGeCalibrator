#include "MainWindow.h"
#include "SpectrumPlotWidget.h"

#include <QApplication>
#include <QComboBox>
#include <QLabel>
#include <QPushButton>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    QApplication application(argc, argv);
    if (argc != 2) {
        std::cerr << "Expected one ROOT sample path\n";
        return 2;
    }

    hpge::MainWindow window;
    if (!window.OpenRootFiles({argv[1]})) {
        std::cerr << "Could not load the multi-histogram ROOT sample\n";
        return 1;
    }
    auto* histogram = window.findChild<QComboBox*>("alignmentHistogramCombo");
    auto* preview = window.findChild<QPushButton*>("showAlignmentButton");
    auto* plot = dynamic_cast<hpge::SpectrumPlotWidget*>(
        window.findChild<QWidget*>("primaryPlot"));
    auto* status = window.findChild<QLabel*>("statusLabel");
    if (!histogram || !preview || !plot || !status || histogram->count() < 2) {
        std::cerr << "Alignment preview controls or multiple histograms are missing\n";
        return 1;
    }

    histogram->setCurrentIndex(0);
    preview->click();
    application.processEvents();
    const std::string firstName = histogram->currentText().toStdString();
    if (plot->Title().find(firstName) == std::string::npos || plot->Series().size() != 2) {
        std::cerr << "First alignment preview does not identify or plot its histogram\n";
        return 1;
    }
    const std::vector<double> firstReference = plot->Series().front().y;

    histogram->setCurrentIndex(1);
    application.processEvents();
    const std::string secondName = histogram->currentText().toStdString();
    if (secondName == firstName || plot->Title().find(secondName) == std::string::npos ||
        status->text().toStdString().find(secondName) == std::string::npos ||
        plot->Series().size() != 2) {
        std::cerr << "Changing the alignment histogram did not refresh the preview\n";
        return 1;
    }
    const auto& secondReference = plot->Series().front().y;
    if (secondReference.size() != firstReference.size()) {
        std::cout << "PASS: alignment-multiple-histogram-preview\n";
        return 0;
    }
    double difference = 0.0;
    for (std::size_t index = 0; index < firstReference.size(); ++index) {
        difference += std::abs(firstReference[index] - secondReference[index]);
    }
    if (difference <= 1e-6) {
        std::cerr << "The refreshed preview still contains the first histogram's data\n";
        return 1;
    }
    std::cout << "PASS: alignment-multiple-histogram-preview\n";
    return 0;
}
