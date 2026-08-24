#include "SpectrumPlotWidget.h"

#include <QApplication>
#include <QImage>
#include <QMouseEvent>
#include <QPainter>
#include <QWheelEvent>

#include <cmath>
#include <iostream>

int main(int argc, char** argv) {
    QApplication application(argc, argv);
    hpge::SpectrumPlotWidget plot;
    plot.resize(900, 520);
    plot.SetPlot("Interaction test", "Charge", "Counts",
                 {{{0.0, 25.0, 50.0, 75.0, 100.0},
                   {0.0, 4.0, 10.0, 4.0, 0.0},
                   "Spectrum", QColor("#2563eb"), 2, false, false}},
                 {{50.0, "centroid", QColor("#dc2626"), false}});
    plot.show();
    application.processEvents();

    const auto initialRange = plot.VisibleXRange();
    if (std::abs(initialRange.first) > 1e-9 || std::abs(initialRange.second - 100.0) > 1e-9) {
        std::cerr << "Initial plot range is incorrect\n";
        return 1;
    }
    plot.ZoomIn();
    const auto zoomedRange = plot.VisibleXRange();
    if (zoomedRange.first < 14.9 || zoomedRange.second > 85.1) {
        std::cerr << "Zoom-in control did not narrow the view\n";
        return 1;
    }

    // Adding a fitted overlay and switching interaction mode must not reset a
    // user's zoomed peak-selection window.
    plot.SetInteractionMode(hpge::SpectrumPlotWidget::InteractionMode::SelectRange);
    plot.SetPlot("Interaction test", "Charge", "Counts",
                 {{{0.0, 25.0, 50.0, 75.0, 100.0},
                   {0.0, 4.0, 10.0, 4.0, 0.0},
                   "Spectrum", QColor("#2563eb"), 2, false, false},
                  {{42.0, 50.0, 58.0}, {3.0, 10.0, 3.0},
                   "Fitted peak", QColor("#dc2626"), 2, false, false}},
                 {{50.0, "centroid", QColor("#dc2626"), false}},
                 {true, 0.0});
    const auto preservedRange = plot.VisibleXRange();
    if (std::abs(preservedRange.first - zoomedRange.first) > 1e-9 ||
        std::abs(preservedRange.second - zoomedRange.second) > 1e-9) {
        std::cerr << "Adding a fit overlay reset the zoomed view\n";
        return 1;
    }

    double selected = -1.0;
    int callbackCount = 0;
    plot.SetRangeSelectedCallback([&](double charge) {
        selected = charge;
        ++callbackCount;
    });
    // The plot rectangle spans x pixels [72, width-24]. Its center remains
    // charge 50 after zooming and overlay refresh.
    const QPointF center((72.0 + 876.0) / 2.0, 250.0);
    QMouseEvent press(QEvent::MouseButtonPress, center, center, center,
                      Qt::LeftButton, Qt::LeftButton, Qt::NoModifier);
    QApplication::sendEvent(&plot, &press);
    QMouseEvent release(QEvent::MouseButtonRelease, center, center, center,
                        Qt::LeftButton, Qt::NoButton, Qt::NoModifier);
    QApplication::sendEvent(&plot, &release);
    if (callbackCount != 1 || std::abs(selected - 50.0) > 0.2) {
        std::cerr << "Selection did not map the click to the expected charge\n";
        return 1;
    }

    const double widthBeforeWheel = preservedRange.second - preservedRange.first;
    QWheelEvent wheel(center, center, QPoint(), QPoint(0, 120), Qt::NoButton,
                      Qt::NoModifier, Qt::NoScrollPhase, false);
    QApplication::sendEvent(&plot, &wheel);
    const auto wheelRange = plot.VisibleXRange();
    if (wheelRange.second - wheelRange.first >= widthBeforeWheel || callbackCount != 1) {
        std::cerr << "Wheel zoom is not available during peak-range selection\n";
        return 1;
    }

    plot.SetInteractionMode(hpge::SpectrumPlotWidget::InteractionMode::ZoomPan);
    QApplication::sendEvent(&plot, &press);
    QApplication::sendEvent(&plot, &release);
    if (callbackCount != 1) {
        std::cerr << "Zoom mode incorrectly emitted a peak-range selection\n";
        return 1;
    }


    plot.SetPlot("Padded fit", "Charge", "Energy",
                 {{{0.0, 100.0}, {0.0, 100.0}, "Fit", QColor("#dc2626"), 2, false, false}},
                 {}, {false, 0.15});
    const auto paddedRange = plot.FullXRange();
    if (std::abs(paddedRange.first + 15.0) > 1e-9 ||
        std::abs(paddedRange.second - 115.0) > 1e-9) {
        std::cerr << "Calibration plot padding was not applied\n";
        return 1;
    }

    QImage image(plot.size(), QImage::Format_ARGB32_Premultiplied);
    image.fill(Qt::transparent);
    QPainter painter(&image);
    plot.render(&painter);
    painter.end();
    if (image.isNull()) {
        std::cerr << "Plot widget did not render\n";
        return 1;
    }
    std::cout << "PASS: persistent ROOT-like zoom, range selection, overlays and padded plots\n";
    return 0;
}
