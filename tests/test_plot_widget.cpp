#include "SpectrumPlotWidget.h"

#include <QApplication>
#include <QImage>
#include <QMouseEvent>
#include <QPainter>

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

    double selected = -1.0;
    int callbackCount = 0;
    plot.SetRangeSelectedCallback([&](double charge) {
        selected = charge;
        ++callbackCount;
    });
    plot.SetInteractionMode(hpge::SpectrumPlotWidget::InteractionMode::SelectRange);

    // The plot rectangle spans x pixels [72, width-24]. Its center maps to charge 50.
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

    plot.SetInteractionMode(hpge::SpectrumPlotWidget::InteractionMode::ZoomPan);
    QApplication::sendEvent(&plot, &press);
    QApplication::sendEvent(&plot, &release);
    if (callbackCount != 1) {
        std::cerr << "Zoom mode incorrectly emitted a peak-range selection\n";
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
    std::cout << "PASS: Qt plot rendering and separated selection/zoom interaction\n";
    return 0;
}
