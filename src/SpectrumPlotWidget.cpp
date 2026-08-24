#include "SpectrumPlotWidget.h"

#include <QFontMetrics>
#include <QMouseEvent>
#include <QPainter>
#include <QPainterPath>
#include <QPalette>
#include <QWheelEvent>

#include <algorithm>
#include <cmath>
#include <limits>

namespace hpge {
namespace {

QString Number(double value) {
    const double magnitude = std::abs(value);
    if ((magnitude > 0.0 && magnitude < 0.01) || magnitude >= 100000.0) {
        return QString::number(value, 'e', 2);
    }
    return QString::number(value, 'f', magnitude < 10.0 ? 2 : magnitude < 1000.0 ? 1 : 0);
}

} // namespace

SpectrumPlotWidget::SpectrumPlotWidget(QWidget* parent) : QWidget(parent) {
    setMinimumSize(520, 320);
    setMouseTracking(true);
    setFocusPolicy(Qt::StrongFocus);
    setAutoFillBackground(true);
}

void SpectrumPlotWidget::SetPlot(std::string title, std::string xLabel,
                                 std::string yLabel, std::vector<PlotSeries> series,
                                 std::vector<PlotMarker> markers) {
    title_ = std::move(title);
    xLabel_ = std::move(xLabel);
    yLabel_ = std::move(yLabel);
    series_ = std::move(series);
    markers_ = std::move(markers);
    emptyMessage_.clear();

    fullMinimum_ = std::numeric_limits<double>::infinity();
    fullMaximum_ = -std::numeric_limits<double>::infinity();
    for (const auto& item : series_) {
        for (double x : item.x) {
            if (!std::isfinite(x)) continue;
            fullMinimum_ = std::min(fullMinimum_, x);
            fullMaximum_ = std::max(fullMaximum_, x);
        }
    }
    if (!std::isfinite(fullMinimum_) || !std::isfinite(fullMaximum_) ||
        fullMaximum_ <= fullMinimum_) {
        fullMinimum_ = 0.0;
        fullMaximum_ = 1.0;
    }
    ResetView();
}

void SpectrumPlotWidget::Clear(const std::string& message) {
    title_.clear();
    series_.clear();
    markers_.clear();
    emptyMessage_ = message;
    fullMinimum_ = viewMinimum_ = 0.0;
    fullMaximum_ = viewMaximum_ = 1.0;
    update();
}

void SpectrumPlotWidget::SetInteractionMode(InteractionMode mode) {
    mode_ = mode;
    dragging_ = false;
    setCursor(mode == InteractionMode::SelectRange ? Qt::CrossCursor : Qt::ArrowCursor);
    update();
}

void SpectrumPlotWidget::ResetView() {
    viewMinimum_ = fullMinimum_;
    viewMaximum_ = fullMaximum_;
    update();
}

void SpectrumPlotWidget::SetRangeSelectedCallback(std::function<void(double)> callback) {
    rangeSelected_ = std::move(callback);
}

QRect SpectrumPlotWidget::PlotRect() const {
    return rect().adjusted(72, 42, -24, -58);
}

double SpectrumPlotWidget::PixelToX(double pixel) const {
    const auto area = PlotRect();
    if (area.width() <= 0) return viewMinimum_;
    const double fraction = (pixel - area.left()) / static_cast<double>(area.width());
    return viewMinimum_ + fraction * (viewMaximum_ - viewMinimum_);
}

double SpectrumPlotWidget::XToPixel(double value) const {
    const auto area = PlotRect();
    return area.left() + (value - viewMinimum_) / (viewMaximum_ - viewMinimum_) * area.width();
}

double SpectrumPlotWidget::YToPixel(double value, double minimum, double maximum) const {
    const auto area = PlotRect();
    return area.bottom() - (value - minimum) / (maximum - minimum) * area.height();
}

std::pair<double, double> SpectrumPlotWidget::VisibleYRange() const {
    double minimum = std::numeric_limits<double>::infinity();
    double maximum = -std::numeric_limits<double>::infinity();
    for (const auto& item : series_) {
        const auto count = std::min(item.x.size(), item.y.size());
        for (std::size_t i = 0; i < count; ++i) {
            if (item.x[i] < viewMinimum_ || item.x[i] > viewMaximum_ ||
                !std::isfinite(item.y[i])) continue;
            minimum = std::min(minimum, item.y[i]);
            maximum = std::max(maximum, item.y[i]);
        }
    }
    if (!std::isfinite(minimum) || !std::isfinite(maximum)) return {0.0, 1.0};
    if (minimum >= 0.0) minimum = 0.0;
    if (maximum <= 0.0) maximum = 0.0;
    if (maximum <= minimum) maximum = minimum + 1.0;
    const double margin = 0.08 * (maximum - minimum);
    return {minimum - (minimum < 0.0 ? margin : 0.0), maximum + margin};
}

void SpectrumPlotWidget::paintEvent(QPaintEvent*) {
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing, true);
    painter.fillRect(rect(), palette().brush(QPalette::Base));
    painter.setPen(palette().color(QPalette::Text));

    if (series_.empty()) {
        painter.drawText(rect(), Qt::AlignCenter,
                         QString::fromStdString(emptyMessage_));
        return;
    }

    const QRect area = PlotRect();
    const auto [yMinimum, yMaximum] = VisibleYRange();
    painter.setFont(QFont(font().family(), font().pointSize() + 1, QFont::DemiBold));
    painter.drawText(QRect(72, 8, width() - 96, 28), Qt::AlignLeft | Qt::AlignVCenter,
                     QString::fromStdString(title_));
    painter.setFont(font());

    const QColor grid = palette().color(QPalette::Midlight);
    for (int tick = 0; tick <= 6; ++tick) {
        const double fraction = tick / 6.0;
        const int x = area.left() + static_cast<int>(fraction * area.width());
        painter.setPen(QPen(grid, 1, Qt::DotLine));
        painter.drawLine(x, area.top(), x, area.bottom());
        painter.setPen(palette().color(QPalette::Text));
        const double value = viewMinimum_ + fraction * (viewMaximum_ - viewMinimum_);
        painter.drawText(QRect(x - 45, area.bottom() + 7, 90, 20), Qt::AlignHCenter,
                         Number(value));
    }
    for (int tick = 0; tick <= 5; ++tick) {
        const double fraction = tick / 5.0;
        const int y = area.bottom() - static_cast<int>(fraction * area.height());
        painter.setPen(QPen(grid, 1, Qt::DotLine));
        painter.drawLine(area.left(), y, area.right(), y);
        painter.setPen(palette().color(QPalette::Text));
        const double value = yMinimum + fraction * (yMaximum - yMinimum);
        painter.drawText(QRect(2, y - 10, 64, 20), Qt::AlignRight | Qt::AlignVCenter,
                         Number(value));
    }
    painter.setPen(QPen(palette().color(QPalette::Text), 1));
    painter.drawRect(area);
    painter.drawText(QRect(area.left(), height() - 29, area.width(), 20), Qt::AlignCenter,
                     QString::fromStdString(xLabel_));
    painter.save();
    painter.translate(18, area.center().y());
    painter.rotate(-90.0);
    painter.drawText(QRect(-area.height() / 2, -10, area.height(), 20), Qt::AlignCenter,
                     QString::fromStdString(yLabel_));
    painter.restore();

    painter.save();
    painter.setClipRect(area.adjusted(1, 1, -1, -1));
    for (const auto& item : series_) {
        QPen pen(item.color, item.width,
                 item.dashed ? Qt::DashLine : Qt::SolidLine,
                 Qt::RoundCap, Qt::RoundJoin);
        painter.setPen(pen);
        const auto count = std::min(item.x.size(), item.y.size());
        QPainterPath path;
        bool started = false;
        for (std::size_t i = 0; i < count; ++i) {
            if (!std::isfinite(item.x[i]) || !std::isfinite(item.y[i])) continue;
            const QPointF point(XToPixel(item.x[i]),
                                YToPixel(item.y[i], yMinimum, yMaximum));
            if (!started) {
                path.moveTo(point);
                started = true;
            } else {
                path.lineTo(point);
            }
            if (item.points && item.x[i] >= viewMinimum_ && item.x[i] <= viewMaximum_) {
                painter.setBrush(item.color);
                painter.drawEllipse(point, 3.5, 3.5);
            }
        }
        if (!item.points) painter.drawPath(path);
    }
    for (const auto& marker : markers_) {
        if (marker.x < viewMinimum_ || marker.x > viewMaximum_) continue;
        const int x = static_cast<int>(XToPixel(marker.x));
        painter.setPen(QPen(marker.color, 2, marker.dashed ? Qt::DashLine : Qt::SolidLine));
        painter.drawLine(x, area.top(), x, area.bottom());
        if (!marker.label.empty()) {
            painter.save();
            painter.translate(x + 4, area.top() + 8);
            painter.rotate(90.0);
            painter.drawText(0, 0, QString::fromStdString(marker.label));
            painter.restore();
        }
    }
    painter.restore();

    int legendX = area.right() - 170;
    int legendY = area.top() + 10;
    for (const auto& item : series_) {
        if (item.name.empty()) continue;
        painter.setPen(QPen(item.color, item.width, item.dashed ? Qt::DashLine : Qt::SolidLine));
        painter.drawLine(legendX, legendY + 7, legendX + 24, legendY + 7);
        painter.setPen(palette().color(QPalette::Text));
        painter.drawText(legendX + 30, legendY + 13, QString::fromStdString(item.name));
        legendY += 19;
    }

    if (dragging_ && mode_ == InteractionMode::ZoomPan && dragButton_ == Qt::LeftButton) {
        QRect selection(dragStart_, dragCurrent_);
        selection = selection.normalized().intersected(area);
        painter.setPen(QPen(QColor("#2563eb"), 1, Qt::DashLine));
        painter.setBrush(QColor(37, 99, 235, 35));
        painter.drawRect(selection);
    }
}

void SpectrumPlotWidget::mousePressEvent(QMouseEvent* event) {
    if (!PlotRect().contains(event->position().toPoint())) return;
    dragStart_ = dragCurrent_ = event->position().toPoint();
    dragButton_ = event->button();
    dragging_ = true;
}

void SpectrumPlotWidget::mouseMoveEvent(QMouseEvent* event) {
    if (!dragging_) return;
    const QPoint previous = dragCurrent_;
    dragCurrent_ = event->position().toPoint();
    if (mode_ == InteractionMode::ZoomPan && dragButton_ == Qt::RightButton) {
        const double shift = -(dragCurrent_.x() - previous.x()) /
                             static_cast<double>(std::max(PlotRect().width(), 1)) *
                             (viewMaximum_ - viewMinimum_);
        viewMinimum_ += shift;
        viewMaximum_ += shift;
        ClampView();
    }
    update();
}

void SpectrumPlotWidget::mouseReleaseEvent(QMouseEvent* event) {
    if (!dragging_) return;
    dragCurrent_ = event->position().toPoint();
    const int distance = (dragCurrent_ - dragStart_).manhattanLength();
    if (mode_ == InteractionMode::SelectRange && event->button() == Qt::LeftButton &&
        PlotRect().contains(dragCurrent_) && distance < 8) {
        if (rangeSelected_) rangeSelected_(PixelToX(dragCurrent_.x()));
    } else if (mode_ == InteractionMode::ZoomPan && event->button() == Qt::LeftButton &&
               std::abs(dragCurrent_.x() - dragStart_.x()) >= 12) {
        const double first = PixelToX(dragStart_.x());
        const double second = PixelToX(dragCurrent_.x());
        viewMinimum_ = std::min(first, second);
        viewMaximum_ = std::max(first, second);
        ClampView();
    }
    dragging_ = false;
    dragButton_ = Qt::NoButton;
    update();
}

void SpectrumPlotWidget::mouseDoubleClickEvent(QMouseEvent*) { ResetView(); }

void SpectrumPlotWidget::wheelEvent(QWheelEvent* event) {
    if (mode_ != InteractionMode::ZoomPan || series_.empty()) return;
    const double anchor = PixelToX(event->position().x());
    const double factor = event->angleDelta().y() > 0 ? 0.80 : 1.25;
    viewMinimum_ = anchor + (viewMinimum_ - anchor) * factor;
    viewMaximum_ = anchor + (viewMaximum_ - anchor) * factor;
    ClampView();
    update();
    event->accept();
}

void SpectrumPlotWidget::ClampView() {
    const double fullWidth = fullMaximum_ - fullMinimum_;
    double width = viewMaximum_ - viewMinimum_;
    const double minimumWidth = fullWidth * 1e-5;
    if (width < minimumWidth) {
        const double center = 0.5 * (viewMinimum_ + viewMaximum_);
        viewMinimum_ = center - 0.5 * minimumWidth;
        viewMaximum_ = center + 0.5 * minimumWidth;
        width = minimumWidth;
    }
    if (width >= fullWidth) {
        viewMinimum_ = fullMinimum_;
        viewMaximum_ = fullMaximum_;
        return;
    }
    if (viewMinimum_ < fullMinimum_) {
        viewMinimum_ = fullMinimum_;
        viewMaximum_ = fullMinimum_ + width;
    }
    if (viewMaximum_ > fullMaximum_) {
        viewMaximum_ = fullMaximum_;
        viewMinimum_ = fullMaximum_ - width;
    }
}

} // namespace hpge
