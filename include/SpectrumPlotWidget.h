#pragma once

#include <QColor>
#include <QPoint>
#include <QRect>
#include <QWidget>

#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace hpge {

struct PlotSeries {
    std::vector<double> x;
    std::vector<double> y;
    std::string name;
    QColor color = QColor("#2563eb");
    int width = 2;
    bool dashed = false;
    bool points = false;
};

struct PlotMarker {
    double x = 0.0;
    std::string label;
    QColor color = QColor("#dc2626");
    bool dashed = false;
};

struct PlotViewOptions {
    bool preserveView = false;
    double horizontalPaddingFraction = 0.0;
};

class SpectrumPlotWidget final : public QWidget {
public:
    enum class InteractionMode { ZoomPan, SelectRange };

    explicit SpectrumPlotWidget(QWidget* parent = nullptr);

    void SetPlot(std::string title, std::string xLabel, std::string yLabel,
                 std::vector<PlotSeries> series,
                 std::vector<PlotMarker> markers = {},
                 PlotViewOptions options = {});
    void Clear(const std::string& message = "No spectrum loaded");
    void SetInteractionMode(InteractionMode mode);
    InteractionMode GetInteractionMode() const { return mode_; }
    void ResetView();
    void ZoomIn();
    void ZoomOut();
    void PreviousView();
    std::pair<double, double> VisibleXRange() const;
    std::pair<double, double> FullXRange() const;
    void SetRangeSelectedCallback(std::function<void(double)> callback);

protected:
    void paintEvent(QPaintEvent* event) override;
    void mousePressEvent(QMouseEvent* event) override;
    void mouseMoveEvent(QMouseEvent* event) override;
    void mouseReleaseEvent(QMouseEvent* event) override;
    void mouseDoubleClickEvent(QMouseEvent* event) override;
    void wheelEvent(QWheelEvent* event) override;

private:
    QRect PlotRect() const;
    double PixelToX(double pixel) const;
    double XToPixel(double value) const;
    double YToPixel(double value, double minimum, double maximum) const;
    std::pair<double, double> VisibleYRange() const;
    void ZoomAround(double factor, double anchor, bool remember);
    void RememberView();
    void ClampView();

    std::string title_;
    std::string xLabel_;
    std::string yLabel_;
    std::string emptyMessage_ = "No spectrum loaded";
    std::vector<PlotSeries> series_;
    std::vector<PlotMarker> markers_;
    InteractionMode mode_ = InteractionMode::SelectRange;
    std::function<void(double)> rangeSelected_;
    double fullMinimum_ = 0.0;
    double fullMaximum_ = 1.0;
    double viewMinimum_ = 0.0;
    double viewMaximum_ = 1.0;
    std::vector<std::pair<double, double>> viewHistory_;
    QPoint dragStart_;
    QPoint dragCurrent_;
    bool dragging_ = false;
    Qt::MouseButton dragButton_ = Qt::NoButton;
};

} // namespace hpge
