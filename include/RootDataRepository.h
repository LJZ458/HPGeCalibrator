#pragma once

#include "CalibrationTypes.h"

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

class TH1D;
class TH2;

namespace hpge {

struct HistogramDescriptor {
    std::string id;
    std::string filePath;
    std::string objectPath;
    std::string displayName;
    int xBins = 0;
    int yBins = 0;
};

class RootDataRepository {
public:
    std::vector<HistogramDescriptor> Discover(const std::string& filePath,
                                              std::string& error) const;
    std::shared_ptr<TH2> Load(const HistogramDescriptor& descriptor,
                             std::string& error);
    std::shared_ptr<TH1D> ProjectCrystal(const HistogramDescriptor& descriptor,
                                        int crystal,
                                        AxisOrientation orientation,
                                        std::string& error);
    void ClearCache();

private:
    std::unordered_map<std::string, std::shared_ptr<TH2>> cache_;
};

} // namespace hpge

