#include "RootDataRepository.h"

#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2.h>
#include <TKey.h>
#include <TROOT.h>

#include <filesystem>
#include <sstream>

namespace hpge {
namespace {

void DiscoverDirectory(TDirectory& directory, const std::string& prefix,
                       const std::string& filePath,
                       std::vector<HistogramDescriptor>& output) {
    TIter next(directory.GetListOfKeys());
    while (auto* key = dynamic_cast<TKey*>(next())) {
        const std::string name = key->GetName();
        const std::string objectPath = prefix.empty() ? name : prefix + "/" + name;
        // Ask the directory directly instead of resolving every key through
        // TClass. This keeps browsing on ROOT's compiled I/O/dictionary path
        // and avoids unnecessarily starting Cling for each object.
        if (auto* child = directory.GetDirectory(name.c_str())) {
            DiscoverDirectory(*child, objectPath, filePath, output);
            continue;
        }
        auto* histogram = dynamic_cast<TH2*>(directory.Get(name.c_str()));
        if (!histogram) continue;
        HistogramDescriptor descriptor;
        descriptor.id = filePath + "::" + objectPath;
        descriptor.filePath = filePath;
        descriptor.objectPath = objectPath;
        descriptor.displayName = std::filesystem::path(filePath).filename().string() +
                                 " :: " + objectPath;
        descriptor.xBins = histogram->GetNbinsX();
        descriptor.yBins = histogram->GetNbinsY();
        output.push_back(std::move(descriptor));
    }
}

} // namespace

std::vector<HistogramDescriptor> RootDataRepository::Discover(
    const std::string& filePath, std::string& error) const {
    error.clear();
    std::unique_ptr<TFile> file(TFile::Open(filePath.c_str(), "READ"));
    if (!file || file->IsZombie()) {
        error = "Could not open ROOT file: " + filePath;
        return {};
    }
    std::vector<HistogramDescriptor> output;
    DiscoverDirectory(*file, "", filePath, output);
    if (output.empty()) error = "No TH2 histograms found in " + filePath;
    return output;
}

std::shared_ptr<TH2> RootDataRepository::Load(const HistogramDescriptor& descriptor,
                                             std::string& error) {
    error.clear();
    if (const auto cached = cache_.find(descriptor.id); cached != cache_.end()) {
        return cached->second;
    }
    std::unique_ptr<TFile> file(TFile::Open(descriptor.filePath.c_str(), "READ"));
    if (!file || file->IsZombie()) {
        error = "Could not open ROOT file: " + descriptor.filePath;
        return {};
    }
    auto* source = dynamic_cast<TH2*>(file->Get(descriptor.objectPath.c_str()));
    if (!source) {
        error = "TH2 no longer exists: " + descriptor.objectPath;
        return {};
    }
    const std::string cloneName = "hpge2_" + std::to_string(std::hash<std::string>{}(descriptor.id));
    auto* clone = dynamic_cast<TH2*>(source->Clone(cloneName.c_str()));
    if (!clone) {
        error = "Could not clone TH2: " + descriptor.objectPath;
        return {};
    }
    clone->SetDirectory(nullptr);
    auto result = std::shared_ptr<TH2>(clone);
    cache_[descriptor.id] = result;
    return result;
}

std::shared_ptr<TH1D> RootDataRepository::ProjectCrystal(
    const HistogramDescriptor& descriptor, int crystal, AxisOrientation orientation,
    std::string& error) {
    auto histogram = Load(descriptor, error);
    if (!histogram) return {};
    const int channelBins = orientation == AxisOrientation::ChargeOnX
                                ? histogram->GetNbinsY()
                                : histogram->GetNbinsX();
    if (crystal < 0 || crystal >= channelBins) {
        std::ostringstream message;
        message << "Crystal " << crystal << " is outside the histogram channel axis ("
                << channelBins << " bins)";
        error = message.str();
        return {};
    }
    const std::string name = "spectrum_" +
                             std::to_string(std::hash<std::string>{}(descriptor.id)) + "_" +
                             std::to_string(crystal) + "_" +
                             (orientation == AxisOrientation::ChargeOnX ? "x" : "y");
    TH1D* projection = nullptr;
    if (orientation == AxisOrientation::ChargeOnX) {
        projection = histogram->ProjectionX(name.c_str(), crystal + 1, crystal + 1, "e");
    } else {
        projection = histogram->ProjectionY(name.c_str(), crystal + 1, crystal + 1, "e");
    }
    if (!projection) {
        error = "ROOT failed to project crystal spectrum";
        return {};
    }
    projection->SetDirectory(nullptr);
    projection->SetTitle((descriptor.displayName + " | crystal " +
                          std::to_string(crystal) + ";Charge;Counts")
                             .c_str());
    return std::shared_ptr<TH1D>(projection);
}

void RootDataRepository::ClearCache() { cache_.clear(); }

} // namespace hpge
