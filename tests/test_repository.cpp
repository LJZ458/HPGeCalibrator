#include "RootDataRepository.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>

namespace {

class TestRootFile {
public:
    explicit TestRootFile(const std::string& suffix = "functions",
                          double chargeX = 42.0, double chargeY = 68.0) {
        path_ = std::filesystem::temp_directory_path() /
                ("hpge_repository_" + suffix + ".root");
        TFile file(path_.string().c_str(), "RECREATE");
        auto* nested = file.mkdir("nested");
        nested->cd();
        TH2D chargeOnX("charge_on_x", "charge on x", 100, 0.0, 100.0, 64, 0.0, 64.0);
        chargeOnX.Fill(chargeX, 7.5, 12.0);
        chargeOnX.Write();

        TH2D chargeOnY("charge_on_y", "charge on y", 64, 0.0, 64.0, 100, 0.0, 100.0);
        chargeOnY.Fill(9.5, chargeY, 15.0);
        chargeOnY.Write();
    }

    ~TestRootFile() { std::filesystem::remove(path_); }
    const std::filesystem::path& Path() const { return path_; }

private:
    std::filesystem::path path_;
};

const hpge::HistogramDescriptor* Find(const std::vector<hpge::HistogramDescriptor>& items,
                                      const std::string& name) {
    for (const auto& item : items) {
        if (item.objectPath == "nested/" + name) return &item;
    }
    return nullptr;
}

bool TestDiscoverAndLoad() {
    TestRootFile testFile("discover");
    hpge::RootDataRepository repository;
    std::string error;
    const auto discovered = repository.Discover(testFile.Path().string(), error);
    if (discovered.size() != 2 || !error.empty()) {
        std::cerr << "Expected two recursively discovered TH2 objects: " << error << '\n';
        return false;
    }
    const auto* descriptor = Find(discovered, "charge_on_x");
    if (!descriptor || descriptor->xBins != 100 || descriptor->yBins != 64 ||
        descriptor->displayName.find("charge_on_x") == std::string::npos) {
        std::cerr << "Discovered descriptor metadata is incorrect\n";
        return false;
    }
    const auto first = repository.Load(*descriptor, error);
    const auto cached = repository.Load(*descriptor, error);
    if (!first || first.get() != cached.get() || !error.empty()) {
        std::cerr << "Load/cache behavior is incorrect: " << error << '\n';
        return false;
    }
    repository.ClearCache();
    const auto reloaded = repository.Load(*descriptor, error);
    if (!reloaded || reloaded.get() == first.get()) {
        std::cerr << "ClearCache did not force a reload\n";
        return false;
    }
    return true;
}

bool TestProjectChargeOnX() {
    TestRootFile testFile("project_x");
    hpge::RootDataRepository repository;
    std::string error;
    const auto discovered = repository.Discover(testFile.Path().string(), error);
    const auto* descriptor = Find(discovered, "charge_on_x");
    if (!descriptor) return false;
    const auto projection = repository.ProjectCrystal(*descriptor, 7,
        hpge::AxisOrientation::ChargeOnX, error);
    if (!projection || !error.empty() || projection->Integral() < 11.9) {
        std::cerr << "Charge-on-X projection failed: " << error << '\n';
        return false;
    }
    return std::abs(projection->GetBinCenter(projection->GetMaximumBin()) - 42.0) < 1.0;
}

bool TestProjectChargeOnY() {
    TestRootFile testFile("project_y");
    hpge::RootDataRepository repository;
    std::string error;
    const auto discovered = repository.Discover(testFile.Path().string(), error);
    const auto* descriptor = Find(discovered, "charge_on_y");
    if (!descriptor) return false;
    const auto projection = repository.ProjectCrystal(*descriptor, 9,
        hpge::AxisOrientation::ChargeOnY, error);
    if (!projection || !error.empty() || projection->Integral() < 14.9) {
        std::cerr << "Charge-on-Y projection failed: " << error << '\n';
        return false;
    }
    return std::abs(projection->GetBinCenter(projection->GetMaximumBin()) - 68.0) < 1.0;
}

bool TestErrors() {
    hpge::RootDataRepository repository;
    std::string error;
    const auto missing = repository.Discover("/tmp/definitely_missing_hpge_file.root", error);
    if (!missing.empty() || error.find("Could not open") == std::string::npos) {
        std::cerr << "Missing-file error was not returned\n";
        return false;
    }

    TestRootFile testFile("errors");
    const auto discovered = repository.Discover(testFile.Path().string(), error);
    const auto* descriptor = Find(discovered, "charge_on_x");
    if (!descriptor) return false;
    const auto invalid = repository.ProjectCrystal(*descriptor, 64,
        hpge::AxisOrientation::ChargeOnX, error);
    if (invalid || error.find("outside") == std::string::npos) {
        std::cerr << "Out-of-range crystal was not rejected\n";
        return false;
    }

    auto stale = *descriptor;
    stale.objectPath = "nested/not_present";
    stale.id += "::missing";
    const auto missingObject = repository.Load(stale, error);
    if (missingObject || error.find("no longer exists") == std::string::npos) {
        std::cerr << "Missing TH2 object was not rejected\n";
        return false;
    }
    return true;
}

bool TestMultipleFilesStability() {
    TestRootFile firstFile("multi_a", 42.0, 68.0);
    TestRootFile secondFile("multi_b", 57.0, 73.0);
    hpge::RootDataRepository repository;
    std::string error;
    const auto firstItems = repository.Discover(firstFile.Path().string(), error);
    if (!error.empty()) return false;
    const auto secondItems = repository.Discover(secondFile.Path().string(), error);
    if (!error.empty()) return false;
    const auto* first = Find(firstItems, "charge_on_x");
    const auto* second = Find(secondItems, "charge_on_x");
    if (!first || !second || first->id == second->id) {
        std::cerr << "Multiple files did not produce independent descriptors\n";
        return false;
    }

    for (int iteration = 0; iteration < 100; ++iteration) {
        const auto& descriptor = iteration % 2 == 0 ? *first : *second;
        const double expected = iteration % 2 == 0 ? 42.0 : 57.0;
        auto projection = repository.ProjectCrystal(
            descriptor, 7, hpge::AxisOrientation::ChargeOnX, error);
        if (!projection || !error.empty() ||
            std::abs(projection->GetBinCenter(projection->GetMaximumBin()) - expected) >= 1.0) {
            std::cerr << "Repeated multi-file projection failed at iteration " << iteration
                      << ": " << error << '\n';
            return false;
        }
        if (iteration % 11 == 10) repository.ClearCache();
    }
    return true;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Expected one test-case name\n";
        return 2;
    }
    const std::string test = argv[1];
    bool passed = false;
    if (test == "discover-and-load") passed = TestDiscoverAndLoad();
    else if (test == "project-charge-on-x") passed = TestProjectChargeOnX();
    else if (test == "project-charge-on-y") passed = TestProjectChargeOnY();
    else if (test == "repository-errors") passed = TestErrors();
    else if (test == "multiple-files-stability") passed = TestMultipleFilesStability();
    else {
        std::cerr << "Unknown test case: " << test << '\n';
        return 2;
    }
    if (passed) std::cout << "PASS: " << test << '\n';
    return passed ? 0 : 1;
}
