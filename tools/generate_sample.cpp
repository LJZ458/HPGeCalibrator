#include <TDirectory.h>
#include <TFile.h>
#include <TH2D.h>
#include <TRandom3.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

double ChargeForEnergy(double energy, double p0, double p1, double p2) {
    if (std::abs(p2) < 1e-15) return (energy - p0) / p1;
    return (-p1 + std::sqrt(p1 * p1 - 4.0 * p2 * (p0 - energy))) / (2.0 * p2);
}

void FillDataset(TH2D& histogram, const std::vector<double>& energies, TRandom3& random) {
    for (int crystal = 0; crystal < 64; ++crystal) {
        const double p0 = -1.5 + 0.04 * crystal;
        const double p1 = 0.68 + 0.0015 * crystal;
        const double p2 = 1.8e-5 + 1.0e-7 * (crystal % 7);
        for (int background = 0; background < 6000; ++background) {
            histogram.Fill(random.Exp(900.0), crystal + 0.5);
        }
        for (double energy : energies) {
            const double center = ChargeForEnergy(energy, p0, p1, p2);
            const double sigma = 2.2 + 0.0007 * center;
            for (int event = 0; event < 2500; ++event) {
                histogram.Fill(random.Gaus(center, sigma), crystal + 0.5);
            }
        }
    }
}

} // namespace

int main(int argc, char** argv) {
    const std::string outputPath = argc > 1 ? argv[1] : "hpge_sample.root";
    TFile output(outputPath.c_str(), "RECREATE");
    if (output.IsZombie()) {
        std::cerr << "Could not create " << outputPath << '\n';
        return 1;
    }
    TRandom3 random(20260824);
    auto* sources = output.mkdir("sources");
    sources->cd();
    TH2D co60("co60_charge_vs_crystal",
              "Co-60 charge versus crystal;Charge (ADC);Crystal index",
              4096, 0.0, 4096.0, 64, 0.0, 64.0);
    FillDataset(co60, {1173.228, 1332.492}, random);
    co60.Write();

    TH2D co56("co56_charge_vs_crystal",
              "Co-56 charge versus crystal;Charge (ADC);Crystal index",
              4096, 0.0, 4096.0, 64, 0.0, 64.0);
    FillDataset(co56, {846.771, 1238.282, 1771.351, 2598.459}, random);
    co56.Write();

    output.Write();
    output.Close();
    std::cout << "Wrote " << outputPath << " with two TH2 datasets and 64 crystals.\n";
    return 0;
}

