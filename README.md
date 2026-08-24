# HPGe Crystal Calibrator

HPGe Crystal Calibrator is a native CERN ROOT desktop application for energy-calibrating up to 64 HPGe crystals from charge-versus-crystal `TH2` histograms. It uses reference-spectrum peak selection, automatic spectrum-to-spectrum peak mapping, and a second-order calibration

\[
E(q) = p_0 + p_1 q + p_2 q^2.
\]

The application runs on macOS and Linux anywhere CERN ROOT with GUI support is available.

## Features

- Recursively browses ROOT files and lists every `TH2`, including objects inside ROOT directories.
- Selects any combination of histograms and any subset of crystals 0–63.
- Supports both common axis layouts: charge-on-X/crystal-on-Y and the inverse.
- Provides commonly used Co-60, Co-56, Cs-137, Na-22, K-40, and Tl-208 lines.
- Accepts arbitrary custom background or contaminant energies.
- Records reference peaks interactively by clicking the first crystal's spectrum.
- Finds corresponding peaks in the other crystals with `TSpectrum` and an affine peak-pattern mapping, without requiring predefined charge windows.
- Combines peak points from multiple source histograms into one quadratic fit per crystal.
- Displays calibration curves and per-line residuals, and flags high-RMS or exactly determined fits for review.
- Lets the user replace any automatically matched point by clicking the correct peak in a problematic crystal and refitting it.
- Exports coefficients, fit statistics, calibration points, manual/automatic status, and residuals to CSV.

## Prerequisites

- CMake 3.18 or newer
- A C++17 compiler
- CERN ROOT 6.26 or newer, built with GUI and Spectrum components

ROOT installed through Homebrew, conda-forge, Spack, or the official ROOT packages is supported. Make sure `root-config` is on `PATH`, or pass ROOT's CMake directory explicitly with `-DROOT_DIR=...`.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

The test suite reports ten focused checks separately: peak discovery, peak refinement,
single- and multi-peak mapping, quadratic fitting, fit validation, recursive histogram
discovery/cache loading, both TH2 axis orientations, and repository error handling.

Run the application:

```bash
./build/hpge-calibrator
```

On macOS, start it from a normal graphical Terminal session. On Linux, an X11 or Wayland display accessible to ROOT is required.

If ROOT prints Cling module-map errors on a newly upgraded macOS/Xcode installation, reinstall or upgrade ROOT so it is built against the active Xcode command-line tools. The application avoids runtime-interpreted calibration and callback code, but ROOT's own optional image/plugin loader can still warn when the ROOT and Xcode standard-library versions do not match.

## Calibration workflow

1. In **Data**, add each ROOT file. Select one or more discovered `TH2` histograms, choose the axis orientation, reference crystal, and crystals to calibrate.
2. In **Reference peaks**, choose the first source histogram and show its reference spectrum. Select a known energy (or add a custom one), then click the corresponding peak. Repeat for all usable lines and all source histograms. Clicking near a peak snaps to its local maximum.
3. In **Calibration & review**, adjust peak-search parameters if necessary and calibrate the selected crystals. At least three total reference points are required for a quadratic fit; four or more provide a meaningful residual-based quality check.
4. Select any `REVIEW` or `FAIL` result. Choose its source histogram and energy, show the spectrum, click the correct peak, and press **Refit crystal**. A manual point replaces the automatic point for that dataset and energy.
5. Inspect **Fit + residuals**, then export the complete result table to CSV.

## Histogram convention

For the default orientation, X contains charge and Y contains crystal/channel. Crystal 0 is projected from Y bin 1, crystal 1 from Y bin 2, and so on. The application validates the number of bins on the selected channel axis. It intentionally uses bin order instead of axis labels so integer axes defined as either `[0,64]` or `[1,65]` work consistently.

Every selected source histogram is handled independently. Reference peak assignments belong to the histogram on which they were picked, so different files can contain different sources or background lines.

## Generate reproducible sample data

The build also creates a small utility that writes Co-60 and Co-56 `TH2` spectra for 64 crystals:

```bash
./build/hpge-generate-sample hpge_sample.root
```

Load both histograms under the `sources` directory. On crystal 0, select the two Co-60 lines in the Co-60 histogram and several Co-56 lines in the Co-56 histogram, then run the full calibration.

## Calibration output

CSV output contains one `fit` row per crystal followed by one `point` row per calibration line. Fit rows include `p0`, `p1`, `p2`, chi-square, degrees of freedom, residual RMS, and review status. Point rows include dataset identity, measured charge, assigned energy, residual, and whether the point was manually overridden.

An exactly three-point quadratic has zero degrees of freedom and is marked `REVIEW`, even when its residuals are numerically zero. This is deliberate: such a fit has no independent information with which to assess quality.
