# HPGe Crystal Calibrator

HPGe Crystal Calibrator is a Qt 6 desktop application for energy-calibrating up to 64 HPGe crystals from charge-versus-crystal ROOT `TH2` histograms. CERN ROOT is used only for file I/O, histograms, peak search, and analysis; the event loop, file browser, controls, and plots do not use ROOT GUI classes. It uses reference-spectrum peak selection, automatic spectrum-to-spectrum peak mapping, and a second-order calibration

\[
E(q) = p_0 + p_1 q + p_2 q^2.
\]

The application runs on macOS and Linux with Qt 6 and a CERN ROOT installation. ROOT itself does not need GUI, X11, Gpad, or Graf support.

## Features

- Uses the native Qt multi-file browser, recursively discovers every `TH2`, and immediately previews the first spectrum.
- Selects any combination of histograms and any subset of crystals 0–63.
- Supports both common axis layouts: charge-on-X/crystal-on-Y and the inverse.
- Separates commonly used Co-60, Co-56, Cs-137, Na-22, background, contaminant, and custom lines by source.
- Accepts arbitrary custom background or contaminant energies.
- Records reference peaks by clicking the lower and upper fit limits on the first crystal's spectrum.
- Fits each selected interval with a RadWare/GF3-style Gaussian, low-energy tail, smoothed step, and quadratic background model; calibration uses the fitted centroid and uncertainty.
- Draws each fitted peak curve and centroid in red, without covering the spectrum with a selected-range band.
- Finds peak patterns independently in the reference and target spectra, without using the user-assigned energy lines or relative peak intensities. Missing peaks and additional contaminant peaks are tolerated.
- Seeds alignment from peak spacings, then refines a monotonic second-order charge mapping for broad sources such as Co-56; two-line Co-60 alignment remains affine.
- Previews the pattern mapping before calibration by overlaying normalized spectra in a common reference-charge coordinate.
- Provides separate **Zoom / pan** and **Select peak-fit range** mouse modes so inspecting a spectrum cannot accidentally create a fit interval. Wheel and left-drag zoom, right-drag pans, and double-click resets the view.
- Uses a custom Qt plot widget, so ROOT object selection, class/editor panels, and canvas callbacks cannot interrupt spectrum interaction.
- Combines peak points from multiple source histograms into one quadratic fit per crystal.
- Displays calibration curves and per-line residuals, and flags high-RMS or exactly determined fits for review.
- Lets the user replace any automatically matched point by selecting a new range in a problematic crystal and refitting it.
- Exports coefficients, fit statistics, calibration points, peak-fit parameters, manual/automatic status, and residuals to CSV.

## Prerequisites

- CMake 3.18 or newer
- A C++17 compiler
- Qt 6.4 or newer with Widgets
- CERN ROOT 6.26 or newer with Core, RIO, Hist, and Spectrum components

ROOT installed through Homebrew, conda-forge, Spack, or an official ROOT binary is supported. Make sure `root-config` is on `PATH`, or pass ROOT's CMake directory explicitly with `-DROOT_DIR=...`.

### Debian 12 and newer

Install the compiler, CMake, and Qt Widgets development package:

```bash
sudo apt update
sudo apt install build-essential cmake qt6-base-dev
```

Install ROOT 6.26+ using a CERN binary or conda environment, then activate/source that ROOT installation before configuring this project. A ROOT build without its classic GUI is sufficient.

### macOS with Homebrew

```bash
brew install cmake qt root
```

CMake automatically checks Homebrew's keg-only Qt prefix if normal Qt discovery does not find it.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

The test suite reports twenty-two focused checks separately: peak discovery and refinement,
single- and multi-peak mapping, RadWare peak fitting and validation, mapped interval fitting,
intensity-independent and two-line pattern alignment, wide-range Co-56 crystal alignment and
alignment validation, quadratic fitting and validation, recursive histogram discovery/cache loading, both TH2 axis
orientations, repeated multi-file projection/cache cycling, repository error handling, and sample
ROOT-file generation, plus Qt plot interaction/rendering, offscreen GUI startup, and a full
GUI open/discover/project/render check using the generated ROOT file.

Run the application:

```bash
./build/hpge-calibrator
```

You can also pass one or more ROOT files at startup; each file is added to the browser before the window is shown:

```bash
./build/hpge-calibrator run60Co.root run56Co.root
```

On Linux, Qt uses the desktop's Wayland or X11 platform plugin. It no longer initializes ROOT graphics or requires ROOT's X11 backend.

If ROOT prints Cling module-map errors on a newly upgraded macOS/Xcode installation, reinstall or upgrade ROOT so it is built against the active Xcode command-line tools. The application avoids runtime-interpreted calibration and callback code, but ROOT's own optional image/plugin loader can still warn when the ROOT and Xcode standard-library versions do not match.

### Debian/Linux startup troubleshooting

Confirm that Qt, ROOT, and their shared libraries resolve:

```bash
ldd ./build/hpge-calibrator | grep 'not found'
QT_DEBUG_PLUGINS=1 ./build/hpge-calibrator
```

If Qt reports that the `xcb` platform plugin cannot be initialized on a minimal Debian install, verify that `qt6-qpa-plugins` is installed (it is a dependency of Debian's `qt6-base-dev`). For automated headless startup checks, use `QT_QPA_PLATFORM=offscreen ./build/hpge-calibrator --check-startup`.

## Calibration workflow

1. In **Data**, add one or several ROOT files in the file browser. Select one or more discovered `TH2` histograms, choose the axis orientation, reference crystal, and crystals to calibrate.
2. In **Reference peaks**, choose the first source histogram and show its reference spectrum. Select the radioactive source first so only that source's energies are listed. Use **Zoom / pan** above the plot to inspect the spectrum, switch to **Select peak-fit range**, select an energy, then click the lower and upper fit limits. The red RadWare fit curve, fitted centroid, and energy label remain visible. Repeat for all usable lines and source histograms.
3. In **Calibration & review**, choose a histogram and target crystal under **Pre-calibration spectrum alignment**, then click **Show aligned spectra**. The blue reference and red target spectra are normalized and overlaid after independent peak-pattern mapping; no user-selected energy lines or peak intensities are used. The status bar reports the matched-pattern count and charge-mapping coefficients. This is only a charge-axis alignment preview—no energy calibration has been applied.
4. Adjust peak-search parameters if necessary and calibrate the selected crystals. At least three total reference points are required for a quadratic fit; four or more provide a meaningful residual-based quality check.
5. Select any `REVIEW` or `FAIL` result. Choose its source histogram and energy, show the spectrum, click the lower and upper limits around the correct peak, and press **Refit crystal**. A manual centroid replaces the automatic point for that dataset and energy. After refitting, the spectrum stays visible with the fitted peak overlays.
6. Inspect **Fit + residuals**, then export the complete result table to CSV.

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

CSV output contains one `fit` row per crystal followed by one `point` row per calibration line. Fit rows include `p0`, `p1`, `p2`, chi-square, degrees of freedom, residual RMS, and review status. Point rows include dataset identity, fitted centroid and uncertainty, selected range, peak sigma, RadWare fit statistics and tail/step parameters, assigned energy, residual, and whether the point was manually overridden.

The peak-shape implementation follows the public RadWare/GF3 model used for gamma-spectrum analysis: a Gaussian photopeak plus a low-energy exponential tail, a smoothed step, and polynomial background. See the [RadWare source repository](https://github.com/radforddc/rw05) and the [GammaSpecAnalysis RadWare fit documentation](https://nucleardata.berkeley.edu/nsd_software/doxy/html/d0/d67/_spectrum_analysis_8cpp.html).

An exactly three-point quadratic has zero degrees of freedom and is marked `REVIEW`, even when its residuals are numerically zero. This is deliberate: such a fit has no independent information with which to assess quality.
