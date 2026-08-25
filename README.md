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
- Accepts arbitrary custom background or contaminant energies with optional labels. Custom lines can be session-only or saved in the platform's per-user application settings for reuse on later launches; saved lines are visibly marked and can be removed.
- Records reference peaks by clicking the lower and upper fit limits on the first crystal's spectrum.
- Fits each selected interval with a RadWare/GF3-style Gaussian, low-energy tail, smoothed step, and quadratic background model; calibration uses the fitted centroid and uncertainty.
- Draws each fitted peak curve and centroid in red, without covering the spectrum with a selected-range band.
- Lists every source histogram that contributed fitted peaks to a selected crystal result. Review can show one source with its stored red fit curves or a vertically stacked, normalized overlay of all fitted source spectra and curves.
- Finds peak patterns independently in the reference and target spectra, without using relative peak intensities. Missing peaks and additional contaminant peaks are tolerated.
- Selects the best full-spectrum alignment by minimizing a one-to-one ordered peak-pattern cost. The cost includes squared positional mismatch plus separate penalties for missing reference lines, extra target lines, extreme gain/offset, and excessive curvature, so a locally plausible subset cannot override the complete spectrum pattern.
- Determines calibration peaks directly from the returned affine or quadratic alignment; there is no second candidate finder or subset-matching fallback. The user sets an aligned reference-charge half-range, both exact boundaries are transformed into the target spectrum, and the RadWare fit returns the centroid only within that fixed interval—without automatic recentering or widening.
- Provides a continuous 0–100% alignment-peak sensitivity. Optional auto-tuning adjusts it independently for every reference and target crystal spectrum using that projection's count level, dynamic range, and isolated-spike content.
- Uses square-root-compressed peak discovery and spatially balanced candidate limits so intense low-energy X-rays cannot hide weaker gamma lines or consume the alignment pattern. When at least three higher-charge candidates exist, auto-tuning excludes the lowest 7% of the charge axis from alignment only; original spectra and peak fits are unchanged.
- Seeds alignment from peak spacings, minimizes the global assignment, then refines the winning monotonic second-order charge mapping for broad sources such as Co-56; two-line Co-60 alignment remains affine.
- Allows the alignment mapping to be selected explicitly as affine, second-order polynomial, or automatic. The second-order alignment is `target charge = a0 + a1*q + a2*q*q`; its seeding uses a broad provisional correspondence followed by tight iterative rematching and curvature regularization.
- Refreshes the alignment plot immediately when another TH2 histogram is selected, includes the ROOT dataset name in the title and status, and clears stale plot data if the newly selected histogram cannot be aligned.
- Previews the pattern mapping before calibration by overlaying normalized spectra in a common reference-charge coordinate.
- Exposes selectable preview coefficients and a per-crystal list of stored `a0`, `a1`, `a2`, minimized cost, matched-pattern count, and effective reference/target sensitivities. Any stored parameter row can be copied to the clipboard.
- Provides separate **Zoom / pan** and **Select peak-fit range** mouse modes so inspecting a spectrum cannot accidentally create a fit interval. Wheel zoom and right-drag pan remain available while selecting fit limits; Back, Zoom −, Zoom +, Reset, left-drag zoom, and double-click reset provide ROOT-like navigation. The zoom window stays fixed when modes change or fitted overlays are added.
- Uses a custom Qt plot widget, so ROOT object selection, class/editor panels, and canvas callbacks cannot interrupt spectrum interaction.
- Saves and restores complete versioned `.hpgecal.json` projects from the **File** menu. A project records ROOT files with portable relative-path fallbacks, selected TH2 datasets/crystals, fitting controls, custom/reference peaks, pending manual peaks, applied manual and automatic peak fits, calibration coefficients/results, and alignment parameters. Restored projects remain editable and accept additional peaks and ROOT files.
- Combines peak points from multiple source histograms into one quadratic fit per crystal.
- Displays calibration curves and per-line residuals versus assigned peak energy, and flags high-RMS or exactly determined fits for review.
- Displays every fitted energy residual versus detector crystal on a fixed 0–63 axis, filterable by source histogram, with a second plot of per-crystal residual RMS and maximum absolute residual.
- Lists every fitted dataset/energy/centroid for the selected crystal. Selecting a point prepares only that corresponding peak for replacement; selecting another known/custom energy adds a new point.
- Applies pending replacements and additions to the existing point set without rerunning automatic peak finding or refitting untouched RadWare peaks. Only the second-order calibration polynomial and residuals are recomputed.
- Binds manual corrections to the single crystal highlighted in the calibration-results list and labels the correction panel with that crystal number.
- Extends calibration and residual plot axes beyond the outermost data points so the second-order trend and edge residuals have visual context.
- During recalibration, asks whether to refit every peak or preserve all stored centroids/fit curves and fit only newly added dataset-energy lines.
- After the individual fits, provides both per-source sums and one all-selected-sources energy spectrum containing every source/crystal spectrum and every assigned line.
- Lets the user select a combined-spectrum line, click its two fit limits, and refit only that peak while preserving the other combined fits; manual combined fits survive project save/restore.
- Exports coefficients, fit statistics, calibration points, peak-fit parameters, manual/automatic status, combined-spectrum residual/FWHM results, and resolution to CSV.
- Writes a companion C++ header containing three 64-element coefficient arrays—`p0`, `p1`, and `p2`—indexed directly by crystal number.

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

The test suite reports thirty-five focused checks separately: peak discovery and refinement,
single- and multi-peak mapping, RadWare peak fitting and validation, mapped interval fitting,
alignment-constrained peak identification that rejects a locally convincing contaminant pattern,
intensity-independent and two-line pattern alignment, wide-range Co-56 crystal alignment,
low-statistics spike rejection, continuous per-spectrum sensitivity tuning, quadratic alignment
in the presence of dominant low-energy X-rays, multi-histogram alignment-preview refresh,
multi-source fitted-spectrum result review, selective single-crystal peak replacement/addition,
preservation of untouched fitted peaks, persistent/session-only custom-line restoration and removal,
exact direct-alignment fit-range mapping and alignment-parameter access,
complete project save/restore/extension including manual fits, 64-crystal energy-residual
overview and per-crystal quality plots, and
alignment validation, quadratic fitting and validation, recursive histogram discovery/cache loading, both TH2 axis
orientations, repeated multi-file projection/cache cycling, repository error handling, and sample
ROOT-file generation, combined calibrated-spectrum/FWHM analysis, three-list C++ export,
Qt plot interaction/rendering, offscreen GUI startup, and a full
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

1. In **Data**, add one or several ROOT files in the file browser. Select one or more discovered `TH2` histograms, choose the axis orientation, reference crystal, and crystals to calibrate. To resume earlier work, use **File → Open project…**; `.hpgecal.json` projects reload their ROOT files and all stored fit state. They can also be passed to `hpge-calibrator` on the command line.
2. In **Reference peaks**, choose the first source histogram and show its reference spectrum. Select the radioactive source first so only that source's energies are listed. For a custom line, enter its energy and optional label, then choose **Add for session** or **Save for future**; persistent entries appear under **Custom** with a `[saved]` marker and can be deleted with **Remove selected custom line**. Use the zoom toolbar, wheel, or **Zoom / pan** mode to inspect the spectrum, then switch to **Select peak-fit range**, select an energy, and click the lower and upper fit limits. The zoomed view remains in place, and the red RadWare fit curve, fitted centroid, and energy label remain visible. Repeat for all usable lines and source histograms.
3. In **Calibration & review**, choose the starting peak sensitivity, **Aligned peak-fit half-range**, charge-mapping model, histogram, and target crystal under **Pre-calibration spectrum alignment**, then click **Show aligned spectra**. The half-range is defined in the displayed reference-charge coordinate; automatic calibration maps both boundaries exactly into each target spectrum and fits only inside that interval. Leave **Auto-tune for each crystal spectrum** enabled to adjust the sensitivity independently for each reference/target projection; disable it when you want the exact 0–100% value. Use **Auto — affine or 2nd order** normally, or force `a0 + a1*q` or the second-order polynomial `a0 + a1*q + a2*q*q`. Selecting another histogram refreshes the preview immediately; its dataset name appears in both the plot title and status. The blue reference and red target spectra are normalized and overlaid after independent peak-pattern mapping; no user-selected energy lines or peak intensities are used to obtain the alignment. The selectable preview text reports the coefficients and quality values. After calibration, select a result to inspect or copy every dataset's stored alignment row. No energy calibration is applied in the preview.
4. Adjust peak-search parameters if necessary and calibrate the selected crystals. At least three total reference points are required for a quadratic fit; four or more provide a meaningful residual-based quality check.
5. Select a calibration result. **Show all fitted sources** displays every contributing source spectrum as a normalized vertical stack with all fitted peak curves in red. Use the fitted-source selector or **Show selected source** to inspect one spectrum at full count scale.
6. For any `REVIEW` or `FAIL` result, inspect **Existing fitted peaks** under Manual correction. Selecting a row chooses its dataset and energy and opens that crystal spectrum for replacement. Select a new range to create a pending replacement, or choose another known/custom energy and range to add a point. Remove any unwanted pending item, then press **Apply pending peaks (keep all others)**. Untouched centroids and RadWare fits are preserved; only the calibration polynomial and residuals are recomputed.
7. Inspect **Fit + residuals**; residual x coordinates are the assigned peak energies rather than charge centroids.
8. Under **Combined calibrated spectrum quality**, choose **All selected source spectra** to inspect every selected histogram in one energy-space sum, or choose an individual source for diagnosis. Select a listed line and click **Refit selected combined peak**, then click the two desired energy limits.
9. When running calibration again, choose **Keep fitted peaks + add new** to retain every prior peak fit and process only newly assigned lines, or **Refit all peaks** to replace the automatic peak fits.
9. Under **Residual overview across detector crystals**, choose all sources or one source and click **Show residuals vs crystals**. The upper plot overlays every energy-line residual on crystal indices 0–63; the lower plot shows each crystal's RMS and largest absolute residual.
10. Use **File → Save project…** whenever you want a resumable snapshot. After restoring it, existing manual and automatic fits remain available for review and selective correction, and **Add ROOT files…** plus the normal reference-peak controls can extend the analysis before saving again.
11. Export the complete CSV. The same action also writes `<csv-name>_coefficients.hpp` with separate `p0`, `p1`, and `p2` C++ arrays for crystals 0–63; unfitted crystals contain a `missing` NaN placeholder.

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

CSV output contains one `fit` row per crystal followed by one `point` row per calibration line. Fit rows include `p0`, `p1`, `p2`, chi-square, degrees of freedom, residual RMS, and review status. Point rows include dataset identity, fitted centroid and uncertainty, selected range, peak sigma, RadWare fit statistics and tail/step parameters, assigned energy, residual, and whether the point was manually overridden. `combined_peak` rows contain the expected and refitted energies, energy residual, FWHM in keV, percentage resolution, and number of contributing crystals. A companion C++ header contains exactly three coefficient arrays named `p0`, `p1`, and `p2`.

Combined-peak FWHM is reported for the fitted Gaussian core as `2.354820045 × sigma`; percentage resolution is `100 × FWHM / expected energy`.

The peak-shape implementation follows the public RadWare/GF3 model used for gamma-spectrum analysis: a Gaussian photopeak plus a low-energy exponential tail, a smoothed step, and polynomial background. See the [RadWare source repository](https://github.com/radforddc/rw05) and the [GammaSpecAnalysis RadWare fit documentation](https://nucleardata.berkeley.edu/nsd_software/doxy/html/d0/d67/_spectrum_analysis_8cpp.html).

An exactly three-point quadratic has zero degrees of freedom and is marked `REVIEW`, even when its residuals are numerically zero. This is deliberate: such a fit has no independent information with which to assess quality.
