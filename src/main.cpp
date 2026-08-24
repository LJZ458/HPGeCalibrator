#include "MainWindow.h"

#include <TApplication.h>
#include <TGClient.h>
#include <TROOT.h>
#include <TVirtualX.h>

#include <cstdlib>
#include <exception>
#include <iostream>

namespace {

bool HasDisplay() {
#if defined(__linux__) || defined(__FreeBSD__)
    const char* display = std::getenv("DISPLAY");
    return display && display[0] != '\0';
#else
    return true;
#endif
}

void PrintHeadlessError() {
    std::cerr
        << "HPGe Calibrator requires a graphical display, but DISPLAY is not set.\n"
        << "Run it from a desktop session, use 'ssh -X'/'ssh -Y', or start it under Xvfb.\n";
}

} // namespace

int main(int argc, char** argv) {
    // ROOT's classic GUI on Linux uses X11 (directly or through XWayland). Avoid
    // entering ROOT's GUI loader when no display is available: some ROOT builds
    // abort or segfault instead of returning a null backend in this situation.
    if (!HasDisplay()) {
        PrintHeadlessError();
        return 2;
    }

    // This is ROOT's supported sequence for a compiled native GUI. It lets ROOT
    // select the correct platform plugin (GX11 on Linux, Cocoa on macOS) instead
    // of replacing the global TVirtualX/TGClient objects ourselves.
    TApplication::NeedGraphicsLibs();
    TApplication application("hpge-calibrator", &argc, argv);
    application.InitializeGraphics();

    if (gROOT->IsBatch() || !gVirtualX || gVirtualX == gGXBatch) {
        std::cerr
            << "HPGe Calibrator could not initialize ROOT's graphical backend.\n"
            << "Install a ROOT build with GUI/X11 support and verify that DISPLAY is reachable.\n";
        return 3;
    }

    TGClient* client = gClient;
    if (!client || !client->GetRoot()) {
        std::cerr
            << "HPGe Calibrator could not initialize ROOT's GUI client.\n"
            << "Check that libGui and the platform graphics plugin are installed.\n";
        return 4;
    }

    try {
        new hpge::MainWindow(client->GetRoot(), 1440, 900);
        application.Run();
    } catch (const std::exception& error) {
        std::cerr << "HPGe Calibrator failed during GUI startup: " << error.what() << '\n';
        return 5;
    } catch (...) {
        std::cerr << "HPGe Calibrator failed during GUI startup with an unknown error.\n";
        return 5;
    }
    return 0;
}
