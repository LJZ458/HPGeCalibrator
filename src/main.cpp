#include "MainWindow.h"

#include <TApplication.h>
#include <TEnv.h>
#include <TGClient.h>
#include <TGPicture.h>
#include <TROOT.h>
#include <TString.h>
#include <TVirtualX.h>

#if defined(HPGE_USE_COCOA)
#include <TGCocoa.h>
#elif defined(HPGE_USE_X11)
#include <TGX11.h>
#endif

#include <memory>

namespace {

class NativeClient final : public TGClient {
public:
    NativeClient(const char* displayName, const char* iconPath) : TGClient(displayName) {
        defaultPicturePool_ = fPicturePool;
        fPicturePool = new TGPicturePool(this, iconPath);
    }

    ~NativeClient() override {
        delete fPicturePool;
        fPicturePool = defaultPicturePool_;
    }

private:
    TGPicturePool* defaultPicturePool_ = nullptr;
};

} // namespace

int main(int argc, char** argv) {
    // Set this before TApplication: TROOT caches the icon search path during
    // application construction.
    gEnv->SetValue("Gui.IconPath", HPGE_ROOT_ICON_PATH);
    const_cast<TString&>(TROOT::GetIconPath()) = HPGE_ROOT_ICON_PATH;
    TApplication application("hpge-calibrator", &argc, argv);

    // Construct the platform graphics backend directly. ROOT normally loads this
    // through Cling's plugin manager; direct construction is deterministic and
    // also works when the installed compiler and ROOT's Cling differ.
    std::unique_ptr<TVirtualX> graphics;
#if defined(HPGE_USE_COCOA)
    graphics = std::make_unique<TGCocoa>("Cocoa", "HPGe Calibrator Cocoa backend");
#elif defined(HPGE_USE_X11)
    graphics = std::make_unique<TGX11>("X11", "HPGe Calibrator X11 backend");
#else
    application.InitializeGraphics();
#endif
    if (graphics) gVirtualX = graphics.get();
    std::unique_ptr<TGClient> client;
    if (!gClient) client = std::make_unique<NativeClient>(nullptr, HPGE_ROOT_ICON_PATH);
    if (!gClient) return 2;

    new hpge::MainWindow(gClient->GetRoot(), 1440, 900);
    application.Run();
    return 0;
}
