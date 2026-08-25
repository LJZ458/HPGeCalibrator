#include "MainWindow.h"

#include <QApplication>
#include <QCoreApplication>
#include <QStyleFactory>

#include <exception>
#include <iostream>
#include <vector>

int main(int argc, char** argv) {
    QApplication application(argc, argv);
    QCoreApplication::setApplicationName("HPGe Calibrator");
    QCoreApplication::setOrganizationName("HPGeCalibrator");
    if (auto* fusion = QStyleFactory::create("Fusion")) application.setStyle(fusion);

    try {
        hpge::MainWindow window;
        std::vector<std::string> rootFiles;
        std::string projectFile;
        const QStringList arguments = QCoreApplication::arguments();
        for (int index = 1; index < arguments.size(); ++index) {
            if (arguments[index] == "--screenshot") {
                ++index;
                continue;
            }
            if (arguments[index].endsWith(".root", Qt::CaseInsensitive)) {
                rootFiles.push_back(arguments[index].toStdString());
            } else if (arguments[index].endsWith(".hpgecal.json", Qt::CaseInsensitive)) {
                projectFile = arguments[index].toStdString();
            }
        }
        if (!projectFile.empty()) {
            std::string error;
            if (!window.OpenProject(projectFile, error)) {
                std::cerr << error << '\n';
                return 4;
            }
        }
        if (!rootFiles.empty()) window.OpenRootFiles(rootFiles);
        window.show();
        const int screenshotArgument = QCoreApplication::arguments().indexOf("--screenshot");
        if (screenshotArgument >= 0 && screenshotArgument + 1 < QCoreApplication::arguments().size()) {
            application.processEvents();
            const bool saved = window.grab().save(QCoreApplication::arguments().at(screenshotArgument + 1));
            window.close();
            return saved ? 0 : 3;
        }
        if (QCoreApplication::arguments().contains("--check-startup")) {
            application.processEvents();
            window.close();
            return 0;
        }
        return application.exec();
    } catch (const std::exception& error) {
        std::cerr << "HPGe Calibrator failed during startup: " << error.what() << '\n';
        return 2;
    } catch (...) {
        std::cerr << "HPGe Calibrator failed during startup with an unknown error.\n";
        return 2;
    }
}
