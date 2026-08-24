cmake_minimum_required(VERSION 3.18)
project(HPGeCalibrator VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(ROOT 6.26 REQUIRED COMPONENTS Core RIO Hist Spectrum)

# Debian packages expose Qt through the normal CMake prefix. Homebrew keeps it
# keg-only, so discover that prefix automatically when no Qt hint was supplied.
find_package(Qt6 6.4 QUIET COMPONENTS Widgets)
if(APPLE AND NOT Qt6_FOUND AND NOT DEFINED Qt6_DIR)
    find_program(HOMEBREW_EXECUTABLE brew)
    if(HOMEBREW_EXECUTABLE)
        execute_process(COMMAND ${HOMEBREW_EXECUTABLE} --prefix qt
                        OUTPUT_VARIABLE HOMEBREW_QT_PREFIX
                        OUTPUT_STRIP_TRAILING_WHITESPACE
                        ERROR_QUIET)
        if(HOMEBREW_QT_PREFIX)
            set(Qt6_DIR "${HOMEBREW_QT_PREFIX}/lib/cmake/Qt6")
        endif()
    endif()
endif()
find_package(Qt6 6.4 REQUIRED COMPONENTS Widgets)

add_library(hpge_core
    src/CalibrationExporter.cpp
    src/CalibrationEngine.cpp
    src/CombinedSpectrumAnalyzer.cpp
    src/RootDataRepository.cpp
)
target_include_directories(hpge_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(hpge_core PUBLIC ROOT::Core ROOT::RIO ROOT::Hist ROOT::Spectrum)

add_executable(hpge-calibrator
    src/main.cpp
    src/MainWindow.cpp
    src/SpectrumPlotWidget.cpp
)
target_include_directories(hpge-calibrator PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(hpge-calibrator PRIVATE
    hpge_core Qt6::Widgets
)

add_executable(hpge-generate-sample tools/generate_sample.cpp)
target_link_libraries(hpge-generate-sample PRIVATE ROOT::Core ROOT::RIO ROOT::Hist)

include(CTest)
if(BUILD_TESTING)
    add_executable(hpge-calibration-tests tests/test_calibration.cpp)
    target_link_libraries(hpge-calibration-tests PRIVATE hpge_core)
    foreach(test_case IN ITEMS
            find-peaks
            refine-peak
            match-multiple-peaks
            match-single-peak
            radware-peak-fit
            radware-fit-validation
            mapped-radware-fit
            spectrum-alignment
            two-peak-pattern-alignment
            co56-crystal-pattern-alignment
            low-statistics-spike-alignment
            adaptive-spectrum-sensitivity
            quadratic-xray-alignment
            pattern-alignment-validation
            fit-quadratic
            fit-validation)
        add_test(NAME calibration-${test_case}
                 COMMAND hpge-calibration-tests ${test_case})
    endforeach()

    add_executable(hpge-repository-tests tests/test_repository.cpp)
    target_link_libraries(hpge-repository-tests PRIVATE hpge_core)
    foreach(test_case IN ITEMS
            discover-and-load
            project-charge-on-x
            project-charge-on-y
            multiple-files-stability
            repository-errors)
        add_test(NAME repository-${test_case}
                 COMMAND hpge-repository-tests ${test_case})
    endforeach()

    add_executable(hpge-combined-spectrum-tests tests/test_combined_spectrum.cpp)
    target_link_libraries(hpge-combined-spectrum-tests PRIVATE hpge_core)
    add_test(NAME calibration-combined-spectrum-quality COMMAND hpge-combined-spectrum-tests)

    add_executable(hpge-export-tests tests/test_export.cpp)
    target_link_libraries(hpge-export-tests PRIVATE hpge_core)
    add_test(NAME calibration-cpp-coefficient-export COMMAND hpge-export-tests)

    add_test(NAME sample-generator
             COMMAND hpge-generate-sample ${CMAKE_CURRENT_BINARY_DIR}/hpge_test_sample.root)

    add_executable(hpge-plot-tests
        tests/test_plot_widget.cpp
        src/SpectrumPlotWidget.cpp
    )
    target_include_directories(hpge-plot-tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(hpge-plot-tests PRIVATE Qt6::Widgets)
    add_test(NAME qt-plot-interaction COMMAND hpge-plot-tests)
    set_tests_properties(qt-plot-interaction PROPERTIES ENVIRONMENT "QT_QPA_PLATFORM=offscreen")

    add_test(NAME qt-gui-startup COMMAND hpge-calibrator --check-startup)
    set_tests_properties(qt-gui-startup PROPERTIES ENVIRONMENT "QT_QPA_PLATFORM=offscreen")

    add_test(NAME qt-gui-root-file-workflow
             COMMAND hpge-calibrator ${CMAKE_CURRENT_BINARY_DIR}/hpge_test_sample.root --check-startup)
    set_tests_properties(qt-gui-root-file-workflow PROPERTIES
        ENVIRONMENT "QT_QPA_PLATFORM=offscreen"
        DEPENDS sample-generator)

    add_executable(hpge-gui-alignment-tests
        tests/test_gui_alignment.cpp
        src/MainWindow.cpp
        src/SpectrumPlotWidget.cpp
    )
    target_include_directories(hpge-gui-alignment-tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(hpge-gui-alignment-tests PRIVATE hpge_core Qt6::Widgets)
    add_test(NAME qt-alignment-multiple-histogram-preview
             COMMAND hpge-gui-alignment-tests ${CMAKE_CURRENT_BINARY_DIR}/hpge_test_sample.root)
    set_tests_properties(qt-alignment-multiple-histogram-preview PROPERTIES
        ENVIRONMENT "QT_QPA_PLATFORM=offscreen"
        DEPENDS sample-generator)
endif()

install(TARGETS hpge-calibrator hpge-generate-sample RUNTIME DESTINATION bin)
