cmake_minimum_required(VERSION 3.18)
project(HPGeCalibrator VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(ROOT 6.26 REQUIRED COMPONENTS Core RIO Hist Graf Gpad Gui Spectrum Tree)
include(${ROOT_USE_FILE})

add_library(hpge_core
    src/CalibrationEngine.cpp
    src/RootDataRepository.cpp
)
target_include_directories(hpge_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(hpge_core PUBLIC ROOT::Core ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Spectrum)

add_executable(hpge-calibrator
    src/main.cpp
    src/MainWindow.cpp
)
target_include_directories(hpge-calibrator PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
list(GET ROOT_INCLUDE_DIRS 0 HPGE_ROOT_INCLUDE_DIR)
get_filename_component(HPGE_ROOT_PREFIX "${HPGE_ROOT_INCLUDE_DIR}/../.." ABSOLUTE)
target_compile_definitions(hpge-calibrator PRIVATE
    HPGE_ROOT_ICON_PATH="${HPGE_ROOT_PREFIX}/share/root/icons"
)
target_link_libraries(hpge-calibrator PRIVATE
    hpge_core ROOT::Gui ROOT::Gpad ROOT::Tree
)
if(APPLE AND TARGET ROOT::GCocoa)
    target_link_libraries(hpge-calibrator PRIVATE ROOT::GCocoa)
    target_compile_definitions(hpge-calibrator PRIVATE HPGE_USE_COCOA=1)
elseif(UNIX AND TARGET ROOT::GX11)
    target_link_libraries(hpge-calibrator PRIVATE ROOT::GX11)
    target_compile_definitions(hpge-calibrator PRIVATE HPGE_USE_X11=1)
endif()

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
            repository-errors)
        add_test(NAME repository-${test_case}
                 COMMAND hpge-repository-tests ${test_case})
    endforeach()
endif()

install(TARGETS hpge-calibrator hpge-generate-sample RUNTIME DESTINATION bin)
