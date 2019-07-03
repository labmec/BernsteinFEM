cmake_minimum_required(VERSION 3.0.0)

project(BernsteinFEM CXX)

find_package(PZ REQUIRED)

set (LIBRARIES ${PZ_LIBRARIES})

#adds directories
include_directories(
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/Moments"
    "${PROJECT_SOURCE_DIR}/Mass"
    "${PROJECT_SOURCE_DIR}/Stiffness"
    "${PROJECT_SOURCE_DIR}/Elements"
    "${PROJECT_SOURCE_DIR}/Tests"
    "${PROJECT_SOURCE_DIR}/Quadra"
    "${PROJECT_SOURCE_DIR}/Derivatives"
    "${PROJECT_SOURCE_DIR}/Evaluation"
    "${PZ_INCLUDE_DIRS}"
)

#adds dependencies directories
add_subdirectory(Moments)
add_subdirectory(Mass)
add_subdirectory(Stiffness)
add_subdirectory(Elements)
add_subdirectory(Quadra)
add_subdirectory(Derivatives)

option(BUILD_TESTS "Whether or not to build project tests" ON)
if (BUILD_TESTS)
    add_subdirectory(Tests)
endif()

option(OLD_TESTS "Add the old test projects to the build")
if (OLD_TESTS)
    add_subdirectory(Tests/old)
endif()

set(lib_headers ${lib_headers} "BernsteinFEM.h")

add_library(BernsteinFEM STATIC ${lib_sources} ${lib_headers})