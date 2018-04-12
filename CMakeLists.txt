cmake_minimum_required(VERSION 3.4.0)

project(BernsteinFEM)

set(CMAKE_BUILD_TYPE Debug)

#adds directories
include_directories(
    "${PROJECT_SOURCE_DIR}/Moments"
    "${PROJECT_SOURCE_DIR}/Mass"
    "${PROJECT_SOURCE_DIR}/Stiffness"
    "${PROJECT_SOURCE_DIR}/Elements"
    "${PROJECT_SOURCE_DIR}/Tests"
    "${PROJECT_SOURCE_DIR}/Quadra"
)

set(SOURCEDIR ${PROJECT_SOURCE_DIR})

#adds dependencies directories
add_subdirectory(Moments)
add_subdirectory(Mass)
add_subdirectory(Stiffness)
add_subdirectory(Elements)
add_subdirectory(Quadra)

option(BUILD_TESTS "Whether or not to build project tests")
if (BUILD_TESTS)
    add_subdirectory(Tests)
endif()

file(GLOB sources *.cpp)

file(GLOB headers *.h)

#find_package(pz)
#target_link_libraries(MAIN Moments_lib Mass_lib Stiff_lib Elem_lib)