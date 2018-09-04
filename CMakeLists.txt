cmake_minimum_required(VERSION 3.4.0)

project(BernsteinFEM CXX)

find_package(Armadillo REQUIRED)

if(WIN32)
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    set (${ARMADILLO_LIBRARIES} ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    #append blas and lapack libraries into armadillo when on windows
endif()

#adds directories
include_directories(
    "${PROJECT_SOURCE_DIR}/Moments"
    "${PROJECT_SOURCE_DIR}/Mass"
    "${PROJECT_SOURCE_DIR}/Stiffness" # deprecated
    "${PROJECT_SOURCE_DIR}/Elements"
    "${PROJECT_SOURCE_DIR}/Tests"
    "${PROJECT_SOURCE_DIR}/Quadra"
    "${PROJECT_SOURCE_DIR}/Derivatives"
    ${ARMADILLO_INCLUDE_DIRS}
    "${PROJECT_SOURCE_DIR}/Armadillo/include"
)

#adds dependencies directories
add_subdirectory(Moments)
add_subdirectory(Mass)
add_subdirectory(Stiffness) #deprecated
add_subdirectory(Elements)
add_subdirectory(Quadra)
add_subdirectory(Derivatives)

option(BUILD_TESTS "Whether or not to build project tests" ON)
if (BUILD_TESTS)
    add_subdirectory(Tests)
endif()

file(GLOB sources *.cpp)

file(GLOB headers *.h)

#find_package(pz)
#target_link_libraries(MAIN Moments_lib Mass_lib Stiff_lib Elem_lib)