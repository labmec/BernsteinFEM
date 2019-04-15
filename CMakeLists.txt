cmake_minimum_required(VERSION 3.0.0)

project(BernsteinFEM CXX)

find_package(Armadillo REQUIRED)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

#append blas and lapack libraries into armadillo
set (${ARMADILLO_LIBRARIES} ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})


#adds directories
include_directories(
    "${PROJECT_SOURCE_DIR}/Moments"
    "${PROJECT_SOURCE_DIR}/Mass"
    "${PROJECT_SOURCE_DIR}/Stiffness"
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
add_subdirectory(Stiffness)
add_subdirectory(Elements)
add_subdirectory(Quadra)
add_subdirectory(Derivatives)

option(BUILD_TESTS "Whether or not to build project tests" ON)
if (BUILD_TESTS)
    add_subdirectory(Tests)
endif()
option(OLD_TESTS "Add the old test projects to the output build")
if (OLD_TESTS)
add_subdirectory(Tests/old)
endif()

add_library(BernsteinFEM STATIC ${lib_sources} ${lib_headers})

#find_package(pz)
#target_link_libraries(MAIN Moments_lib Mass_lib Stiff_lib Elem_lib)
