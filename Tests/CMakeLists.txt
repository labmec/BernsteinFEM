cmake_minimum_required(VERSION 3.0.0)

#some tests

add_executable(MASS_TEST mas.cpp)
add_dependencies(MASS_TEST BernsteinFEM)
target_link_libraries(MASS_TEST BernsteinFEM ${ARMADILLO_LIBRARIES})
target_include_directories(MASS_TEST INTERFACE  ${ARMADILLO_INCLUDE_DIRS}  "../Armadillo/include")
set(allTests ${allTests} MASS_TEST)

add_executable(MOMENTS_TEST mom.cpp)
add_dependencies(MOMENTS_TEST BernsteinFEM)
target_link_libraries(MOMENTS_TEST BernsteinFEM ${ARMADILLO_LIBRARIES})
target_include_directories(MOMENTS_TEST INTERFACE  ${ARMADILLO_INCLUDE_DIRS}  "../Armadillo/include")
set(allTests ${allTests} MOMENTS_TEST)

add_executable(STIFF_TEST stif.cpp)
add_dependencies(STIFF_TEST BernsteinFEM)
target_link_libraries(STIFF_TEST BernsteinFEM ${ARMADILLO_LIBRARIES})
target_include_directories(STIFF_TEST INTERFACE  ${ARMADILLO_INCLUDE_DIRS}  "../Armadillo/include")
set(allTests ${allTests} STIFF_TEST)

add_executable(MOMENTS_QUAD_TEST mom-quad.cpp)
add_dependencies(MOMENTS_QUAD_TEST BernsteinFEM)
target_link_libraries(MOMENTS_QUAD_TEST BernsteinFEM ${ARMADILLO_LIBRARIES})
target_include_directories(MOMENTS_QUAD_TEST INTERFACE  ${ARMADILLO_INCLUDE_DIRS}  "../Armadillo/include")
set(allTests ${allTests} MOMENTS_QUAD_TEST)

add_executable(MOMENTS_TRI_TEST mom-tri.cpp)
add_dependencies(MOMENTS_TRI_TEST BernsteinFEM)
target_link_libraries(MOMENTS_TRI_TEST BernsteinFEM ${ARMADILLO_LIBRARIES})
target_include_directories(MOMENTS_TRI_TEST INTERFACE  ${ARMADILLO_INCLUDE_DIRS}  "../Armadillo/include")
set(allTests ${allTests} MOMENTS_TRI_TEST)

add_executable(QUAD_STIF_CPU_TEST cpu_test.cpp)
add_dependencies(QUAD_STIF_CPU_TEST BernsteinFEM)
target_link_libraries(QUAD_STIF_CPU_TEST BernsteinFEM ${ARMADILLO_LIBRARIES})
target_include_directories(QUAD_STIF_CPU_TEST INTERFACE ${ARMADILLO_INCLUDE_DIRS} "../Armadillo/include")

add_executable(QUAD_CONVERGENCE_TEST convergence_test.cpp)
add_dependencies(QUAD_CONVERGENCE_TEST BernsteinFEM)
target_link_libraries(QUAD_CONVERGENCE_TEST BernsteinFEM ${ARMADILLO_LIBRARIES})
target_include_directories(QUAD_CONVERGENCE_TEST INTERFACE ${ARMADILLO_INCLUDE_DIRS} "../Armadillo/include")

add_executable(TRI_STIF_CPU_TEST cpu_test_tri.cpp)
add_dependencies(TRI_STIF_CPU_TEST BernsteinFEM)
target_link_libraries(TRI_STIF_CPU_TEST BernsteinFEM ${ARMADILLO_LIBRARIES})
target_include_directories(TRI_STIF_CPU_TEST INTERFACE ${ARMADILLO_INCLUDE_DIRS} "../Armadillo/include")
