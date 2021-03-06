cmake_minimum_required(VERSION 3.0.0)

# this file specifies projects to test and run the routines defined by this library

# Sample that shows how to compute an element's Load Vector
add_executable(LoadVecSample LoadVecSample.cpp)
add_dependencies(LoadVecSample BernsteinFEM)
target_link_libraries(LoadVecSample BernsteinFEM ${ARMADILLO_LIBRARIES})
target_include_directories(LoadVecSample INTERFACE ${ARMADILLO_INCLUDE_DIRS})

# Sample that shows how to compute an element's Mass Matrix
add_executable(MassMatrixSample MassMatrixSample.cpp)
add_dependencies(MassMatrixSample BernsteinFEM)
target_link_libraries(MassMatrixSample BernsteinFEM ${ARMADILLO_LIBRARIES})
target_include_directories(MassMatrixSample INTERFACE  ${ARMADILLO_INCLUDE_DIRS})

# Sample that shows how to compute an element's Stiffness Matrix
add_executable(StiffnessMatrixSample StiffnessMatrixSample.cpp)
add_dependencies(StiffnessMatrixSample BernsteinFEM)
target_link_libraries(StiffnessMatrixSample BernsteinFEM ${ARMADILLO_LIBRARIES})
target_include_directories(StiffnessMatrixSample INTERFACE  ${ARMADILLO_INCLUDE_DIRS})

