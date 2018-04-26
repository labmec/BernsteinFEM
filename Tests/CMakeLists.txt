cmake_minimum_required(VERSION 3.4.0)

file(GLOB sources ./*.cpp)

file(GLOB headers *.h)

set(test ${sources} PARENT_SCOPE)

#some tests
add_executable(ELEM_TEST el.cpp)
add_dependencies(ELEM_TEST Elem_lib)
target_link_libraries(ELEM_TEST Elem_lib)
target_include_directories(ELEM_TEST INTERFACE  ../../externallibs/include)
set(allTests ${allTests} ELEM_TEST)

add_executable(MASS_TEST mas.cpp)
add_dependencies(MASS_TEST Mass_lib)
target_link_libraries(MASS_TEST Mass_lib)
target_include_directories(MASS_TEST INTERFACE  ../../externallibs/include)
set(allTests ${allTests} MASS_TEST)

add_executable(MOMENTS_TEST mom.cpp)
add_dependencies(MOMENTS_TEST Moments_lib)
target_link_libraries(MOMENTS_TEST Moments_lib)
target_include_directories(MOMENTS_TEST INTERFACE  ../../externallibs/include)
set(allTests ${allTests} MOMENTS_TEST)

add_executable(STIFF_TEST stif.cpp)
add_dependencies(STIFF_TEST Stiff_lib)
target_link_libraries(STIFF_TEST Stiff_lib)
target_include_directories(STIFF_TEST INTERFACE  ../../externallibs/include)
set(allTests ${allTests} STIFF_TEST)

add_executable(TEST_1D test1d.cpp)
add_dependencies(TEST_1D Elem_lib)
target_link_libraries(TEST_1D Elem_lib Mass_lib Stiff_lib Moments_lib)
target_include_directories(TEST_1D INTERFACE  ../../externallibs/include)
set(allTests ${allTests} TEST_1D)