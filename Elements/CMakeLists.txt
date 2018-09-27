cmake_minimum_required(VERSION 3.4.0)

file(GLOB sources *.cpp)

file(GLOB headers Elem.h)

set(elem ${sources} PARENT_SCOPE)

add_library(Elem_lib STATIC ${sources} ${headers})

set(LIBRARIES ${LIBRARIES} Elem_lib)

add_dependencies(Elem_lib Mass_lib Stiff_lib)
target_link_libraries(Elem_lib Mass_lib Stiff_lib Moments_lib ${ARMADILLO_LIBRARIES})