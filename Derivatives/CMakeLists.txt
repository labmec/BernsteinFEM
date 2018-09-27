cmake_minimum_required(VERSION 3.4.0)

file(GLOB sources *.cpp)

file(GLOB headers Elem.h)

set(deriv ${sources} PARENT_SCOPE)

add_library(Deriv_lib STATIC ${sources} ${headers})

set(LIBRARIES ${LIBRARIES} Elem_lib)

add_dependencies(Deriv_lib Moments_lib Quadra_lib)
target_link_libraries(Deriv_lib Quadra_lib Moments_lib ${ARMADILLO_LIBRARIES})