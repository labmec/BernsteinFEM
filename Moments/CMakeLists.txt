cmake_minimum_required(VERSION 3.4.0)

#project()

file(GLOB sources *.cpp)

file(GLOB headers Moments.h)

add_library(Moments_lib STATIC ${sources} ${headers})

add_dependencies(Moments_lib Quadra_lib)
target_link_libraries(Moments_lib Quadra_lib)