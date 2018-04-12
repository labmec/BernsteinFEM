cmake_minimum_required(VERSION 3.4.0)

file(GLOB sources ./*.cpp)

file(GLOB headers *.h)

set(stiff ${sources} PARENT_SCOPE)

add_library(Stiff_lib STATIC ${sources} ${headers})

add_dependencies(Stiff_lib Moments_lib)
target_link_libraries(Stiff_lib Moments_lib)