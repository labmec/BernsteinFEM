cmake_minimum_required(VERSION 3.4.0)

file(GLOB sources ./*.cpp)

file(GLOB headers *.h)

set(mass ${sources} PARENT_SCOPE)

add_library(Mass_lib STATIC ${sources} ${headers})

add_dependencies(Mass_lib Moments_lib)
target_link_libraries(Mass_lib Moments_lib)