cmake_minimum_required(VERSION 3.4.0)

file(GLOB sources ./*.cpp)

file(GLOB headers *.h)

set(quadra ${sources} PARENT_SCOPE)

add_library(Quadra_lib STATIC ${sources} ${headers})

set(LIBRARIES ${LIBRARIES} Quadra_lib)