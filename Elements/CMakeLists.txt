cmake_minimum_required(VERSION 3.4.0)

file(GLOB sources *.cpp)

file(GLOB headers Elements.h)

set(deriv ${sources} PARENT_SCOPE)

add_library(Elements_lib STATIC ${sources} ${headers})