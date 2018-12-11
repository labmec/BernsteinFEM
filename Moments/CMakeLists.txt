cmake_minimum_required(VERSION 3.0.0)

file(GLOB sources *.cpp)

file(GLOB headers Moments.h)

set(lib_sources ${lib_sources} ${sources} PARENT_SCOPE)
set(lib_headers ${lib_headers} ${headers} PARENT_SCOPE)