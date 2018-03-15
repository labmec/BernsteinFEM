
project(BernsteinFEM)

file(GLOB sources *.cpp)

file(GLOB headers *.h)

add_executable(main_exe ${sources} ${headers})

#find_package(pz)
#target_link_libraries(main_exe pz)