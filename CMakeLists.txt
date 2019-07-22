# @file neopz/BernsteinFEM/CMakeLists.txt  -- CMake file for the BernsteinFEM module


target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
						   $<INSTALL_INTERFACE:.> 
						  )

#adds dependencies directories
add_subdirectory(Moments)
add_subdirectory(Mass)
add_subdirectory(Stiffness)
add_subdirectory(Elements)
add_subdirectory(Quadra)
add_subdirectory(Derivatives)

file(GLOB headers *.h)

set(lib_headers ${lib_headers} ${headers})

install(FILES ${lib_headers} DESTINATION ${NEOPZ_INCLUDE}/BernsteinFEM)

set(pzbernstein_STAT_HDRS ${lib_headers} PARENT_SCOPE)
set(pzbernstein_STAT_SRCS ${lib_sources} PARENT_SCOPE)