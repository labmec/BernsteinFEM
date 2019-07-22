target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
						   $<INSTALL_INTERFACE:.> 
						  )

file(GLOB sources *.cpp)

file(GLOB headers *.h)

set(lib_sources ${lib_sources} ${sources} PARENT_SCOPE)
set(lib_headers ${lib_headers} ${headers} PARENT_SCOPE)