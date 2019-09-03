# - Basic CMake setup
# Check version meets ou requirements
# Declare project, which will configure compiler for us
cmake_minimum_required(VERSION 3.3)
project(VertexExtrapolator)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
find_package(Falaise REQUIRED)

# Need ROOT library separately
find_package(ROOT REQUIRED)

# find the headers
include_directories(.)

# Build a dynamic library from our sources
#add_library(Vertex SHARED vertex_module.h 
#		   vertex_library.h 
#		   vertex_library.cpp 
#		   vertex_module.cpp)
#target_link_libraries(Vertex PUBLIC Falaise::FalaiseModule ${ROOT_LIBRARIES})

# Build a dynamic library from our sources for algorithm only
# required for testing executables.
include_directories(${ROOT_INCLUDE_DIRS})
add_library(VertexAlgorithm SHARED vertex_library.h 
			  vertex_library.cpp)
target_link_libraries(VertexAlgorithm PUBLIC ${ROOT_LIBRARIES})


# Configure example pipeline script for use from the build dir
configure_file("ve.conf.in" "ve.conf" @ONLY)

# Build the testing code, tell CTest about it
enable_testing()
set(CMAKE_CXX_STANDARD 11)

# main catch library compile only once
add_library(Catch testing/catch-main.cpp)
target_include_directories(Catch INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/testing)

# add test executable
include_directories(. testing)
add_executable(vfroml testing/vfroml.cpp)
target_link_libraries(vfroml PUBLIC Catch VertexAlgorithm)

# add test executable
include_directories(. testing)
add_executable(vfromh testing/vfromh.cpp)
target_link_libraries(vfromh PUBLIC Catch VertexAlgorithm)

# add test executable
include_directories(. testing)
add_executable(wirev testing/wirev.cpp)
target_link_libraries(wirev PUBLIC Catch VertexAlgorithm)


# CTest comand
add_test(NAME vfroml
  COMMAND vfroml -s
)

# CTest comand
add_test(NAME vfromh
  COMMAND vfromh -s
)

# CTest comand
add_test(NAME wirev
  COMMAND wirev -s
)

