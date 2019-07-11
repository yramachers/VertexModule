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
add_library(Fitter SHARED vertex_module.h 
		   vertex_library.h 
		   vertex_library.cpp 
		   vertex_module.cpp)
target_link_libraries(Vertex PUBLIC Falaise::FalaiseModule ${ROOT_LIBRARIES})

# Build a dynamic library from our sources for algorithm only
# required for testing executables.
include_directories(${ROOT_INCLUDE_DIRS})
add_library(FitAlgorithm SHARED vertex_library.h 
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
#include_directories(. testing)
#add_executable(linetest testing/linetest.cpp)
#target_link_libraries(linetest PUBLIC Catch FitAlgorithm)

# CTest comand
#add_test(NAME linetest
#  COMMAND linetest -s
#)

