### Nedler Mead Project ###
# Copyright Paul Tsouchlos
# Author: Paul T
# Last Update: December 9, 2015

cmake_minimum_required(VERSION 3.3)

project(OptimumLib)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

option(OPTIMUM_TESTS "Build all test cases." OFF)
option(OPTIMUM_EXAMPLES "Build examples." OFF)

include(Dependencies)

add_subdirectory(optimum)

if(OPTIMUM_EXAMPLES)
	add_executable(OptimumExamples examples/main.cpp)
	target_link_libraries(OptimumExamples Optimum)
	include_directories(${EIGEN_DIR})
endif(OPTIMUM_EXAMPLES)
#Add target to build documents from visual studio.
set(DOXYGEN_INPUT ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile)
set(DOXYGEN_OUTPUT_DIR "${CMAKE_BINARY_DIR}/docs")
find_package(Doxygen)

if(DOXYGEN_FOUND)
	file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR})
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
	${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

ADD_CUSTOM_TARGET(DOCUMENTATION
  COMMAND ${CMAKE_COMMAND} -E echo_append "Building API Documentation..."
  COMMAND ${CMAKE_COMMAND} -E make_directory ${DOXYGEN_OUTPUT_DIR}
  COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND ${CMAKE_COMMAND} -E echo "Done."
  WORKING_DIRECTORY ${DOXYGEN_OUTPUT_DIR})

endif(DOXYGEN_FOUND)
