### Optimum Project ###
# Copyright Paul Tsouchlos

cmake_minimum_required(VERSION 3.14)

project(optimum)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

option(OPTIMUM_TESTS "Build all test cases." OFF)
option(OPTIMUM_EXAMPLES "Build examples." OFF)

include(Dependencies)

add_subdirectory(optimum)

if(OPTIMUM_EXAMPLES)
  add_subdirectory(examples)
endif(OPTIMUM_EXAMPLES)

# set up documentation
set(DOXYGEN_OUTPUT_DIR "${CMAKE_BINARY_DIR}/docs")
