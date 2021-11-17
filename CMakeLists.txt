### Optimum Project ###
# Copyright Paul Tsouchlos

cmake_minimum_required(VERSION 3.21)

project(optimum)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

option(OPTIMUM_TESTS "Build all test cases." OFF)
option(OPTIMUM_EXAMPLES "Build examples." OFF)

include(cmake/CPM.cmake)

if(OPTIMUM_TESTS)
  enable_testing()
endif()

add_subdirectory(optimum)

if(OPTIMUM_EXAMPLES)
  add_subdirectory(examples)
endif(OPTIMUM_EXAMPLES)
