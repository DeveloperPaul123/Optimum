### Nedler Mead Project ###
# Copyright Paul Tsouchlos
# Author: Paul T
# Last Update: December 9, 2015

cmake_minimum_required(VERSION 3.3)

set(OPTIMUM_VERSION_MAJOR "0")
set(OPTIMUM_VERSION_MINOR "1")
set(OPTIMUM_VERSION_PATCH "0")

project(Optimum)

set(HEADER_FILES optimizer.h neldermead.h geneticalgorithm.h)
set(SOURCE_FILES main.cpp)

set(EIGEN_DIR CACHE PATH "Eigen library directory.")
if(NOT EIGEN_DIR) 
	message(FATAL_ERROR "Define the library directory.")
endif()

include_directories(${EIGEN_DIR})

add_executable(Optimum WIN32 ${HEADER_FILES} ${SOURCE_FILES})

set_source_files_properties(Optimum PROPERTIES LINER_LANGUAGE "CXX")
set_target_properties(Optimum PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")

#Add target to build documents from visual studio.
set(DOXYGEN_INPUT ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile)
set(DOXYGEN_OUTPUT_DIR CACHE PATH "Path for document generation output")
find_package(Doxygen)

if(DOXYGEN_FOUND)
	file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR})

ADD_CUSTOM_TARGET(DOCUMENTATION
  COMMAND ${CMAKE_COMMAND} -E echo_append "Building API Documentation..."
  COMMAND ${CMAKE_COMMAND} -E make_directory ${DOXYGEN_OUTPUT_DIR}
  COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT}
  COMMAND ${CMAKE_COMMAND} -E echo "Done."
  WORKING_DIRECTORY ${DOXYGEN_OUTPUT_DIR})

endif(DOXYGEN_FOUND)