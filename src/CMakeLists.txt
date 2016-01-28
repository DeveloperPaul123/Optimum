### Nedler Mead Project ###
# Copyright Paul Tsouchlos
# Author: Paul T
# Last Update: December 9, 2015

cmake_minimum_required(VERSION 3.3)

set(OPTIMUM_VERSION_MAJOR "0")
set(OPTIMUM_VERSION_MINOR "1")
set(OPTIMUM_VERSION_PATCH "0")

set(OPTIMUM_VERSION "${OPTIMUM_VERSION_MAJOR}.${OPTIMUM_VERSION_MINOR}.${OPTIMUM_VERSION_PATCH}")

project(Optimum)


set(CMAKE_CXX_FLAGS "-Wall")

set(HEADER_FILES neldermead.h geneticalgorithm.h utils.h icp.h)
set(SOURCE_FILES 
	main.cpp)

option(test "Build all test cases." OFF)

set(EIGEN_DIR CACHE PATH "Eigen library directory.")
if(NOT EIGEN_DIR) 
	message(FATAL_ERROR "Define the library directory.")
endif()

include_directories(${EIGEN_DIR})

add_executable(Optimum WIN32 ${HEADER_FILES} ${SOURCE_FILES})

set_source_files_properties(Optimum PROPERTIES LINER_LANGUAGE "CXX")
set_target_properties(Optimum PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")

if(test)
	set( gtest_force_shared_crt ON CACHE BOOL "Always use msvcrt.dll" )
	add_subdirectory(lib/googletest-release-1.7.0)
	enable_testing()
	include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
	
	add_executable(runUnitTests test.cpp)
	
	target_link_libraries(runUnitTests gtest gtest_main)
	
	add_test(runUnitTests runUnitTests)
	
endif(test)

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