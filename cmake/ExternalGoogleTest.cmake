#### Muri Glo Project #####
# Copyright Xstrahl Inc.
# Author: Paul Tsouchlos

# This file will download and configure google testing
# for use within MuriGlo. This is required if you want 
# to do unit testing within MuriGlo. 
cmake_minimum_required(VERSION 2.8.8)

include(cmake/Externals.cmake)

ExternalProject_Add(googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    CMAKE_ARGS -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG:PATH=DebugLibs
               -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE:PATH=ReleaseLibs
               -DCMAKE_CXX_FLAGS=${MSVC_COMPILER_DEFS}
               -Dgtest_force_shared_crt=ON
               -DBUILD_GTEST=ON
    PREFIX "${EXTERNAL_BUILD_DIR}"
	# Disable install step
    INSTALL_COMMAND ""
)

# Specify include dir
ExternalProject_Get_Property(googletest source_dir)
set(GTEST_INCLUDE_DIRS ${source_dir}/googletest/include)

# Specify MainTest's link libraries
ExternalProject_Get_Property(googletest binary_dir)
set(GTEST_LIBS_DIR ${binary_dir}/googlemock/gtest)