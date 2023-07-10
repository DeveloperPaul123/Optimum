
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(Optimum)

set(project_headers
    include/optimum/icp.h
    include/optimum/geneticalgorithm.h
    include/optimum/neldermead.h
    include/optimum/optimizer.h
    include/optimum/utils.h)

add_library(Optimum INTERFACE)
target_include_directories(Optimum 
    INTERFACE 
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>" 
)
target_link_libraries(Optimum INTERFACE Eigen3::Eigen)

set_source_files_properties(Optimum PROPERTIES LINER_LANGUAGE "CXX")
target_compile_features(Optimum INTERFACE cxx_auto_type)
add_custom_target(optimumIde SOURCES ${project_headers})
include_directories(${EIGEN_DIR})
if(OPTIMUM_TESTS)
    set(PROJECT_TEST_NAME "OptimumTests")
    enable_testing()	
    add_executable(${PROJECT_TEST_NAME} test/test.cpp)
    target_link_libraries(${PROJECT_TEST_NAME} Optimum GTest::gtest GTest::gtest_main)
    add_test(${PROJECT_TEST_NAME} ${PROJECT_TEST_NAME})
endif(OPTIMUM_TESTS)
