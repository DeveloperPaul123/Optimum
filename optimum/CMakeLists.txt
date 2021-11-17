cmake_minimum_required(VERSION 3.21 FATAL_ERROR)

project(optimum
    VERSION 0.10.0
    LANGUAGES CXX
)

CPMAddPackage(
    NAME Eigen
    VERSION 3.4.0
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    # Eigen's CMakelists are not intended for library use
    DOWNLOAD_ONLY YES 
  )
  
if(Eigen_ADDED)
  add_library(Eigen3 INTERFACE IMPORTED)
  target_include_directories(Eigen3 INTERFACE ${Eigen_SOURCE_DIR})
endif()

# PackageProject.cmake will be used to make our target installable
CPMAddPackage("gh:TheLartians/PackageProject.cmake@1.6.0")

set(project_headers
    include/optimum/icp.h
    include/optimum/geneticalgorithm.h
    include/optimum/neldermead.h
    include/optimum/optimizer.h
    include/optimum/utils.h
    
    include/optimum/evolutionary/genetic.h
)

add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME} 
    INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>" 
)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_20)
# being a cross-platform target, we enforce standards conformance on MSVC
target_compile_options(${PROJECT_NAME} INTERFACE "$<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/permissive->")
target_link_libraries(${PROJECT_NAME} INTERFACE Eigen3)

# the location where the project's version header will be placed should match the project's regular
# header paths
string(TOLOWER ${PROJECT_NAME}/version.h VERSION_HEADER_LOCATION)

packageProject(
    NAME ${PROJECT_NAME}
    VERSION ${PROJECT_VERSION}
    NAMESPACE "dp"
    BINARY_DIR ${PROJECT_BINARY_DIR}
    INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include
    INCLUDE_DESTINATION include/${PROJECT_NAME}-${PROJECT_VERSION}
    VERSION_HEADER "${VERSION_HEADER_LOCATION}"
    COMPATIBILITY SameMajorVersion
    DEPENDENCIES "Eigen3"
)

if(OPTIMUM_TESTS)
    CPMAddPackage(
        NAME googletest
        GITHUB_REPOSITORY google/googletest
        GIT_TAG release-1.11.0
        VERSION 1.11.0
        OPTIONS
            "INSTALL_GTEST OFF"
            "gtest_force_shared_crt ON"
    )

    set(PROJECT_TEST_NAME "${PROJECT_NAME}.tests")
    file(GLOB_RECURSE test_sources "${CMAKE_CURRENT_SOURCE_DIR}/test/*.cpp")
    add_executable(${PROJECT_TEST_NAME} ${test_sources})
    target_link_libraries(${PROJECT_TEST_NAME} 
        PUBLIC
            ${PROJECT_NAME}
            gtest
            gtest_main
    )

    add_test(${PROJECT_TEST_NAME} COMMAND ${PROJECT_TEST_NAME})

endif(OPTIMUM_TESTS)
