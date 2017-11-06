
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

set_source_files_properties(Optimum PROPERTIES LINER_LANGUAGE "CXX")
target_compile_features(Optimum INTERFACE cxx_auto_type)
add_custom_target(optimumIde SOURCES ${project_headers})
include_directories(${EIGEN_DIR})
if(OPTIMUM_TESTS)
	set(PROJECT_TEST_NAME "OptimumTests")
	enable_testing()	
	add_executable(${PROJECT_TEST_NAME} test/test.cpp)
	target_include_directories(${PROJECT_TEST_NAME} PRIVATE ${GTEST_INCLUDE_DIRS})
	if(NOT WIN32)
	  target_link_libraries(${PROJECT_TEST_NAME}
	    ${GTEST_LIBS_DIR}/libgtest.a
	    ${GTEST_LIBS_DIR}/libgtest_main.a
	    )
	else()
	  target_link_libraries(${PROJECT_TEST_NAME}
	    debug ${GTEST_LIBS_DIR}/DebugLibs/${CMAKE_FIND_LIBRARY_PREFIXES}gtestd${CMAKE_FIND_LIBRARY_SUFFIXES}
	    optimized ${GTEST_LIBS_DIR}/ReleaseLibs/${CMAKE_FIND_LIBRARY_PREFIXES}gtest${CMAKE_FIND_LIBRARY_SUFFIXES}
	    )
	  target_link_libraries(${PROJECT_TEST_NAME}
	    debug ${GTEST_LIBS_DIR}/DebugLibs/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_maind${CMAKE_FIND_LIBRARY_SUFFIXES}
	    optimized ${GTEST_LIBS_DIR}/ReleaseLibs/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main${CMAKE_FIND_LIBRARY_SUFFIXES}
	    )
	endif()
	target_link_libraries(${PROJECT_TEST_NAME} Optimum)
	add_test(${PROJECT_TEST_NAME} ${PROJECT_TEST_NAME})
	add_dependencies(${PROJECT_TEST_NAME} googletest)
endif(OPTIMUM_TESTS)