
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(optimum)

set(project_headers
	include/optimum/icp.h
	include/optimum/geneticalgorithm.h
	include/optimum/neldermead.h
	include/optimum/optimizer.h
	include/optimum/utils.h)

add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME} 
	INTERFACE 
	"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
	"$<INSTALL_INTERFACE:include>" 
)

set_source_files_properties(${PROJECT_NAME} PROPERTIES LINER_LANGUAGE "CXX")
target_compile_features(${PROJECT_NAME} INTERFACE cxx_auto_type)

add_custom_target(${PROJECT_NAME}.headers SOURCES ${project_headers})

target_link_libraries(${PROJECT_NAME} INTERFACE Eigen3::Eigen)

if(OPTIMUM_TESTS)
	set(PROJECT_TEST_NAME "${PROJECT_NAME}.tests")
	enable_testing()	
	add_executable(${PROJECT_TEST_NAME} test/test.cpp)
	target_link_libraries(${PROJECT_TEST_NAME} 
		PUBLIC ${PROJECT_NAME}
		PRIVATE  GTest::GTest GTest::Main)
	add_test(${PROJECT_TEST_NAME} ${PROJECT_TEST_NAME})

	if(NOT GTest_FOUND)
		add_dependencies(${PROJECT_TEST_NAME} GoogleTestExternal)
	endif(NOT GTest_FOUND)
endif(OPTIMUM_TESTS)

if(DOXYGEN_FOUND)
	doxygen_add_docs(documentation
		${project_headers} ALL
		WORKING_DIRECTORY docs
		COMMENT "Generating documentation"
	)
endif(DOXYGEN_FOUND)
