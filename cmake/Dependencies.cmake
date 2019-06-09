find_package(Eigen3 CONFIG REQUIRED)


if(OPTIMUM_TESTS)
	find_package(GTest MODULE)
	if(NOT GTest_FOUND)
		include(Externals)	
		include(ExternalGoogleTest)
	endif(NOT GTest_FOUND)
endif(OPTIMUM_TESTS)

find_package(Doxygen)
