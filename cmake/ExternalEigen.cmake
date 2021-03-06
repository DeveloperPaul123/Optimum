#### Muri Glo Project #####
# Copyright Xstrahl Inc.
# Author: Paul Tsouchlos

# Download and set up Eigen
include(cmake/Externals.cmake)

ExternalProject_Add(eigen
        PREFIX ${EXTERNAL_BUILD_DIR}/eigen
        BINARY_DIR ${EXTERNAL_BUILD_DIR}/eigen
        #GIT_REPOSITORY "https://github.com/RLovelett/eigen.git"
        #GIT_TAG "a46d2e7337c4656f00abe54a8115f6d76153a048"
        URL "https://bitbucket.org/eigen/eigen/get/f3a22f35b044.tar.gz"
        INSTALL_DIR ${EXTERNAL_INSTALL_DIR}
        CMAKE_CACHE_ARGS
            -DCMAKE_BUILD_TYPE:STRING=Release
            -DCMAKE_VERBOSE_MAKEFILE:BOOL=OFF
            -DCMAKE_INSTALL_MESSAGE:BOOL=LAZY
            -DCMAKE_INSTALL_PREFIX:STRING=${EXTERNAL_INSTALL_DIR}
            -DBUILD_TESTING:BOOL=OFF
)

set(EIGEN_DIR ${EXTERNAL_INSTALL_DIR}/include/eigen3)