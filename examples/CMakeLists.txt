cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(examples)

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} optimum)
