cmake_minimum_required(VERSION 3.12)
project(toyexamples LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
endif(NOT CMAKE_BUILD_TYPE)

add_subdirectory(GridLib)
add_subdirectory(Examples)