cmake_minimum_required(VERSION 3.12)
project(toyexamples LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)

add_library(meshLib STATIC unstructured_toylib.cpp)

add_executable(laplace lap.cpp)
target_link_libraries(laplace meshLib)