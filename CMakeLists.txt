cmake_minimum_required(VERSION 3.12)
project(toyexamples LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo." FORCE)
endif()

add_subdirectory(StructuredLib)
add_subdirectory(Examples)