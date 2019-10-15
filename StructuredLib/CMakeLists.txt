project(toyLibrary)
add_library(toylib
            grid.cpp
            grid.h
            )
target_include_directories(toylib PUBLIC .)