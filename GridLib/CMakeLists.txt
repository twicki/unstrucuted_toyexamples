project(meshLibrary)
add_library(meshLib STATIC
                vertex.h
                cell.h
                cell.cpp
                edge.h
                edge.cpp
                grid.h
                grid.cpp
)
target_include_directories(meshLib PUBLIC .)