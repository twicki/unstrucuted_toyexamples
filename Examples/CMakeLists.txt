include_directories(${meshLibrary_SOURCE_DIR})

add_executable(diffusion_explicit diffExplcit.cpp)
target_link_libraries(diffusion_explicit meshLib)

add_executable(diffusion_on_edges diffOnEdges.cpp)
target_link_libraries(diffusion_on_edges meshLib)
