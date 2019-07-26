include_directories(${meshLibrary_SOURCE_DIR})

add_executable(diffusion_explicit diffExplicit.cpp)
target_link_libraries(diffusion_explicit meshLib)

add_executable(diffusion_on_edges diffOnEdges.cpp)
target_link_libraries(diffusion_on_edges meshLib)

add_executable(laplace_explicit lapExplcit.cpp)
target_link_libraries(laplace_explicit meshLib)

#SET_PROPERTY(DIRECTORY PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${CMAKE_BINARY_DIR}/Examples/*.vtk")
