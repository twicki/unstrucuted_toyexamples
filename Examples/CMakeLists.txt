#examples based on our mesh-library
add_executable(test simple.cpp)
target_link_libraries(test toylib)

add_executable(edgetest edge_handeling.cpp)
target_link_libraries(edgetest toylib)