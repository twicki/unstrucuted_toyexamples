#examples based on atlas
find_package(eckit)
find_package(atlas)

add_executable(autogenMesh atlasMesh.cpp)
target_link_libraries(autogenMesh atlas)

add_executable(printData atlasPrint.cpp)
target_link_libraries(printData atlas)

add_executable(neighborAccess atlasGauss.cpp)
target_link_libraries(neighborAccess atlas)

add_executable(diffAtlas atlasDiff.cpp)
target_link_libraries(diffAtlas atlas)
