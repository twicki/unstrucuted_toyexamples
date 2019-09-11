#examples based on atlas
find_package(eckit)
find_package(atlas)

add_executable(meshGen atlasMeshGen.cpp)
target_link_libraries(meshGen atlas)

add_executable(printer atlasPrint.cpp)
target_link_libraries(printer atlas)

add_executable(gaussian atlasGauss.cpp)
target_link_libraries(gaussian atlas)

add_executable(delauny atlasDelauny.cpp)
target_link_libraries(delauny atlas)

add_executable(diffusion atlasDiff.cpp)
target_link_libraries(diffusion atlas)
