set(BUILD_OUTPUT_DIR ${CMAKE_BINARY_DIR}/bin)
file(MAKE_DIRECTORY ${BUILD_OUTPUT_DIR})

add_library(baseline_common STATIC
    src/polygon_io.cpp
)
target_include_directories(baseline_common
    PUBLIC
      ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Earcut CLI
add_executable(earcut_cli src/earcut_cli.cpp)
target_link_libraries(earcut_cli PRIVATE baseline_common)
target_include_directories(earcut_cli PRIVATE
    ${CMAKE_SOURCE_DIR}/external/earcut.hpp/include
)

# PolyTri CLI (Seidel)
add_executable(polytri_cli src/polytri_cli.cpp)
target_link_libraries(polytri_cli PRIVATE baseline_common)
target_include_directories(polytri_cli PRIVATE
    ${CMAKE_SOURCE_DIR}/external/polytri/include
)
target_compile_definitions(polytri_cli PRIVATE POLYTRI_IMPLEMENTATION)

# CGAL CLI (Constrained Delaunay example)
find_package(CGAL REQUIRED COMPONENTS Core)
add_executable(cgal_polygon_cli src/cgal_polygon_cli.cpp)
target_link_libraries(cgal_polygon_cli PRIVATE baseline_common CGAL::CGAL CGAL::CGAL_Core)

set_target_properties(earcut_cli polytri_cli cgal_polygon_cli
    PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY ${BUILD_OUTPUT_DIR}
)
