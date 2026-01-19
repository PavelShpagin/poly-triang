set(BUILD_OUTPUT_DIR ${CMAKE_BINARY_DIR}/bin)
file(MAKE_DIRECTORY ${BUILD_OUTPUT_DIR})

set(PAPER_BASELINES_DIR ${CMAKE_SOURCE_DIR}/paper/CGAT/code/baselines)

add_library(baseline_common STATIC
    src/polygon_io.cpp
)
target_include_directories(baseline_common
    PUBLIC
      ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# PolyTri CLI (Seidel)
add_executable(polytri_cli src/polytri_cli.cpp)
target_link_libraries(polytri_cli PRIVATE baseline_common)
target_include_directories(polytri_cli PRIVATE
    ${PAPER_BASELINES_DIR}/polytri/include
)

# PolyPartition MONO triangulation baseline CLI (Garey-style monotone partition)
add_executable(polypartition_mono_cli src/polypartition_mono_cli.cpp
    ${PAPER_BASELINES_DIR}/polypartition/src/polypartition.cpp
)
target_include_directories(polypartition_mono_cli PRIVATE
    ${PAPER_BASELINES_DIR}/polypartition/src
)

# PolyPartition Hertel–Mehlhorn convex partition baseline CLI
add_executable(polypartition_hm_cli src/polypartition_hm_cli.cpp
    ${PAPER_BASELINES_DIR}/polypartition/src/polypartition.cpp
)
target_include_directories(polypartition_hm_cli PRIVATE
    ${PAPER_BASELINES_DIR}/polypartition/src
)

# Seidel baseline CLI (palmerc/Seidel C implementation)
add_executable(seidel_cli
    src/seidel_cli.c
    ${PAPER_BASELINES_DIR}/Seidel/Seidel/construct.c
    ${PAPER_BASELINES_DIR}/Seidel/Seidel/misc.c
    ${PAPER_BASELINES_DIR}/Seidel/Seidel/monotone.c
    ${PAPER_BASELINES_DIR}/Seidel/Seidel/tri.c
)
target_include_directories(seidel_cli PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${PAPER_BASELINES_DIR}/Seidel/Seidel
)
target_link_libraries(seidel_cli PRIVATE m)
target_compile_options(seidel_cli PRIVATE
    -include ${CMAKE_CURRENT_SOURCE_DIR}/include/seidel.h
)

# CGAL CLI (Constrained Delaunay example) - optional (not required for paper benchmarks)
find_package(CGAL QUIET COMPONENTS Core)
if (CGAL_FOUND)
  add_executable(cgal_polygon_cli src/cgal_polygon_cli.cpp)
  target_link_libraries(cgal_polygon_cli PRIVATE baseline_common CGAL::CGAL CGAL::CGAL_Core)
  add_executable(cgal_gen_cli src/cgal_gen_cli.cpp)
  target_link_libraries(cgal_gen_cli PRIVATE baseline_common CGAL::CGAL CGAL::CGAL_Core)
  # CGAL Hertel–Mehlhorn approximate convex partition baseline CLI (paper/CGAT)
  add_executable(cgal_hm_cli ${CMAKE_SOURCE_DIR}/paper/CGAT/code/cli/cgal_hm_cli.cpp)
  target_link_libraries(cgal_hm_cli PRIVATE baseline_common CGAL::CGAL CGAL::CGAL_Core)
  set_target_properties(polytri_cli cgal_polygon_cli
      PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${BUILD_OUTPUT_DIR}
  )
  set_target_properties(cgal_gen_cli cgal_hm_cli PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BUILD_OUTPUT_DIR})
else()
  message(STATUS "CGAL not found; skipping cgal_polygon_cli")
  set_target_properties(polytri_cli seidel_cli
      PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${BUILD_OUTPUT_DIR}
  )
endif()

# Ensure seidel_cli output directory in CGAL_FOUND case too
if (CGAL_FOUND)
  set_target_properties(seidel_cli PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BUILD_OUTPUT_DIR})
endif()

# Ensure polypartition_mono_cli output directory always
set_target_properties(polypartition_mono_cli PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BUILD_OUTPUT_DIR})
set_target_properties(polypartition_hm_cli PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BUILD_OUTPUT_DIR})
