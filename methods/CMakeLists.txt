cmake_minimum_required(VERSION 3.14)
project(reflex_triangulation C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optimization flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")

# Include directories
include_directories(include)

# Seidel library
set(SEIDEL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../paper/CGATv2/code/baselines/Seidel/Seidel")
add_library(seidel STATIC
    "${SEIDEL_DIR}/construct.c"
    "${SEIDEL_DIR}/misc.c"
    "${SEIDEL_DIR}/monotone.c"
    "${SEIDEL_DIR}/tri.c"
)
target_include_directories(seidel PUBLIC "${SEIDEL_DIR}")
target_compile_definitions(seidel PUBLIC SEGSIZE=100000)

# Reflex triangulation CLI
add_executable(reflex_cli src/reflex_cli.cpp)
target_link_libraries(reflex_cli PRIVATE seidel)

# Bucket triangulation CLI
add_executable(bucket_cli src/bucket_cli.cpp)

# Set output directory
set_target_properties(reflex_cli bucket_cli PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

