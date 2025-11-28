cmake_minimum_required(VERSION 3.14)
project(reflex_triangulation CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optimization flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")

# Include directories
include_directories(include)

# Reflex triangulation CLI
add_executable(reflex_cli src/reflex_cli.cpp)

# Set output directory
set_target_properties(reflex_cli PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

