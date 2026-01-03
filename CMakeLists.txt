cmake_minimum_required(VERSION 3.28)

project(lfmc
    VERSION 0.1.0
    DESCRIPTION "Lock-free Monte Carlo"
    LANGUAGES CXX
)

# Set C++23 standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Export compile commands for clangd and other tools
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Options
option(LFMC_BUILD_TESTS "Build tests" ON)
option(LFMC_BUILD_EXAMPLES "Build examples" ON)
option(LFMC_BUILD_DOCS "Build documentation" OFF)

# Add subdirectories
add_subdirectory(src)

if(LFMC_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(LFMC_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(LFMC_BUILD_DOCS)
    add_subdirectory(docs)
endif()
