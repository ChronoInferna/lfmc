# Source files for the lfmc library
set(LFMC_SOURCES
    timing.cpp
)

# Create library
add_library(lfmc ${LFMC_SOURCES})

# Set target properties
target_include_directories(lfmc
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_compile_features(lfmc PUBLIC cxx_std_23)

add_library(project_warnings INTERFACE)

# Comiler-specific warning flags
target_compile_options(project_warnings INTERFACE
  $<$<CXX_COMPILER_ID:GNU,Clang>:
    -Wall
    -Wextra
    -Wpedantic
    -Wconversion
    -Wshadow
  >
  $<$<CXX_COMPILER_ID:MSVC>:
    /W4
    /permissive-
  >
)

target_link_libraries(lfmc PRIVATE project_warnings)

# Config-specific flags
target_compile_options(lfmc
    PRIVATE
        # GCC/Clang
        $<$<AND:$<CONFIG:Debug>,$<CXX_COMPILER_ID:GNU,Clang>>:-O0 -g>
        $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:GNU,Clang>>:-O3 -DNDEBUG>

        # MSVC
        $<$<AND:$<CONFIG:Debug>,$<CXX_COMPILER_ID:MSVC>>:/Od /Zi>
        $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:MSVC>>:/O2 /DNDEBUG>
)


# Installation rules
install(TARGETS lfmc
    EXPORT lfmc-targets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/lfmc
    DESTINATION include
)
