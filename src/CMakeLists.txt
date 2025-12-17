# Source files for the lfmc library
set(LFMC_SOURCES
    lfmc.cpp
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
