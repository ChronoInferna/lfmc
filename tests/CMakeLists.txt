# Tests for lfmc

include(FetchContent)

# TODO CI has to install the entirety of Catch2 which takes a lot of time - is there an optimization?
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.11.0   # Specify the version you want
)
FetchContent_MakeAvailable(Catch2)

# Find Catch2 (if locally installed, otherwise use FetchContent above)
# find_package(Catch2 3 REQUIRED)

include(CTest)
include(Catch)

set(TEST_SOURCES
    test_timing.cpp
    test_process.cpp
    test_scheme.cpp
    test_pipeline.cpp
    test_adaptive.cpp
    test_strategies.cpp
    test_engine.cpp
    # test_convergency.cpp
)

add_executable(tests ${TEST_SOURCES})
target_link_libraries(tests PRIVATE lfmc Catch2::Catch2WithMain)

# Standalone benchmark executable (not a Catch2 test — has its own main)
add_executable(bench bench_asvr_vs_fixed.cpp)
target_link_libraries(bench PRIVATE lfmc)

# Quick comparison: Engine vs fixed strategies (same budget)
add_executable(quick_compare quick_compare.cpp)
target_link_libraries(quick_compare PRIVATE lfmc)
add_test(NAME quick_compare COMMAND quick_compare)

# Copy UCRT64 C++ runtime DLLs so the executables find the right runtime on Windows.
# Without this, Git's bundled MinGW64 DLLs (on PATH earlier) take precedence and
# cause exit code 0xc0000139 (STATUS_ENTRYPOINT_NOT_FOUND).
if(MSYS OR MINGW)
    set(UCRT64_BIN_DIR "C:/msys64/ucrt64/bin")
    foreach(target tests bench)
        foreach(dll libgcc_s_seh-1.dll libstdc++-6.dll libwinpthread-1.dll)
            if(EXISTS "${UCRT64_BIN_DIR}/${dll}")
                add_custom_command(TARGET ${target} POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E copy
                        "${UCRT64_BIN_DIR}/${dll}"
                        "$<TARGET_FILE_DIR:${target}>/${dll}"
                    VERBATIM
                )
            endif()
        endforeach()
    endforeach()
endif()

catch_discover_tests(
    tests
    ADD_TAGS_AS_LABELS
)
