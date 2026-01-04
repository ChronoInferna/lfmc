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
    test_timer.cpp
    test_process.cpp
    test_scheme.cpp
)

add_executable(tests ${TEST_SOURCES})
target_link_libraries(tests PRIVATE lfmc Catch2::Catch2WithMain)

catch_discover_tests(
    tests
    ADD_TAGS_AS_LABELS
)
