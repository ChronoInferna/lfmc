# Tests for lfmc

include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.11.0   # Specify the version you want
)
FetchContent_MakeAvailable(Catch2)

# Find Catch2 (if locally installed, otherwise use FetchContent)
# find_package(Catch2 3 REQUIRED)

include(CTest)
include(Catch)

# Register test executables
add_executable(test_timer test_timer.cpp)
target_link_libraries(test_timer PRIVATE lfmc Catch2::Catch2WithMain)

# Register entire file
add_test(NAME test_timer COMMAND test_timer)

# Register individual test cases
# catch_discover_tests(
#     test_timer
# )
