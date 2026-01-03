# lfmc

## Description

A modern C++23 library for lock-free Monte Carlo simulations.

## Project Structure

```
lfmc/
├── include/lfmc/    # Public header files
├── src/             # Source files
├── tests/           # Unit tests
├── examples/        # Example applications
├── docs/            # Documentation
├── cmake/           # CMake modules
└── CMakeLists.txt   # Main CMake configuration
```

## Requirements

- C++23 compatible compiler (GCC 12+, Clang 16+, or MSVC 2022+)
- CMake 3.25 or higher

## Building

There are currently 6 presets: base, debug, release, clang-debug, gcc-debug, and msvc-debug.

```bash
cmake --preset <preset-name>
cmake --build --preset <preset-name>
```

### Build Options

- `LFMC_BUILD_TESTS`: Build tests (default: ON)
- `LFMC_BUILD_EXAMPLES`: Build examples (default: ON)
- `LFMC_BUILD_DOCS`: Build documentation (default: OFF)

### Testing

```bash
ctest --preset <preset-name>
```

## Usage

Include the library in your CMake project:

```cmake
find_package(lfmc REQUIRED)
target_link_libraries(your_target PRIVATE lfmc)
```

## Tech Stack

| Technology     | Purpose                                     |
| -------------- | ------------------------------------------- |
| C++23          | C++23 standard for modern tools             |
| Ninja          | Build system generator for efficient builds |
| CMake          | Build configuration and management          |
| Catch2         | Unit testing framework                      |
| Doxygen        | Documentation generation                    |
| GitHub Actions | CI/CD                                       |

## License

See LICENSE file for details.
