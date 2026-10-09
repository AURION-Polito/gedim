# CMake Utilities

This directory contains CMake utility scripts used to support the development, maintenance, and build processes of the **PolyDiM** and **GeDiM** libraries.

> **Note:** These files are primarily intended for PolyDiM and GeDiM developers. Modify them with care, as changes may affect code formatting, static analysis, license management, or cross-platform builds.

## Available Utilities

### Code Formatting - `clang-format.cmake`

Configures the automatic formatting of C++ source (`.cpp`) and header (`.hpp`) files using **ClangFormat**.

This utility helps maintain a consistent coding style across the libraries.

### License Management - `prepend_license.cmake` and `prepend_license.sh`

Provide utilities for automatically prepending the required license header to source files.

These scripts help ensure consistent license information throughout the codebase.

### Static Code Analysis - `cppcheck.cmake`

Configures **Cppcheck** to perform static analysis on the library's C++ source files.

This utility helps developers identify potential bugs, coding issues, and other problems during development.

### Windows Cross-Compilation - `windows_mingw.cmake`

Provides a CMake toolchain configuration for building the project for Windows using the **MinGW** compiler toolchain.

To use this configuration, specify the toolchain file when running CMake:

```bash
cmake -DCMAKE_TOOLCHAIN_FILE=<path_to_windows_mingw.cmake> <other_cmake_options>
```

Replace `<path_to_windows_mingw.cmake>` with the path to the toolchain file and `<other_cmake_options>` with any additional CMake configuration options.

> **Important:** The toolchain file should be specified during the initial CMake configuration, before the compiler and target platform are detected.