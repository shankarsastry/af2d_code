cmake_minimum_required(VERSION 3.24)
project(af2d LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)

# Triangle (libigl mirror)
FetchContent_Declare(
    triangle
    GIT_REPOSITORY https://github.com/libigl/triangle.git
    GIT_TAG        master
)
FetchContent_MakeAvailable(triangle)

# Expose Triangle's headers (its CMakeLists.txt doesn't set include dirs)
target_include_directories(triangle PUBLIC ${triangle_SOURCE_DIR})

# GoogleTest
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.15.2
)
FetchContent_MakeAvailable(googletest)

find_package(OpenMP)

enable_testing()

add_subdirectory(src)
add_subdirectory(tests)

add_executable(af2d_main main.cpp)
target_link_libraries(af2d_main PRIVATE af2d_lib)
