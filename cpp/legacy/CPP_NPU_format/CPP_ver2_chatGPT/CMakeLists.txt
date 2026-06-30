cmake_minimum_required(VERSION 3.14)
project(npu_writer LANGUAGES CXX)

# Find pybind11 (either pip install pybind11 and use the installed CMake files,
# or vendor pybind11 in a subdirectory).
find_package(pybind11 CONFIG REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

pybind11_add_module(npu_writer npu_writer.cpp)
