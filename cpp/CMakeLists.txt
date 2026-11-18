# Sapling-CN scoring core (sapling.cn._core).
#
#   cmake -S cpp -B build -DFASTPPM_CN_DIR=<fastppm-cn checkout> -DPython_EXECUTABLE=<python>
#   cmake --build build -j
#
# fastppm-CN (the copy-number branch of fastppm) is built from source as a library; only the
# targets the module needs are built. The module is written to sapling/cn/.
cmake_minimum_required(VERSION 3.16)
project(sapling_cn_core LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# spdlog (bundled with fastppm-cn) builds under C++20 only with std::format
set(SPDLOG_USE_STD_FORMAT ON CACHE BOOL "" FORCE)

set(FASTPPM_CN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../fastppm-cn" CACHE PATH "fastppm-cn source directory")
add_subdirectory(${FASTPPM_CN_DIR} fastppm-cn EXCLUDE_FROM_ALL)  # also provides pybind11

pybind11_add_module(_core MODULE core.cpp)
target_include_directories(_core PRIVATE ${FASTPPM_CN_DIR}/src)
target_link_libraries(_core PRIVATE fastppm-lib)
set_target_properties(_core PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../sapling/cn)
