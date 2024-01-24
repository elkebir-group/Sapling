cmake_minimum_required(VERSION 3.4)
project(BBT)




set(CMAKE_CXX_STANDARD 20)
set(CXX 1)

set (CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/")
find_package(GUROBI REQUIRED)
MESSAGE(STATUS "Include dir of Gurobi: ${GUROBI_INCLUDE_DIRS}")
MESSAGE(STATUS "C library for Gurobi: ${GUROBI_LIBRARY}")
MESSAGE(STATUS "CXX library for Gurobi: ${GUROBI_CXX_DEBUG_LIBRARY}")

find_package(Boost)

include_directories(${BOOST_INCLUDE_DIRS} ${GUROBI_INCLUDE_DIRS})

set( CommonLibs m ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY} )


link_libraries(${CommonLibs} Boost::boost pybind11::pybind11)

add_executable(BBT src/Input.h src/Input.cpp src/AncestryGraph.cpp src/AncestryGraph.h src/BBT.cpp src/BBT.h src/main.cpp src/LLH.cpp src/LLH.h)

add_subdirectory(pybind11)
pybind11_add_module(pyBBT src/Input.h src/Input.cpp src/AncestryGraph.cpp src/AncestryGraph.h src/BBT.cpp src/BBT.h  src/pymain.cpp)

## EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
## define (VERSION_INFO) here.
#target_compile_definitions(pyBBT
#        PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})