# CMakeLists for GCF_Generator_Suite

# Set up
cmake_minimum_required (VERSION 2.6)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
find_package(ROOT REQUIRED)

project (GCF_Generator_Suite)

set(CMAKE_CXX_FLAGS "-std=c++11")

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include)
include_directories(${INCLUDES})
