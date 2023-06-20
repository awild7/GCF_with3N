# CMakeLists for GCF_Generator_Suite

# Set up
cmake_minimum_required (VERSION 2.6)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
find_package(ROOT REQUIRED)

project (GCF_Generator_Suite)

if(NOT CMAKE_CXX_FLAGS)
       set(CMAKE_CXX_FLAGS "-std=c++11")
endif()

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include)
set(LIBRARIES ${CMAKE_SOURCE_DIR}/libraries)

include_directories(${INCLUDES} ${LIBRARIES})

add_subdirectory(programs)
