cmake_minimum_required(VERSION 3.10)
project(tfdnoise VERSION 1.0 LANGUAGES CXX)

# Compiler flags [[3]]
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -ffast-math -std=c++17")

# Include directories [[2]]
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Find required libraries [[1]]
find_library(FFTW3F_LIBRARY fftw3f REQUIRED)
find_library(SEGYIO_LIBRARY segyio REQUIRED)

# Executable target [[1]]
add_executable(tfdnoise 
    tfdnoise.cpp
    include/tfd.cpp
)

# Link libraries [[1]]
target_link_libraries(tfdnoise
    ${FFTW3F_LIBRARY}
    ${SEGYIO_LIBRARY}
)

# Installation rules (optional)
install(TARGETS tfdnoise DESTINATION bin)
