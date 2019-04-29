cmake_minimum_required(VERSION 3.13)
project(qm_monte_carlo)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
-fno-exceptions \
-Wall \
-Wextra \
-Wshadow \
-Wnon-virtual-dtor \
-Wold-style-cast \
-Wcast-align \
-Wunused \
-Woverloaded-virtual \
-Wpedantic \
-Wsign-conversion \
-Wmisleading-indentation \
-Wduplicated-cond \
-Wlogical-op \
-Wnull-dereference \
-Wuseless-cast \
-Wdouble-promotion \
-Wformat=2 \
-Weffc++ \
-pedantic \
")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O4 -mavx -ftree-vectorize -g")

add_subdirectory(src)
