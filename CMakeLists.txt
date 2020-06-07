cmake_minimum_required(VERSION 3.15)
project(qm_monte_carlo)

set(CMAKE_CXX_STANDARD 20)

file(GLOB_RECURSE files "src/*.cpp")

add_executable(qm_monte_carlo ${files})
