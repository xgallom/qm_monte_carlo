cmake_minimum_required(VERSION 3.15)
project(qm_monte_carlo)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
-Wall -Wextra -Werror -Wmultiple-inheritance -Wuseless-cast -Wimplicit-fallthrough=5 \
-march=native \
")

set(RELEASE_FLAGS "\
-Ofast \
-flto  -fuse-linker-plugin \
-fno-exceptions -fno-rtti -fno-signed-zeros -fno-trapping-math -frename-registers -ftree-parallelize-loops=4 \
")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${RELEASE_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} ${RELEASE_FLAGS}")

set(CMAKE_AR "gcc-ar")
set(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> qcs <TARGET> <LINK_FLAGS> <OBJECTS>")
set(CMAKE_CXX_ARCHIVE_FINISH true)

file(GLOB_RECURSE files "src/*.cpp")

add_executable(qm_monte_carlo ${files})
