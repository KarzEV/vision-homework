cmake_minimum_required(VERSION 3.10.2)
project(LU_decomposition VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 11)

include_directories(include)

add_executable(
    ${PROJECT_NAME}
        src/main.cpp
        src/my_matrix.cpp
)