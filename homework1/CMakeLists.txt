cmake_minimum_required(VERSION 3.10.2)
project(LU_decomposition VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)

find_package(Boost REQUIRED)

include_directories(include)

add_executable(
    ${PROJECT_NAME}
        src/main.cpp
        src/my_matrix.cpp
)

target_link_libraries(
    ${PROJECT_NAME}
        PRIVATE
        ${Boost_INCLUDE_DIRS}
)