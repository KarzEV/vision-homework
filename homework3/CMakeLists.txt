cmake_minimum_required(VERSION 3.10.2)
project(data_processing VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED)

include_directories(include)

add_executable(
        ${PROJECT_NAME}
        src/main.cpp
        src/UniversalSolver.cpp
        src/test.cpp
)

target_link_libraries(
        ${PROJECT_NAME}
        PRIVATE
        Eigen3::Eigen
)