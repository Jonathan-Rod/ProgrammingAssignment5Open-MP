cmake_minimum_required(VERSION 3.10)

project(ProgrammingAssignment1Open-MP C)

set(CMAKE_C_STANDARD 17)

set(CMAKE_C_FLAGS_DEBUG "-g -O0 -DDEBUG")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -s")

include_directories(include)

file(GLOB SRC_FILES src/*.c)
add_executable(${PROJECT_NAME} ${SRC_FILES})

find_package(OpenMP REQUIRED)

target_link_libraries(
    ${PROJECT_NAME}
    OpenMP::OpenMP_C
)