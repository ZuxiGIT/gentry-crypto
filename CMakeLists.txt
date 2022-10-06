cmake_minimum_required(VERSION 3.2)

project(gentry-cryptosystem LANGUAGES CXX)

add_compile_options(-Wall -Wextra -Werror -fno-omit-frame-pointer)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-ggdb -fsanitize=address)
    add_link_options(-fsanitize=address)
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-O3)
endif()

add_executable(server
    server/server.cpp)

add_executable(client
    client/client.cpp)
