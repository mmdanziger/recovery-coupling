cmake_minimum_required(VERSION 2.6)
project(adjlistgen)

add_executable(adjlistgen main.cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Ofast -Wall -Wno-unused-local-typedefs")

install(TARGETS adjlistgen RUNTIME DESTINATION bin)
