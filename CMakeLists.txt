cmake_minimum_required(VERSION 3.22)
project(YadroTest)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) 

#find_package(Qt5 REQUIRED COMPONENTS Charts Widgets)

add_executable(${PROJECT_NAME} main.cpp classFFT.h)

#target_link_libraries(${PROJECT_NAME} Qt5::Charts Qt5::Widgets)