cmake_minimum_required(VERSION 3.22)
project(FFT)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) 

add_executable(${PROJECT_NAME} main.cpp classFFT.h)
