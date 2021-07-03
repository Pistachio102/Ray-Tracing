cmake_minimum_required(VERSION 3.19)
project(Ray_Tracing)

set(CMAKE_CXX_STANDARD 14)

add_executable(Ray_Tracing 1605102_main.cpp 1605102_classes.h)
set(OpenGlLinkers -lglut -lGLU -lGL)
target_link_libraries(Ray_Tracing ${OpenGlLinkers})