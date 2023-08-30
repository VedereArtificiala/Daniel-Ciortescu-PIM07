cmake_minimum_required(VERSION 3.0)
project(imageProject)

find_package(OpenCV REQUIRED)

add_executable(imageProject main.cpp functions.cpp)
    
target_link_libraries(imageProject ${OpenCV_LIBS})
