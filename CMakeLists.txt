cmake_minimum_required(VERSION 3.0)

project(SOMA_CLP_CEC21) 

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
find_package( Threads )

add_executable(SOMA_CLP_CEC21 SOMA_CLP_CEC21.cpp kmeans.cpp cec21_test_func.cpp)
target_link_libraries( SOMA_CLP_CEC21 ${CMAKE_THREAD_LIBS_INIT} )