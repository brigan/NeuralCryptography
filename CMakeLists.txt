cmake_minimum_required(VERSION 3.3)
project(NeuralCryptography)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES helper.h runTPM.cpp TPM.h geneticAttackInfrastructure.h)
add_executable(NeuralCryptography ${SOURCE_FILES})

#include_directories(calendars)

add_subdirectory(nc_tests)