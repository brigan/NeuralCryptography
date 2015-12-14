cmake_minimum_required(VERSION 3.3)
project(NeuralCryptography)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES helper.h runTPM.cpp TPM.cpp geneticAttackInfrastructure.h)
add_executable(RunTPM ${SOURCE_FILES})

set(GA_SOURCE_FILES helper.h TPM.cpp geneticAttackInfrastructure.cpp runGeneticAttack.cpp)
add_executable(RunGeneticAttack ${GA_SOURCE_FILES})

#include_directories(calendars)

add_subdirectory(nc_tests)