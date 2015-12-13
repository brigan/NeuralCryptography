include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(runDevLangTests basic_check.cpp)

target_link_libraries(runDevLangTests gtest gtest_main)