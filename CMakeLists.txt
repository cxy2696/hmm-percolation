cmake_minimum_required(VERSION 3.16)
project(hp VERSION 0.1.2 LANGUAGES C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB HP_C_SOURCES
    src/core/*.c src/hmm/*.c src/perc/*.c src/market/*.c src/cross/*.c)
file(GLOB HP_APP_SOURCES src/app/*.cpp)

add_library(hp STATIC ${HP_C_SOURCES})
target_include_directories(hp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(hp PUBLIC m)

add_library(hp_shared SHARED ${HP_C_SOURCES})
set_target_properties(hp_shared PROPERTIES OUTPUT_NAME hp)
target_include_directories(hp_shared PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(hp_shared PUBLIC m)

add_library(hpapp STATIC ${HP_APP_SOURCES})
target_include_directories(hpapp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(hpapp PUBLIC hp)

foreach(app hp_synth hp_regime hp_connectedness hp_perc_lab)
    add_executable(${app} apps/${app}.cpp)
    target_link_libraries(${app} PRIVATE hpapp)
endforeach()

foreach(ex discrete_hmm gaussian_hmm percolation)
    add_executable(ex_${ex} examples/${ex}.c)
    target_link_libraries(ex_${ex} PRIVATE hp)
endforeach()

include(CTest)
if(BUILD_TESTING)
    file(GLOB HP_TEST_SOURCES tests/test_*.c)
    add_executable(hp_tests ${HP_TEST_SOURCES})
    target_include_directories(hp_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    target_link_libraries(hp_tests PRIVATE hp)
    add_test(NAME hp_tests COMMAND hp_tests)
    set_tests_properties(hp_tests PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_app.cpp)
        add_executable(hp_tests_app tests/test_app.cpp)
        target_include_directories(hp_tests_app PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
        target_link_libraries(hp_tests_app PRIVATE hpapp)
        add_test(NAME hp_tests_app COMMAND hp_tests_app)
        set_tests_properties(hp_tests_app PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    endif()
endif()
