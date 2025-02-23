cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)

# Explicitly list the test source code and headers. The Catch header-only unit
# test framework is stored in with the test source.

set(RT_SIM_VIS_TEST_SRC
    test_main.cpp
)
set(RT_SIM_VIS_TEST_HEADER
    catch.hpp
)

PREPEND(CRT_SIM_VIS_TEST_HEADER)

add_executable(TestRTSimVis ${RT_SIM_VIS_TEST_SRC} ${RT_SIM_VIS_TEST_SRC})

enable_testing()

# Add tests manually here
#add_test(NAME TestRTSimVis  COMMAND TestRTVimSis)

target_link_libraries(TestRTSimVis RTSimVis)