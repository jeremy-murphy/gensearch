cmake_minimum_required(VERSION 2.6)
project(gensearch)
include_directories(counting_components)
add_executable(test_search test_search.cpp)
# add_executable(count_long_search count_long_search.cpp)
add_executable(experimental_search experimental_search.cpp)
add_executable(test_long_search test_long_search.cpp)
add_executable(test_word_search test_word_search.cpp)
add_executable(time_long_search time_long_search.cpp)
add_executable(time_word_search time_word_search.cpp)

set(Boost_NO_SYSTEM_PATHS ON)
find_package(Boost REQUIRED filesystem system test_exec_monitor)
set(BOOST_TEST_LIBS ${Boost_LIBRARIES})
include_directories(${Boost_INCLUDE_DIRS})

add_executable(search_test2 "$ENV{BOOST_ROOT}/libs/algorithm/test/search_test2.cpp")
target_link_libraries(search_test2 ${BOOST_TEST_LIBS})
add_executable(search_test3 "$ENV{BOOST_ROOT}/libs/algorithm/test/search_test3.cpp")
target_link_libraries(search_test3 ${BOOST_TEST_LIBS})
add_executable(search_test5 "$ENV{BOOST_ROOT}/libs/algorithm/test/search_test5.cpp")
target_link_libraries(search_test5 ${BOOST_TEST_LIBS})
