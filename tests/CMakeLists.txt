cmake_minimum_required(VERSION 3.16)
project(pairhmm-test)

file(GLOB_RECURSE SOURCE_FILES ${CMAKE_SOURCE_DIR}/tests/*.cpp ${CMAKE_SOURCE_DIR}/src/pairhmm/*.cpp ${CMAKE_SOURCE_DIR}/src/table/*.cpp ${CMAKE_SOURCE_DIR}/src/utils/*.cpp)

add_executable(pairhmm-test ${SOURCE_FILES})
target_link_libraries(pairhmm-test ${Boost_LIBRARIES} )
target_link_libraries(pairhmm-test biovoltron)
target_include_directories(pairhmm-test PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_include_directories(pairhmm-test PRIVATE ${CMAKE_SOURCE_DIR}/submodules/biovoltron/include)

add_test(NAME pairhmm COMMAND pairhmm-test)
