cmake_minimum_required (VERSION 3.16)

# Define project name
project(pairhmm)

# Define cpp standard
set(CMAKE_CXX_STANDARD 20)

# Include boost
find_package( Boost 1.71.0 COMPONENTS program_options REQUIRED )
if(NOT Boost_FOUND)
  message(FATAL_ERROR "Boost program_options library not found")
endif()
include_directories(${Boost_INCLUDE_DIRS} )

# Include biovoltron
add_subdirectory(${PROJECT_SOURCE_DIR}/submodules/biovoltron)

# Link your source code to your project
file(GLOB_RECURSE SOURCE_FILES src/*.cpp)
add_executable(pairhmm ${SOURCE_FILES})
target_compile_options(pairhmm PRIVATE -g)
target_include_directories(pairhmm PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(pairhmm PRIVATE ${PROJECT_SOURCE_DIR}/submodules/biovoltron/include)
target_link_libraries(pairhmm PRIVATE ${Boost_LIBRARIES} )
target_link_libraries(pairhmm PRIVATE biovoltron)

#include(CTest)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/tests)