project(error-matrix-calculation)
cmake_minimum_required(VERSION 2.8)

# put binary output in bin/ subdirectory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
set(LIBRARY_OUTPUT_PATH ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# select release build type as default
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "'${CMAKE_BUILD_TYPE}' built type was selected.")

# add include directories
include_directories(SYSTEM "${PROJECT_SOURCE_DIR}")

# add compiler flags to define C version
# and to compile position independent code
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c89 -fPIC")

# add compiler flags to enable strict compilation in debug mode
#set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -O0 --coverage")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O2 ")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -Wextra -pedantic")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wmissing-prototypes")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wstrict-prototypes")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wfatal-errors")

# add compiler flags to define C++ version
# and to compile position independent code
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fPIC ")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 --coverage")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread ")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")

# add C++ compiler flag to enable code coverage testing
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} --coverage")

set(HEADERS
    ${PROJECT_SOURCE_DIR}/src/tasks.hh
    ${PROJECT_SOURCE_DIR}/src/include.hh
    ${PROJECT_SOURCE_DIR}/src/util/split.hh
    ${PROJECT_SOURCE_DIR}/src/classes/BinnedTypedMatrix.hh
)

# add source files of this module
set(SOURCES
    ${PROJECT_SOURCE_DIR}/src/main.cpp
    ${PROJECT_SOURCE_DIR}/src/tasks/assemble.cpp
    ${PROJECT_SOURCE_DIR}/src/tasks/rebin.cpp
    ${PROJECT_SOURCE_DIR}/src/tasks/invert.cpp
    ${PROJECT_SOURCE_DIR}/src/tasks/triangulate.cpp
    ${PROJECT_SOURCE_DIR}/src/tasks/print_stats.cpp
    ${PROJECT_SOURCE_DIR}/src/tasks/apply.cpp
    ${PROJECT_SOURCE_DIR}/src/tasks/fit_values.cpp
    ${PROJECT_SOURCE_DIR}/src/tasks/rebin_invert.cpp
    ${PROJECT_SOURCE_DIR}/src/tasks/de_pileup.cpp
    ${PROJECT_SOURCE_DIR}/src/tasks/subtract.cpp
    ${PROJECT_SOURCE_DIR}/src/tasks/multiply_resp_curve.cpp
    ${PROJECT_SOURCE_DIR}/src/classes/BinnedTypedMatrix.cpp
)

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

#find_package(dlib REQUIRED)
#include_directories(${dlib_INCLUDE_DIRS})

# add object files to be compiled
add_executable(error-matrix-calculation ${SOURCES})
target_link_libraries(error-matrix-calculation "stdc++fs")# "c++experimental")
#target_link_libraries(timepix-data-crunching ${Boost_LIBRARIES}) #${Boost_ublas_LIBRARY} )
#target_link_libraries(timepix-data-crunching dlib::dlib)

