cmake_minimum_required(VERSION 2.6)

project(optim)

#Some Paths config
set(OPTIM_SRC ${CMAKE_SOURCE_DIR}/src)

include_directories(${OPTIM_SRC})
include_directories(${OPTIM_SRC}/SRES)
add_subdirectory(${OPTIM_SRC}/SRES)
include_directories(${OPTIM_SRC}/odepack)
add_subdirectory(${OPTIM_SRC}/odepack)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

set(CMAKE_BUILD_TYPE Release)

add_definitions(-std=c++11)

# Build the program
add_executable(test src/test.cxx src/lsoda.cxx src/optim.cxx)
target_link_libraries(test SRES odepack gfortran ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(test PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(test PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
