find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

set(CMAKE_BUILD_TYPE Release)

add_library(SRES ESES.c ESSRSort.c sharefunc.c ESES.h ESSRSort.h sharefunc.h)
target_link_libraries(SRES ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(SRES PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(SRES PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
