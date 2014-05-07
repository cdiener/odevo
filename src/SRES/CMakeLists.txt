find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

set(CMAKE_BUILD_TYPE Release)

# Set addtitional dependencies and flags
set(LIBS ${LIBS} ${MPI_LIBRARIES})
set(LIBS ${LIBS} -Wl,-O1,--sort-common,--as-needed,-z,relro)
add_definitions(-Wall -Wextra 
				-fstack-protector --param=ssp-buffer-size=4 -D_FORTIFY_SOURCE=2)

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
