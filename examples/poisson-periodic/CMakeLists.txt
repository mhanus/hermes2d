if(NOT H2D_REAL)
	return()
endif(NOT H2D_REAL)

project(poisson-periodic)

set(HERMES ${HERMES_REAL_BIN})

add_executable(${PROJECT_NAME} main.cpp H1PeriodicSpace.cpp)
include (../CMake.common)
