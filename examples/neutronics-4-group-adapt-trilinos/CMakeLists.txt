if(NOT H2D_REAL)
	return()
endif(NOT H2D_REAL)

project(neutronics-4-group-adapt-trilinos)
set(HERMES ${HERMES_REAL_BIN})

add_executable(${PROJECT_NAME} main.cpp h1_adapt_norm_maxim.cpp)
include (../CMake.common)

