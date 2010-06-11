if(NOT H2D_REAL)
	return()
endif(NOT H2D_REAL)

project(neutronics-2-group-adapt-exact)
set(HERMES ${HERMES_REAL_BIN})

add_executable(${PROJECT_NAME} main.cpp h1_adapt_custom_norm.cpp TikZGraph.cpp)
include (../CMake.common)

