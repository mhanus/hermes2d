project(hermes_common)
cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(PythonLibs REQUIRED)
find_package(NumPy REQUIRED)

include(UseCython)
# if the user specified PYTHON_INCLUDE_PATH, let's use that, otherwise let's
# try to guess using cmake/FindPythonLibs.cmake
if (NOT PYTHON_INCLUDE_PATH)
    find_package(PythonLibs REQUIRED)
endif (NOT PYTHON_INCLUDE_PATH)
find_package(NumPy REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})

enable_testing()

#PYTHONPATH=${hermes_common_SOURCE_DIR}
set(HERMES_COMMON hermes_common)





set(SRC
    _hermes_common_api_new.cpp
    matrix.cpp
    python_solvers.cpp
    python_api.cpp
    )


if(MSVC)
	add_library(${HERMES_COMMON} STATIC ${SRC})
else(MSVC)
	add_library(${HERMES_COMMON} SHARED ${SRC})
endif(MSVC)
target_link_libraries(${HERMES_COMMON} ${PYTHON_LIBRARIES})
set_source_files_properties(matrix.cpp PROPERTIES
    OBJECT_DEPENDS ${hermes_common_SOURCE_DIR}/_hermes_common_api_new.h
    )
add_custom_command(
    OUTPUT _hermes_common_api.h
    COMMAND cython _hermes_common.pyx
    )
add_subdirectory(tests)

cython_add_module(_hermes_common stdcython.cpp)

if(MSVC)
    IF(NOT PYTHON_VERSION)
        execute_process(COMMAND python -c "import sys;print '%d.%d' % sys.version_info[:2]" OUTPUT_VARIABLE PYTHON_VERSION)
        string(STRIP ${PYTHON_VERSION} PYTHON_VERSION)
    ENDIF(NOT PYTHON_VERSION)

    get_filename_component(THIS_DIR . ABSOLUTE) 

    STRING(REPLACE "/" "\\" THIS_DIR ${THIS_DIR})   #function get_filename_component returns the path in Unix format, so we need to adjust it

    add_custom_command(
	TARGET _hermes_common
	POST_BUILD

	#needed so that the function PyImport_Import can find the _hermes_common library
	COMMAND if exist Debug\\_hermes_common.dll copy /Y Debug\\_hermes_common.dll Debug\\_hermes_common.pyd 
	COMMAND if exist Release\\_hermes_common.dll copy /Y Release\\_hermes_common.dll Release\\_hermes_common.pyd

	#this is a very tricky thing on windows, we set the varible after build, but the variable's value is not updated in the current process
	#and for the same reason it is not possible to set it completely using cmake
        COMMAND setx PYTHONPATH \"%PYTHONPATH%\;${THIS_DIR}\\Debug\;${THIS_DIR}\\Release\" /M
    )
    #if the variable exists, set it empty
    execute_process(COMMAND if defined PYTHONPATH set PYTHONPATH=)
    #fill the variable from windows registry
    execute_process(COMMAND setx PYTHONPATH /K HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\${PYTHON_VERSION}\\PythonPath /M)
endif(MSVC)

target_link_libraries(_hermes_common ${HERMES_COMMON})
set_property(TARGET _hermes_common PROPERTY COMPILE_FLAGS "-g")
