find_package(JUDY REQUIRED)
include_directories(${JUDY_INCLUDE_DIR})

# examples
add_subdirectory(smooth-iso)
add_subdirectory(smooth-aniso-x)
add_subdirectory(smooth-aniso-y)
add_subdirectory(lshape)
add_subdirectory(layer)
add_subdirectory(layer-2)
add_subdirectory(line-singularity)
add_subdirectory(kellogg)
add_subdirectory(bessel)
add_subdirectory(screen)
add_subdirectory(neutronics-heat-conduction)
add_subdirectory(neutronics-heat-conduction-adapt)
add_subdirectory(neutronics-2-group-adapt)

