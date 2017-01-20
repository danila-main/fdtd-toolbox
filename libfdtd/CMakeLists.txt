cmake_minimum_required(VERSION 3.0)

include_directories(./include)
SET (LIBFDTD_SOLVER_cpp solver/debye.c
                        solver/drude.c
                        solver/fdtd.c
                        solver/fdtd_cpml.c
    )

add_library(fdtd STATIC ${LIBFDTD_SOLVER_cpp})
