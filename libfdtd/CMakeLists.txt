cmake_minimum_required(VERSION 3.0)

include_directories(./include)
SET (LIBFDTD_SOLVER_cpp solver/debye.c
                        solver/drude.c
                        solver/fdtd.c
                        solver/fdtd_cpml.c
                        solver/fdtd_mur.c
                        solver/fdtd_nfff.c
                        solver/fdtd_tfsf.c
                        utils/calc_sig_cur.c
                        utils/cells.c
                        utils/globals.c
                        utils/medium.c
                        utils/reports.c
                        utils/time.c
                        utils/load/setpar.c
                        utils/save/auto_save.c
                        utils/save/save_dat.c
                        #utils/save/save_vtr.cpp
    )

add_library(fdtd STATIC ${LIBFDTD_SOLVER_cpp})
