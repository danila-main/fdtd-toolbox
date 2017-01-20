
/*! \file 
* \brief Memory allocation routines  
*
* Выделяет память для массивов электромагнитных полей, проводимости, токов, источников вторичных электронов, концентрации ионов.
* \author A.Kryukov anton.krv@gmail.com
*/
#include <stdlib.h>

#include "globals.h"
#include "fdtd.h"
#include "data.h"

/*! Выделение памяти для массивов физических величин. */
void allocateCells()
{
    int i, j, k;
    grid _grd = &fdtd_grd;

    PARMS_malloc(fdtd_gi, _grd->nx+2, int**);
    PARMS_malloc(fdtd_Sig, _grd->nx+2, double**);

    for (i = 0; i < _grd->nx+2; ++i) {
        PARMS_malloc(fdtd_gi[i], _grd->ny+2, int*);
        PARMS_malloc(fdtd_Sig[i], _grd->ny+2, double*);
        for (j = 0; j < _grd->ny+2; ++j) {
            PARMS_malloc(fdtd_gi[i][j], _grd->nz+2, int);
            PARMS_malloc(fdtd_Sig[i][j], _grd->nz+2, double);

            memset(fdtd_gi[i][j], 0, (_grd->nz+2)*sizeof(int));
            memset(fdtd_Sig[i][j], 0, (_grd->nz+2)*sizeof(double));
        }
    }

    PARMS_malloc(fdtd_Ex, _grd->nx, double**);
    PARMS_malloc(fdtd_Jx, _grd->nx, double**);
    for (i = 0; i < _grd->nx; ++i) {
        PARMS_malloc(fdtd_Ex[i], _grd->ny+1, double*);
        PARMS_malloc(fdtd_Jx[i], _grd->ny+1, double*);
        for (j = 0; j < _grd->ny+1; ++j) {
            PARMS_malloc(fdtd_Ex[i][j], _grd->nz+1, double);
            PARMS_malloc(fdtd_Jx[i][j], _grd->nz+1, double);

            memset(fdtd_Ex[i][j], 0, (_grd->nz+1)*sizeof(double));
            memset(fdtd_Jx[i][j], 0, (_grd->nz+1)*sizeof(double));
        }
    }

    PARMS_malloc(fdtd_Ey, _grd->ny, double**);
    PARMS_malloc(fdtd_Jy, _grd->ny, double**);
    for (j = 0; j < _grd->ny; ++j) {
        PARMS_malloc(fdtd_Ey[j], _grd->nz+1, double*);
        PARMS_malloc(fdtd_Jy[j], _grd->nz+1, double*);
        for (k = 0; k < _grd->nz+1; ++k) {
            PARMS_malloc(fdtd_Ey[j][k], _grd->nx+1, double);
            PARMS_malloc(fdtd_Jy[j][k], _grd->nx+1, double);

            memset(fdtd_Ey[j][k], 0, (_grd->nx+1)*sizeof(double));
            memset(fdtd_Jy[j][k], 0, (_grd->nx+1)*sizeof(double));
        }
    }

    PARMS_malloc(fdtd_Ez, _grd->nz, double**);
    PARMS_malloc(fdtd_Jz, _grd->nz, double**);
    for (k = 0; k < _grd->nz; ++k) {
        PARMS_malloc(fdtd_Ez[k], _grd->nx+1, double*);
        PARMS_malloc(fdtd_Jz[k], _grd->nx+1, double*);
        for (i = 0; i < _grd->nx+1; ++i) {
            PARMS_malloc(fdtd_Ez[k][i], _grd->ny+1, double);
            PARMS_malloc(fdtd_Jz[k][i], _grd->ny+1, double);

            memset(fdtd_Ez[k][i], 0, (_grd->ny+1)*sizeof(double));
            memset(fdtd_Jz[k][i], 0, (_grd->ny+1)*sizeof(double));
        }
    }

    PARMS_malloc(fdtd_Hx, _grd->nx+1, double**);
    for (i = 0; i < _grd->nx+1; ++i) {
        PARMS_malloc(fdtd_Hx[i], _grd->ny, double*);
        for (j = 0; j < _grd->ny; ++j) {
            PARMS_malloc(fdtd_Hx[i][j], _grd->nz, double);

            memset(fdtd_Hx[i][j], 0, (_grd->nz)*sizeof(double));
        }
    }

    PARMS_malloc(fdtd_Hy, _grd->ny+1, double**);
    for (j = 0; j < _grd->ny+1; ++j) {
        PARMS_malloc(fdtd_Hy[j], _grd->nz, double*);
        for (k = 0; k < _grd->nz; ++k) {
            PARMS_malloc(fdtd_Hy[j][k], _grd->nx, double);

            memset(fdtd_Hy[j][k], 0, (_grd->nx)*sizeof(double));
        }
    }

    PARMS_malloc(fdtd_Hz, _grd->nz+1, double**);
    for (k = 0; k < _grd->nz+1; ++k) {
        PARMS_malloc(fdtd_Hz[k], _grd->nx, double*);
        for (i = 0; i < _grd->nx; ++i) {
            PARMS_malloc(fdtd_Hz[k][i], _grd->ny, double);

            memset(fdtd_Hz[k][i], 0, (_grd->ny)*sizeof(double));
        }
    }

    if (fdtd_comp.l_gran == BND_CPML) {
        allocatePML();
        set_PML_profile();
    }

    if (fdtd_comp.l_tfsf) {
        init_plane_wave();
    }

    if (fdtd_comp.l_nfff) {
        init_nfff_domain();
    }

    if (fdtd_comp.l_disp) {
        init_debye_data();
    }
}

void zeroCurrents()
{
    int i, j, k;
    grid _grd = &fdtd_grd;

    for (i = 0; i < _grd->nx; ++i) {
        for (j = 0; j < _grd->ny+1; ++j) {
            memset(fdtd_Jx[i][j], 0, (_grd->nz+1)*sizeof(double));
        }
    }

    for (j = 0; j < _grd->ny; ++j) {
        for (k = 0; k < _grd->nz+1; ++k) {
            memset(fdtd_Jy[j][k], 0, (_grd->nx+1)*sizeof(double));
        }
    }

    for (k = 0; k < _grd->nz; ++k) {
        for (i = 0; i < _grd->nx+1; ++i) {
            memset(fdtd_Jz[k][i], 0, (_grd->ny+1)*sizeof(double));
        }
    }
}
