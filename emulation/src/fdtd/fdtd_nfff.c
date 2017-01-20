
/*! \file 
* \brief Функции диаграммы направленности.
* \author Крюков А.А. anton.krv@gmail.com
*/

#include <math.h>
#include <hdf5.h>
#include "../../include/globals.h"
#include "../../include/constants.h"
#include "../../include/functions.h"
#include "../../include/fdtd.h"

extern void append_nfff_data(const char *file_name, const double time);

#define NFFF_GRID_FILE "nfff_grid.h5"

static double *nfff_ex_y_min;/*!< массив z- компоненты поля Ex. */
static double *nfff_ex_y_max;/*!< массив z- компоненты поля Ex. */
static double *nfff_ex_z_min;/*!< массив y- компоненты поля Ex. */
static double *nfff_ex_z_max;/*!< массив y- компоненты поля Ex. */

static double *nfff_ey_x_min;/*!< массив z- компоненты поля Ey. */
static double *nfff_ey_x_max;/*!< массив z- компоненты поля Ey. */
static double *nfff_ey_z_min;/*!< массив x- компоненты поля Ey. */
static double *nfff_ey_z_max;/*!< массив x- компоненты поля Ey. */

static double *nfff_ez_x_min;/*!< массив y- компоненты поля Ez. */
static double *nfff_ez_x_max;/*!< массив y- компоненты поля Ez. */
static double *nfff_ez_y_min;/*!< массив x- компоненты поля Ez. */
static double *nfff_ez_y_max;/*!< массив x- компоненты поля Ez. */

static double *nfff_hx_y_min;/*!< массив z- компоненты поля Hx. */
static double *nfff_hx_y_max;/*!< массив z- компоненты поля Hx. */
static double *nfff_hx_z_min;/*!< массив y- компоненты поля Hx. */
static double *nfff_hx_z_max;/*!< массив y- компоненты поля Hx. */

static double *nfff_hy_x_min;/*!< массив z- компоненты поля Hy. */
static double *nfff_hy_x_max;/*!< массив z- компоненты поля Hy. */
static double *nfff_hy_z_min;/*!< массив x- компоненты поля Hy. */
static double *nfff_hy_z_max;/*!< массив x- компоненты поля Hy. */

static double *nfff_hz_x_min;/*!< массив y- компоненты поля Hz. */
static double *nfff_hz_x_max;/*!< массив y- компоненты поля Hz. */
static double *nfff_hz_y_min;/*!< массив x- компоненты поля Hz. */
static double *nfff_hz_y_max;/*!< массив x- компоненты поля Hz. */

static int nfff_width_xmin = 12;
static int nfff_width_xmax = 12;

static int nfff_width_ymin = 12;
static int nfff_width_ymax = 12;

static int nfff_width_zmin = 5;
static int nfff_width_zmax = 15;

static int nfff_nx;
static int nfff_ny;
static int nfff_nz;

static void _zero_nfff_domen()
{
    //Xmin
    memset(nfff_ey_x_min, 0, nfff_ny*nfff_nz*sizeof(double));
    memset(nfff_ez_x_min, 0, nfff_ny*nfff_nz*sizeof(double));
    memset(nfff_hy_x_min, 0, nfff_ny*nfff_nz*sizeof(double));
    memset(nfff_hz_x_min, 0, nfff_ny*nfff_nz*sizeof(double));

    //Xmax
    memset(nfff_ey_x_max, 0, nfff_ny*nfff_nz*sizeof(double));
    memset(nfff_ez_x_max, 0, nfff_ny*nfff_nz*sizeof(double));
    memset(nfff_hy_x_max, 0, nfff_ny*nfff_nz*sizeof(double));
    memset(nfff_hz_x_max, 0, nfff_ny*nfff_nz*sizeof(double));

    //Ymin
    memset(nfff_ex_y_min, 0, nfff_nx*nfff_nz*sizeof(double));
    memset(nfff_ez_y_min, 0, nfff_nx*nfff_nz*sizeof(double));
    memset(nfff_hx_y_min, 0, nfff_nx*nfff_nz*sizeof(double));
    memset(nfff_hz_y_min, 0, nfff_nx*nfff_nz*sizeof(double));

    //Ymax
    memset(nfff_ex_y_max, 0, nfff_nx*nfff_nz*sizeof(double));
    memset(nfff_ez_y_max, 0, nfff_nx*nfff_nz*sizeof(double));
    memset(nfff_hx_y_max, 0, nfff_nx*nfff_nz*sizeof(double));
    memset(nfff_hz_y_max, 0, nfff_nx*nfff_nz*sizeof(double));

    //Zmin
    memset(nfff_ex_z_min, 0, nfff_nx*nfff_ny*sizeof(double));
    memset(nfff_ey_z_min, 0, nfff_nx*nfff_ny*sizeof(double));
    memset(nfff_hx_z_min, 0, nfff_nx*nfff_ny*sizeof(double));
    memset(nfff_hy_z_min, 0, nfff_nx*nfff_ny*sizeof(double));

    //Zmax
    memset(nfff_ex_z_max, 0, nfff_nx*nfff_ny*sizeof(double));
    memset(nfff_ey_z_max, 0, nfff_nx*nfff_ny*sizeof(double));
    memset(nfff_hx_z_max, 0, nfff_nx*nfff_ny*sizeof(double));
    memset(nfff_hy_z_max, 0, nfff_nx*nfff_ny*sizeof(double));
}

void init_nfff_domain()
{
    grid _grd = &fdtd_grd;

    nfff_nx = _grd->nx - (nfff_width_xmax + nfff_width_xmin) - 1;
    nfff_ny = _grd->ny - (nfff_width_ymax + nfff_width_ymin) - 1;
    nfff_nz = _grd->nz - (nfff_width_zmax + nfff_width_zmin) - 1;

    //Xmin
    PARMS_malloc(nfff_ey_x_min, nfff_ny*nfff_nz, double);
    PARMS_malloc(nfff_ez_x_min, nfff_ny*nfff_nz, double);
    PARMS_malloc(nfff_hy_x_min, nfff_ny*nfff_nz, double);
    PARMS_malloc(nfff_hz_x_min, nfff_ny*nfff_nz, double);

    //Xmax
    PARMS_malloc(nfff_ey_x_max, nfff_ny*nfff_nz, double);
    PARMS_malloc(nfff_ez_x_max, nfff_ny*nfff_nz, double);
    PARMS_malloc(nfff_hy_x_max, nfff_ny*nfff_nz, double);
    PARMS_malloc(nfff_hz_x_max, nfff_ny*nfff_nz, double);

    //Ymin
    PARMS_malloc(nfff_ex_y_min, nfff_nx*nfff_nz, double);
    PARMS_malloc(nfff_ez_y_min, nfff_nx*nfff_nz, double);
    PARMS_malloc(nfff_hx_y_min, nfff_nx*nfff_nz, double);
    PARMS_malloc(nfff_hz_y_min, nfff_nx*nfff_nz, double);

    //Ymax
    PARMS_malloc(nfff_ex_y_max, nfff_nx*nfff_nz, double);
    PARMS_malloc(nfff_ez_y_max, nfff_nx*nfff_nz, double);
    PARMS_malloc(nfff_hx_y_max, nfff_nx*nfff_nz, double);
    PARMS_malloc(nfff_hz_y_max, nfff_nx*nfff_nz, double);

    //Zmin
    PARMS_malloc(nfff_ex_z_min, nfff_nx*nfff_ny, double);
    PARMS_malloc(nfff_ey_z_min, nfff_nx*nfff_ny, double);
    PARMS_malloc(nfff_hx_z_min, nfff_nx*nfff_ny, double);
    PARMS_malloc(nfff_hy_z_min, nfff_nx*nfff_ny, double);

    //Zmax
    PARMS_malloc(nfff_ex_z_max, nfff_nx*nfff_ny, double);
    PARMS_malloc(nfff_ey_z_max, nfff_nx*nfff_ny, double);
    PARMS_malloc(nfff_hx_z_max, nfff_nx*nfff_ny, double);
    PARMS_malloc(nfff_hy_z_max, nfff_nx*nfff_ny, double);

    _zero_nfff_domen();
}

static void _preload_h_field()
{
    int i, j, k;
    grid _grd = &fdtd_grd;
    double val;

    for (j = 0; j < nfff_ny; ++j) {
        for (k = 0; k < nfff_nz; ++k) {
            //Xmin
            val = 0.0;
            val += fdtd_Hy[j + 1 + nfff_width_ymin][k + 1 + nfff_width_zmin][nfff_width_xmin];
            val += fdtd_Hy[j + 1 + nfff_width_ymin][k + nfff_width_zmin][nfff_width_xmin];
            val *= 0.25;
            nfff_hy_x_min[j*nfff_nz + k] += val;

            val = 0.0;
            val += fdtd_Hz[k + 1 + nfff_width_zmin][nfff_width_xmin][j + 1 + nfff_width_ymin];
            val += fdtd_Hz[k + 1 + nfff_width_zmin][nfff_width_xmin][j + nfff_width_ymin];
            val *= 0.25;
            nfff_hz_x_min[j*nfff_nz + k] += val;


            //Xmax
            val = 0.0;
            val += fdtd_Hy[j + 1 + nfff_width_ymin][k + 1 + nfff_width_zmin][_grd->nx - nfff_width_xmax - 1];
            val += fdtd_Hy[j + 1 + nfff_width_ymin][k + nfff_width_zmin][_grd->nx - nfff_width_xmax - 1];
            val *= 0.25;
            nfff_hy_x_max[j*nfff_nz + k] += val;

            val = 0.0;
            val += fdtd_Hz[k + 1 + nfff_width_zmin][_grd->nx - nfff_width_xmax - 1][j + 1 + nfff_width_ymin];
            val += fdtd_Hz[k + 1 + nfff_width_zmin][_grd->nx - nfff_width_xmax - 1][j + nfff_width_ymin];
            val *= 0.25;
            nfff_hz_x_max[j*nfff_nz + k] += val;
        }
    }

    for (i = 0; i < nfff_nx; ++i) {
        for (k = 0; k < nfff_nz; ++k) {
            //Ymin
            val = 0.0;
            val += fdtd_Hx[i + 1 + nfff_width_xmin][nfff_width_ymin][k + 1 + nfff_width_zmin];
            val += fdtd_Hx[i + 1 + nfff_width_xmin][nfff_width_ymin][k + nfff_width_zmin];
            val *= 0.25;
            nfff_hx_y_min[i*nfff_nz + k] += val;

            val = 0.0;
            val += fdtd_Hz[k + 1 + nfff_width_zmin][i + 1 + nfff_width_xmin][nfff_width_ymin];
            val += fdtd_Hz[k + 1 + nfff_width_zmin][i + nfff_width_xmin][nfff_width_ymin];
            val *= 0.25;
            nfff_hz_y_min[i*nfff_nz + k] += val;


            //Ymax
            val = 0.0;
            val += fdtd_Hx[i + 1 + nfff_width_xmin][_grd->ny - nfff_width_ymax - 1][k + 1 + nfff_width_zmin];
            val += fdtd_Hx[i + 1 + nfff_width_xmin][_grd->ny - nfff_width_ymax - 1][k + nfff_width_zmin];
            val *= 0.25;
            nfff_hx_y_max[i*nfff_nz + k] += val;

            val = 0.0;
            val += fdtd_Hz[k + 1 + nfff_width_zmin][i + 1 + nfff_width_xmin][_grd->ny - nfff_width_ymax - 1];
            val += fdtd_Hz[k + 1 + nfff_width_zmin][i + nfff_width_xmin][_grd->ny - nfff_width_ymax - 1];
            val *= 0.25;
            nfff_hz_y_max[i*nfff_nz + k] += val;
        }
    }

    for (i = 0; i < nfff_nx; ++i) {
        for (j = 0; j < nfff_ny; ++j) {
            //Zmin
            val = 0.0;
            val += fdtd_Hx[i + 1 + nfff_width_xmin][j + 1 + nfff_width_ymin][nfff_width_zmin];
            val += fdtd_Hx[i + 1 + nfff_width_xmin][j + nfff_width_ymin][nfff_width_zmin];
            val *= 0.25;
            nfff_hx_z_min[i*nfff_ny + j] += val;

            val = 0.0;
            val += fdtd_Hy[j + 1 + nfff_width_ymin][nfff_width_zmin][i + 1 + nfff_width_xmin];
            val += fdtd_Hy[j + 1 + nfff_width_ymin][nfff_width_zmin][i + nfff_width_xmin];
            val *= 0.25;
            nfff_hy_z_min[i*nfff_ny + j] += val;


            //Zmax
            val = 0.0;
            val += fdtd_Hx[i + 1 + nfff_width_xmin][j + 1 + nfff_width_ymin][_grd->nz - nfff_width_zmax - 1];
            val += fdtd_Hx[i + 1 + nfff_width_xmin][j + nfff_width_ymin][_grd->nz - nfff_width_zmax - 1];
            val *= 0.25;
            nfff_hx_z_max[i*nfff_ny + j] += val;

            val = 0.0;
            val += fdtd_Hy[j + 1 + nfff_width_ymin][_grd->nz - nfff_width_zmax - 1][i + 1 + nfff_width_xmin];
            val += fdtd_Hy[j + 1 + nfff_width_ymin][_grd->nz - nfff_width_zmax - 1][i + nfff_width_xmin];
            val *= 0.25;
            nfff_hy_z_max[i*nfff_ny + j] += val;
        }
    }
}

static void _preload_e_field()
{
    int i, j, k;
    grid _grd = &fdtd_grd;
    double val;

    for (j = 0; j < nfff_ny; ++j) {
        for (k = 0; k < nfff_nz; ++k) {
            //Xmin
            val = 0.0;
            val += fdtd_Ey[j + nfff_width_ymin][k + 1 + nfff_width_zmin][nfff_width_xmin];
            val += fdtd_Ey[j + nfff_width_ymin][k + 1 + nfff_width_zmin][nfff_width_xmin + 1];
            val += fdtd_Ey[j + 1 + nfff_width_ymin][k + 1 + nfff_width_zmin][nfff_width_xmin];
            val += fdtd_Ey[j + 1 + nfff_width_ymin][k + 1 + nfff_width_zmin][nfff_width_xmin + 1];
            val *= 0.25;
            nfff_ey_x_min[j*nfff_nz + k] += val;

            val = 0.0;
            val += fdtd_Ez[k + nfff_width_zmin][nfff_width_xmin][j + 1 + nfff_width_ymin];
            val += fdtd_Ez[k + nfff_width_zmin][nfff_width_xmin + 1][j + 1 + nfff_width_ymin];
            val += fdtd_Ez[k + 1 + nfff_width_zmin][nfff_width_xmin][j + 1 + nfff_width_ymin];
            val += fdtd_Ez[k + 1 + nfff_width_zmin][nfff_width_xmin + 1][j + 1 + nfff_width_ymin];
            val *= 0.25;
            nfff_ez_x_min[j*nfff_nz + k] += val;


            //Xmax
            val = 0.0;
            val += fdtd_Ey[j + nfff_width_ymin][k + 1 + nfff_width_zmin][_grd->nx - nfff_width_xmax - 1];
            val += fdtd_Ey[j + nfff_width_ymin][k + 1 + nfff_width_zmin][_grd->nx - nfff_width_xmax];
            val += fdtd_Ey[j + 1 + nfff_width_ymin][k + 1 + nfff_width_zmin][_grd->nx - nfff_width_xmax - 1];
            val += fdtd_Ey[j + 1 + nfff_width_ymin][k + 1 + nfff_width_zmin][_grd->nx - nfff_width_xmax];
            val *= 0.25;
            nfff_ey_x_max[j*nfff_nz + k] += val;

            val = 0.0;
            val += fdtd_Ez[k + nfff_width_zmin][_grd->nx - nfff_width_xmax - 1][j + 1 + nfff_width_ymin];
            val += fdtd_Ez[k + nfff_width_zmin][_grd->nx - nfff_width_xmax][j + 1 + nfff_width_ymin];
            val += fdtd_Ez[k + 1 + nfff_width_zmin][_grd->nx - nfff_width_xmax - 1][j + 1 + nfff_width_ymin];
            val += fdtd_Ez[k + 1 + nfff_width_zmin][_grd->nx - nfff_width_xmax][j + 1 + nfff_width_ymin];
            val *= 0.25;
            nfff_ez_x_max[j*nfff_nz + k] += val;
        }
    }

    for (i = 0; i < nfff_nx; ++i) {
        for (k = 0; k < nfff_nz; ++k) {
            //Ymin
            val = 0.0;
            val += fdtd_Ex[i + nfff_width_xmin][nfff_width_ymin][k + 1 + nfff_width_zmin];
            val += fdtd_Ex[i + nfff_width_xmin][nfff_width_ymin + 1][k + 1 + nfff_width_zmin];
            val += fdtd_Ex[i + 1 + nfff_width_xmin][nfff_width_ymin][k + 1 + nfff_width_zmin];
            val += fdtd_Ex[i + 1 + nfff_width_xmin][nfff_width_ymin + 1][k + 1 + nfff_width_zmin];
            val *= 0.25;
            nfff_ex_y_min[i*nfff_nz + k] += val;

            val = 0.0;
            val += fdtd_Ez[k + nfff_width_zmin][i + 1 + nfff_width_xmin][nfff_width_ymin];
            val += fdtd_Ez[k + nfff_width_zmin][i + 1 + nfff_width_xmin][nfff_width_ymin + 1];
            val += fdtd_Ez[k + 1 + nfff_width_zmin][i + 1 + nfff_width_xmin][nfff_width_ymin];
            val += fdtd_Ez[k + 1 + nfff_width_zmin][i + 1 + nfff_width_xmin][nfff_width_ymin + 1];
            val *= 0.25;
            nfff_ez_y_min[i*nfff_nz + k] += val;


            //Ymax
            val = 0.0;
            val += fdtd_Ex[i + nfff_width_xmin][_grd->ny - nfff_width_ymax - 1][k + 1 + nfff_width_zmin];
            val += fdtd_Ex[i + nfff_width_xmin][_grd->ny - nfff_width_ymax][k + 1 + nfff_width_zmin];
            val += fdtd_Ex[i + 1 + nfff_width_xmin][_grd->ny - nfff_width_ymax - 1][k + 1 + nfff_width_zmin];
            val += fdtd_Ex[i + 1 + nfff_width_xmin][_grd->ny - nfff_width_ymax][k + 1 + nfff_width_zmin];
            val *= 0.25;
            nfff_ex_y_max[i*nfff_nz + k] += val;

            val = 0.0;
            val += fdtd_Ez[k + nfff_width_zmin][i + 1 + nfff_width_xmin][_grd->ny - nfff_width_ymax - 1];
            val += fdtd_Ez[k + nfff_width_zmin][i + 1 + nfff_width_xmin][_grd->ny - nfff_width_ymax];
            val += fdtd_Ez[k + 1 + nfff_width_zmin][i + 1 + nfff_width_xmin][_grd->ny - nfff_width_ymax - 1];
            val += fdtd_Ez[k + 1 + nfff_width_zmin][i + 1 + nfff_width_xmin][_grd->ny - nfff_width_ymax];
            val *= 0.25;
            nfff_ez_y_max[i*nfff_nz + k] += val;
        }
    }

    for (i = 0; i < nfff_nx; ++i) {
        for (j = 0; j < nfff_ny; ++j) {
            //Zmin
            val = 0.0;
            val += fdtd_Ex[i + nfff_width_xmin][j + 1 + nfff_width_ymin][nfff_width_zmin];
            val += fdtd_Ex[i + nfff_width_xmin][j + 1 + nfff_width_ymin][nfff_width_zmin + 1];
            val += fdtd_Ex[i + 1 + nfff_width_xmin][j + 1 + nfff_width_ymin][nfff_width_zmin];
            val += fdtd_Ex[i + 1 + nfff_width_xmin][j + 1 + nfff_width_ymin][nfff_width_zmin + 1];
            val *= 0.25;
            nfff_ex_z_min[i*nfff_ny + j] += val;

            val = 0.0;
            val += fdtd_Ey[j + nfff_width_ymin][nfff_width_zmin][i + 1 + nfff_width_xmin];
            val += fdtd_Ey[j + nfff_width_ymin][nfff_width_zmin + 1][i + 1 + nfff_width_xmin];
            val += fdtd_Ey[j + 1 + nfff_width_ymin][nfff_width_zmin][i + 1 + nfff_width_xmin];
            val += fdtd_Ey[j + 1 + nfff_width_ymin][nfff_width_zmin + 1][i + 1 + nfff_width_xmin];
            val *= 0.25;
            nfff_ey_z_min[i*nfff_ny + j] += val;


            //Zmax
            val = 0.0;
            val += fdtd_Ex[i + nfff_width_xmin][j + 1 + nfff_width_ymin][_grd->nz - nfff_width_zmax - 1];
            val += fdtd_Ex[i + nfff_width_xmin][j + 1 + nfff_width_ymin][_grd->nz - nfff_width_zmax];
            val += fdtd_Ex[i + 1 + nfff_width_xmin][j + 1 + nfff_width_ymin][_grd->nz - nfff_width_zmax - 1];
            val += fdtd_Ex[i + 1 + nfff_width_xmin][j + 1 + nfff_width_ymin][_grd->nz - nfff_width_zmax];
            val *= 0.25;
            nfff_ex_z_max[i*nfff_ny + j] += val;

            val = 0.0;
            val += fdtd_Ey[j + nfff_width_ymin][_grd->nz - nfff_width_zmax - 1][i + 1 + nfff_width_xmin];
            val += fdtd_Ey[j + nfff_width_ymin][_grd->nz - nfff_width_zmax][i + 1 + nfff_width_xmin];
            val += fdtd_Ey[j + 1 + nfff_width_ymin][_grd->nz - nfff_width_zmax - 1][i + 1 + nfff_width_xmin];
            val += fdtd_Ey[j + 1 + nfff_width_ymin][_grd->nz - nfff_width_zmax][i + 1 + nfff_width_xmin];
            val *= 0.25;
            nfff_ey_z_max[i*nfff_ny + j] += val;
        }
    }
}

void load_nfff_domain()
{
    _zero_nfff_domen();

    _preload_e_field();

    _preload_h_field();
}

static void _write_axis(double *x, double *dx, int n, char *id, hid_t file_id)
{
    hid_t dataspace_id;
    hid_t dataset1_id;
    hid_t dataset2_id;
    hsize_t dim;
    herr_t status;
    char name[32];

    //запись границ домена
    //создали пространство
    dim = n;
    dataspace_id = H5Screate_simple(1, &dim, NULL);

    sprintf(name, "/axis_%s", id);
    //создали набор
    dataset1_id = H5Dcreate2(file_id, name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    //записала данные
    status = H5Dwrite(dataset1_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);

    //закрыли набор
    status = H5Dclose(dataset1_id);

    sprintf(name, "/steps_%s", id);
    //создали набор
    dataset2_id = H5Dcreate2(file_id, name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    //записала данные
    status = H5Dwrite(dataset2_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dx);

    //закрыли набор
    status = H5Dclose(dataset2_id);

    //закрыли пространство
    status = H5Sclose(dataspace_id);
}

void write_nfff_grid()
{
    grid _grd = &fdtd_grd;
    hid_t file_id, dataset_id, dataspace_id;
    hsize_t dims[2];
    herr_t status;
    int dim_data[3];
    double bnd_data[2][3];

    //открытие файла на запись
    file_id = H5Fcreate(NFFF_GRID_FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    //запись размерности домена
    //создали пространство
    dims[0] = 3;
    dataspace_id = H5Screate_simple(1, dims, NULL);

    //создали набор
    dataset_id = H5Dcreate2(file_id, "/dimensions", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    dim_data[0] = nfff_nx;
    dim_data[1] = nfff_ny;
    dim_data[2] = nfff_nz;
    //записала данные
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dim_data);

    //закрыли набор
    status = H5Dclose(dataset_id);

    //закрыли пространство
    status = H5Sclose(dataspace_id);

    //запись границ домена
    //создали пространство
    dims[0] = 2;
    dims[1] = 3;
    dataspace_id = H5Screate_simple(2, dims, NULL);

    //создали набор
    dataset_id = H5Dcreate2(file_id, "/boundaries", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    bnd_data[0][0] = _grd->xi05[nfff_width_xmin + 1];//Xmin
    bnd_data[1][0] = _grd->xi05[_grd->nx - nfff_width_xmax];//Xmax
    bnd_data[0][1] = _grd->yi05[nfff_width_ymin + 1];//Ymin
    bnd_data[1][1] = _grd->yi05[_grd->ny - nfff_width_ymax];//Ymax
    bnd_data[0][2] = _grd->zi05[nfff_width_zmin + 1];//Zmin
    bnd_data[1][2] = _grd->zi05[_grd->nz - nfff_width_zmax];//Zmax
    //записала данные
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, bnd_data);

    //закрыли набор
    status = H5Dclose(dataset_id);

    //закрыли пространство
    status = H5Sclose(dataspace_id);

    //запись осей
    _write_axis(&(_grd->xi[nfff_width_xmin + 1]), &(_grd->dxi[nfff_width_xmin + 1]), dim_data[0], "x", file_id);
    _write_axis(&(_grd->yi[nfff_width_ymin + 1]), &(_grd->dyi[nfff_width_ymin + 1]), dim_data[1], "y", file_id);
    _write_axis(&(_grd->zi[nfff_width_zmin + 1]), &(_grd->dzi[nfff_width_zmin + 1]), dim_data[2], "z", file_id);

    //закрыли файл
    status = H5Fclose(file_id);
}

void write_nfff_domain(int n, double time)
{
    grid _grd = &fdtd_grd;
    hid_t file_id, dataset_id, dataspace_id;
    hsize_t dims[2];
    herr_t status;
    char file_name[256];

    _preload_h_field();

    sprintf(file_name, "%snfff_%d.h5", fdtd_start.result_dir, n);

    //открытие файла на запись
    file_id = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    //Xmin, Xmax
    dims[0] = nfff_ny;
    dims[1] = nfff_nz;

    //создали пространство
    dataspace_id = H5Screate_simple(2, dims, NULL);

    //Xmin
    //создали набор
    dataset_id = H5Dcreate2(file_id, "/ey_x_min", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //записала данные
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nfff_ey_x_min);
    //закрыли набор
    status = H5Dclose(dataset_id);

    //создали набор
    dataset_id = H5Dcreate2(file_id, "/ez_x_min", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //записала данные
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nfff_ez_x_min);
    //закрыли набор
    status = H5Dclose(dataset_id);

    //создали набор
    dataset_id = H5Dcreate2(file_id, "/hy_x_min", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //записала данные
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nfff_hy_x_min);
    //закрыли набор
    status = H5Dclose(dataset_id);

    //создали набор
    dataset_id = H5Dcreate2(file_id, "/hz_x_min", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //записала данные
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nfff_hz_x_min);
    //закрыли набор
    status = H5Dclose(dataset_id);

    //Xmax
    //создали набор
    dataset_id = H5Dcreate2(file_id, "/ey_x_max", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //записала данные
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nfff_ey_x_max);
    //закрыли набор
    status = H5Dclose(dataset_id);

    //создали набор
    dataset_id = H5Dcreate2(file_id, "/ez_x_max", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //записала данные
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nfff_ez_x_max);
    //закрыли набор
    status = H5Dclose(dataset_id);

    //создали набор
    dataset_id = H5Dcreate2(file_id, "/hy_x_max", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //записала данные
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nfff_hy_x_max);
    //закрыли набор
    status = H5Dclose(dataset_id);

    //создали набор
    dataset_id = H5Dcreate2(file_id, "/hz_x_max", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //записала данные
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nfff_hz_x_max);
    //закрыли набор
    status = H5Dclose(dataset_id);

    //закрыли пространство
    status = H5Sclose(dataspace_id);

    //Ymin, Ymax
    dims[0] = nfff_nx;
    dims[1] = nfff_nz;

    //создали пространство
    dataspace_id = H5Screate_simple(2, dims, NULL);

    //Xmin
    //создали набор
    dataset_id = H5Dcreate2(file_id, "/ex_y_min", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //записала данные
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nfff_ex_y_min);
    //закрыли набор
    status = H5Dclose(dataset_id);

    //создали набор
    dataset_id = H5Dcreate2(file_id, "/ez_y_min", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //записала данные
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nfff_ez_y_min);
    //закрыли набор
    status = H5Dclose(dataset_id);

    //создали набор
    dataset_id = H5Dcreate2(file_id, "/hx_y_min", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //записала данные
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nfff_hx_y_min);
    //закрыли набор
    status = H5Dclose(dataset_id);

    //создали набор
    dataset_id = H5Dcreate2(file_id, "/hz_y_min", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //записала данные
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nfff_hz_y_min);
    //закрыли набор
    status = H5Dclose(dataset_id);

    //Xmax
    //создали набор
    dataset_id = H5Dcreate2(file_id, "/ex_y_max", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //записала данные
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nfff_ex_y_max);
    //закрыли набор
    status = H5Dclose(dataset_id);

    //создали набор
    dataset_id = H5Dcreate2(file_id, "/ez_y_max", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //записала данные
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nfff_ez_y_max);
    //закрыли набор
    status = H5Dclose(dataset_id);

    //создали набор
    dataset_id = H5Dcreate2(file_id, "/hx_y_max", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //записала данные
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nfff_hx_y_max);
    //закрыли набор
    status = H5Dclose(dataset_id);

    //создали набор
    dataset_id = H5Dcreate2(file_id, "/hz_y_max", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //записала данные
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nfff_hz_y_max);
    //закрыли набор
    status = H5Dclose(dataset_id);

    //закрыли пространство
    status = H5Sclose(dataspace_id);

    //Zmin, Zmax
    dims[0] = nfff_nx;
    dims[1] = nfff_ny;

    //создали пространство
    dataspace_id = H5Screate_simple(2, dims, NULL);

    //Zmin
    //создали набор
    dataset_id = H5Dcreate2(file_id, "/ex_z_min", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //записала данные
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nfff_ex_z_min);
    //закрыли набор
    status = H5Dclose(dataset_id);

    //создали набор
    dataset_id = H5Dcreate2(file_id, "/ey_z_min", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //записала данные
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nfff_ey_z_min);
    //закрыли набор
    status = H5Dclose(dataset_id);

    //создали набор
    dataset_id = H5Dcreate2(file_id, "/hx_z_min", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //записала данные
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nfff_hx_z_min);
    //закрыли набор
    status = H5Dclose(dataset_id);

    //создали набор
    dataset_id = H5Dcreate2(file_id, "/hy_z_min", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //записала данные
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nfff_hy_z_min);
    //закрыли набор
    status = H5Dclose(dataset_id);

    //Zmax
    //создали набор
    dataset_id = H5Dcreate2(file_id, "/ex_z_max", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //записала данные
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nfff_ex_z_max);
    //закрыли набор
    status = H5Dclose(dataset_id);

    //создали набор
    dataset_id = H5Dcreate2(file_id, "/ey_z_max", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //записала данные
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nfff_ey_z_max);
    //закрыли набор
    status = H5Dclose(dataset_id);

    //создали набор
    dataset_id = H5Dcreate2(file_id, "/hx_z_max", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //записала данные
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nfff_hx_z_max);
    //закрыли набор
    status = H5Dclose(dataset_id);

    //создали набор
    dataset_id = H5Dcreate2(file_id, "/hy_z_max", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //записала данные
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nfff_hy_z_max);
    //закрыли набор
    status = H5Dclose(dataset_id);

    //закрыли пространство
    status = H5Sclose(dataspace_id);

    //закрыли файл
    status = H5Fclose(file_id);

    append_nfff_data(file_name, time);
}
