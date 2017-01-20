
/*! \file
* \brief Функции расчета полей с MUR на границе.
* \author Крюков А.А. anton.krv@gmail.com
*/

#include "../../include/globals.h"
#include "../../include/constants.h"
#include "../../include/functions.h"
#include "../../include/fdtd.h"

mur_double **mur_ex_y_min;/*!< массив y- компоненты для поля Ex. */
mur_double **mur_ex_y_max;/*!< массив y- компоненты для поля Ex. */
mur_double **mur_ex_z_min;/*!< массив z- компоненты для поля Ex. */
mur_double **mur_ex_z_max;/*!< массив z- компоненты для поля Ex. */

mur_double **mur_ey_x_min;/*!< массив x- компоненты для поля Ey. */
mur_double **mur_ey_x_max;/*!< массив x- компоненты для поля Ey. */
mur_double **mur_ey_z_min;/*!< массив z- компоненты для поля Ey. */
mur_double **mur_ey_z_max;/*!< массив z- компоненты для поля Ey. */

mur_double **mur_ez_x_min;/*!< массив x- компоненты для поля Ez. */
mur_double **mur_ez_x_max;/*!< массив x- компоненты для поля Ez. */
mur_double **mur_ez_y_min;/*!< массив y- компоненты для поля Ez. */
mur_double **mur_ez_y_max;/*!< массив y- компоненты для поля Ez. */

void allocateMUR()
{
    int i, j, k;
    grid _grd = &fdtd_grd;

    PARMS_malloc(mur_ex_y_min, _grd->nx, mur_double*);
    PARMS_malloc(mur_ex_y_max, _grd->nx, mur_double*);
    PARMS_malloc(mur_ex_z_min, _grd->nx, mur_double*);
    PARMS_malloc(mur_ex_z_max, _grd->nx, mur_double*);
    for (i = 0; i < _grd->nx; ++i) {
        PARMS_malloc(mur_ex_y_min[i], (_grd->nz + 1), mur_double);
        PARMS_malloc(mur_ex_y_max[i], (_grd->nz + 1), mur_double);

        memset(mur_ex_y_min[i], 0, (_grd->nz + 1) * sizeof(mur_double));
        memset(mur_ex_y_max[i], 0, (_grd->nz + 1) * sizeof(mur_double));

        PARMS_malloc(mur_ex_z_min[i], (_grd->ny + 1), mur_double);
        PARMS_malloc(mur_ex_z_max[i], (_grd->ny + 1), mur_double);

        memset(mur_ex_z_min[i], 0, (_grd->ny + 1) * sizeof(mur_double));
        memset(mur_ex_z_max[i], 0, (_grd->ny + 1) * sizeof(mur_double));
    }

    PARMS_malloc(mur_ey_x_min, _grd->ny, mur_double*);
    PARMS_malloc(mur_ey_x_max, _grd->ny, mur_double*);
    PARMS_malloc(mur_ey_z_min, _grd->ny, mur_double*);
    PARMS_malloc(mur_ey_z_max, _grd->ny, mur_double*);
    for (j = 0; j < _grd->ny; ++j) {
        PARMS_malloc(mur_ey_x_min[j], (_grd->nz + 1), mur_double);
        PARMS_malloc(mur_ey_x_max[j], (_grd->nz + 1), mur_double);

        memset(mur_ey_x_min[j], 0, (_grd->nz + 1) * sizeof(mur_double));
        memset(mur_ey_x_max[j], 0, (_grd->nz + 1) * sizeof(mur_double));

        PARMS_malloc(mur_ey_z_min[j], (_grd->nx + 1), mur_double);
        PARMS_malloc(mur_ey_z_max[j], (_grd->nx + 1), mur_double);

        memset(mur_ey_z_min[j], 0, (_grd->nx + 1) * sizeof(mur_double));
        memset(mur_ey_z_max[j], 0, (_grd->nx + 1) * sizeof(mur_double));
    }

    PARMS_malloc(mur_ez_x_min, _grd->nz, mur_double*);
    PARMS_malloc(mur_ez_x_max, _grd->nz, mur_double*);
    PARMS_malloc(mur_ez_y_min, _grd->nz, mur_double*);
    PARMS_malloc(mur_ez_y_max, _grd->nz, mur_double*);
    for (k = 0; k < _grd->nz; ++k) {
        PARMS_malloc(mur_ez_x_min[k], (_grd->ny + 1), mur_double);
        PARMS_malloc(mur_ez_x_max[k], (_grd->ny + 1), mur_double);

        memset(mur_ez_x_min[k], 0, (_grd->ny + 1) * sizeof(mur_double));
        memset(mur_ez_x_max[k], 0, (_grd->ny + 1) * sizeof(mur_double));

        PARMS_malloc(mur_ez_y_min[k], (_grd->nx + 1), mur_double);
        PARMS_malloc(mur_ez_y_max[k], (_grd->nx + 1), mur_double);

        memset(mur_ez_y_min[k], 0, (_grd->nx + 1) * sizeof(mur_double));
        memset(mur_ez_y_max[k], 0, (_grd->nx + 1) * sizeof(mur_double));
    }
}
