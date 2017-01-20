
/*! \file 
* \brief Convolutional PML model .
* \author A. Kryukov anton.krv@gmail.com
*/

#include <math.h>
#include "globals.h"
#include "constants.h"
#include "functions.h"
#include "fdtd.h"

pml_double **phu_ex_y_min;/*!< массив y- компоненты свертки для поля Ex. */
pml_double **phu_ex_y_max;/*!< массив y- компоненты свертки для поля Ex. */
pml_double **phu_ex_z_min;/*!< массив z- компоненты свертки для поля Ex. */
pml_double **phu_ex_z_max;/*!< массив z- компоненты свертки для поля Ex. */

pml_double **phu_ey_x_min;/*!< массив x- компоненты свертки для поля Ey. */
pml_double **phu_ey_x_max;/*!< массив x- компоненты свертки для поля Ey. */
pml_double **phu_ey_z_min;/*!< массив z- компоненты свертки для поля Ey. */
pml_double **phu_ey_z_max;/*!< массив z- компоненты свертки для поля Ey. */

pml_double **phu_ez_x_min;/*!< массив x- компоненты свертки для поля Ez. */
pml_double **phu_ez_x_max;/*!< массив x- компоненты свертки для поля Ez. */
pml_double **phu_ez_y_min;/*!< массив y- компоненты свертки для поля Ez. */
pml_double **phu_ez_y_max;/*!< массив y- компоненты свертки для поля Ez. */

pml_double **phu_hx_y_min;/*!< массив y- компоненты свертки для поля Hx. */
pml_double **phu_hx_y_max;/*!< массив y- компоненты свертки для поля Hx. */
pml_double **phu_hx_z_min;/*!< массив z- компоненты свертки для поля Hx. */
pml_double **phu_hx_z_max;/*!< массив z- компоненты свертки для поля Hx. */

pml_double **phu_hy_x_min;/*!< массив x- компоненты свертки для поля Hy. */
pml_double **phu_hy_x_max;/*!< массив x- компоненты свертки для поля Hy. */
pml_double **phu_hy_z_min;/*!< массив z- компоненты свертки для поля Hy. */
pml_double **phu_hy_z_max;/*!< массив z- компоненты свертки для поля Hy. */

pml_double **phu_hz_x_min;/*!< массив x- компоненты свертки для поля Hz. */
pml_double **phu_hz_x_max;/*!< массив x- компоненты свертки для поля Hz. */
pml_double **phu_hz_y_min;/*!< массив y- компоненты свертки для поля Hz. */
pml_double **phu_hz_y_max;/*!< массив y- компоненты свертки для поля Hz. */

pml_double sig_e_x_min;/*!< массив проводимости PML слоя по x для электрического поля. */
pml_double sig_e_y_min;/*!< массив проводимости PML слоя по y для электрического поля. */
pml_double sig_e_z_min;/*!< массив проводимости PML слоя по z для электрического поля. */

pml_double sig_e_x_max;/*!< массив проводимости PML слоя по x для электрического поля. */
pml_double sig_e_y_max;/*!< массив проводимости PML слоя по y для электрического поля. */
pml_double sig_e_z_max;/*!< массив проводимости PML слоя по z для электрического поля. */

pml_double sig_h_x_min;/*!< массив проводимости PML слоя по x для магнитного поля. */
pml_double sig_h_y_min;/*!< массив проводимости PML слоя по y для магнитного поля. */
pml_double sig_h_z_min;/*!< массив проводимости PML слоя по z для магнитного поля. */

pml_double sig_h_x_max;/*!< массив проводимости PML слоя по x для магнитного поля. */
pml_double sig_h_y_max;/*!< массив проводимости PML слоя по y для магнитного поля. */
pml_double sig_h_z_max;/*!< массив проводимости PML слоя по z для магнитного поля. */

pml_double k_x_min;/*!< массив коэффициентов маштабирования PML слоя по x. */
pml_double k_y_min;/*!< массив коэффициентов маштабирования PML слоя по y. */
pml_double k_z_min;/*!< массив коэффициентов маштабирования PML слоя по z. */

pml_double k_x_max;/*!< массив коэффициентов маштабирования PML слоя по x. */
pml_double k_y_max;/*!< массив коэффициентов маштабирования PML слоя по y. */
pml_double k_z_max;/*!< массив коэффициентов маштабирования PML слоя по z. */

pml_double alpha_x_min;/*!< массив коэффициентов ослабления PML слоя по x. */
pml_double alpha_y_min;/*!< массив коэффициентов ослабления PML слоя по y. */
pml_double alpha_z_min;/*!< массив коэффициентов ослабления PML слоя по z. */

pml_double alpha_x_max;/*!< массив коэффициентов ослабления PML слоя по x. */
pml_double alpha_y_max;/*!< массив коэффициентов ослабления PML слоя по y. */
pml_double alpha_z_max;/*!< массив коэффициентов ослабления PML слоя по z. */

void allocatePML()
{
    int i, j, k;
    grid _grd = &fdtd_grd;

    PARMS_malloc(phu_ex_y_min, _grd->nx, pml_double*);
    PARMS_malloc(phu_ex_y_max, _grd->nx, pml_double*);
    PARMS_malloc(phu_ex_z_min, _grd->nx, pml_double*);
    PARMS_malloc(phu_ex_z_max, _grd->nx, pml_double*);
    for (i = 0; i < _grd->nx; ++i) {
        PARMS_malloc(phu_ex_y_min[i], (_grd->nz+1), pml_double);
        PARMS_malloc(phu_ex_y_max[i], (_grd->nz+1), pml_double);

        memset(phu_ex_y_min[i], 0, (_grd->nz+1)*sizeof(pml_double));
        memset(phu_ex_y_max[i], 0, (_grd->nz+1)*sizeof(pml_double));

        PARMS_malloc(phu_ex_z_min[i], (_grd->ny+1), pml_double);
        PARMS_malloc(phu_ex_z_max[i], (_grd->ny+1), pml_double);

        memset(phu_ex_z_min[i], 0, (_grd->ny+1)*sizeof(pml_double));
        memset(phu_ex_z_max[i], 0, (_grd->ny+1)*sizeof(pml_double));
    }

    PARMS_malloc(phu_ey_x_min, _grd->ny, pml_double*);
    PARMS_malloc(phu_ey_x_max, _grd->ny, pml_double*);
    PARMS_malloc(phu_ey_z_min, _grd->ny, pml_double*);
    PARMS_malloc(phu_ey_z_max, _grd->ny, pml_double*);
    for (j = 0; j < _grd->ny; ++j) {
        PARMS_malloc(phu_ey_x_min[j], (_grd->nz+1), pml_double);
        PARMS_malloc(phu_ey_x_max[j], (_grd->nz+1), pml_double);

        memset(phu_ey_x_min[j], 0, (_grd->nz+1)*sizeof(pml_double));
        memset(phu_ey_x_max[j], 0, (_grd->nz+1)*sizeof(pml_double));

        PARMS_malloc(phu_ey_z_min[j], (_grd->nx+1), pml_double);
        PARMS_malloc(phu_ey_z_max[j], (_grd->nx+1), pml_double);

        memset(phu_ey_z_min[j], 0, (_grd->nx+1)*sizeof(pml_double));
        memset(phu_ey_z_max[j], 0, (_grd->nx+1)*sizeof(pml_double));
    }

    PARMS_malloc(phu_ez_x_min, _grd->nz, pml_double*);
    PARMS_malloc(phu_ez_x_max, _grd->nz, pml_double*);
    PARMS_malloc(phu_ez_y_min, _grd->nz, pml_double*);
    PARMS_malloc(phu_ez_y_max, _grd->nz, pml_double*);
    for (k = 0; k < _grd->nz; ++k) {
        PARMS_malloc(phu_ez_x_min[k], (_grd->ny+1), pml_double);
        PARMS_malloc(phu_ez_x_max[k], (_grd->ny+1), pml_double);

        memset(phu_ez_x_min[k], 0, (_grd->ny+1)*sizeof(pml_double));
        memset(phu_ez_x_max[k], 0, (_grd->ny+1)*sizeof(pml_double));

        PARMS_malloc(phu_ez_y_min[k], (_grd->nx+1), pml_double);
        PARMS_malloc(phu_ez_y_max[k], (_grd->nx+1), pml_double);

        memset(phu_ez_y_min[k], 0, (_grd->nx+1)*sizeof(pml_double));
        memset(phu_ez_y_max[k], 0, (_grd->nx+1)*sizeof(pml_double));
    }

    PARMS_malloc(phu_hx_y_min, (_grd->nx+1), pml_double*);
    PARMS_malloc(phu_hx_y_max, (_grd->nx+1), pml_double*);
    PARMS_malloc(phu_hx_z_min, (_grd->nx+1), pml_double*);
    PARMS_malloc(phu_hx_z_max, (_grd->nx+1), pml_double*);
    for (i = 0; i < (_grd->nx+1); ++i) {
        PARMS_malloc(phu_hx_y_min[i], (_grd->nz), pml_double);
        PARMS_malloc(phu_hx_y_max[i], (_grd->nz), pml_double);

        memset(phu_hx_y_min[i], 0, (_grd->nz)*sizeof(pml_double));
        memset(phu_hx_y_max[i], 0, (_grd->nz)*sizeof(pml_double));

        PARMS_malloc(phu_hx_z_min[i], (_grd->ny), pml_double);
        PARMS_malloc(phu_hx_z_max[i], (_grd->ny), pml_double);

        memset(phu_hx_z_min[i], 0, (_grd->ny)*sizeof(pml_double));
        memset(phu_hx_z_max[i], 0, (_grd->ny)*sizeof(pml_double));
    }

    PARMS_malloc(phu_hy_x_min, (_grd->ny+1), pml_double*);
    PARMS_malloc(phu_hy_x_max, (_grd->ny+1), pml_double*);
    PARMS_malloc(phu_hy_z_min, (_grd->ny+1), pml_double*);
    PARMS_malloc(phu_hy_z_max, (_grd->ny+1), pml_double*);
    for (j = 0; j < (_grd->ny+1); ++j) {
        PARMS_malloc(phu_hy_x_min[j], (_grd->nz), pml_double);
        PARMS_malloc(phu_hy_x_max[j], (_grd->nz), pml_double);

        memset(phu_hy_x_min[j], 0, (_grd->nz)*sizeof(pml_double));
        memset(phu_hy_x_max[j], 0, (_grd->nz)*sizeof(pml_double));

        PARMS_malloc(phu_hy_z_min[j], (_grd->nx), pml_double);
        PARMS_malloc(phu_hy_z_max[j], (_grd->nx), pml_double);

        memset(phu_hy_z_min[j], 0, (_grd->nx)*sizeof(pml_double));
        memset(phu_hy_z_max[j], 0, (_grd->nx)*sizeof(pml_double));
    }

    PARMS_malloc(phu_hz_x_min, (_grd->nz+1), pml_double*);
    PARMS_malloc(phu_hz_x_max, (_grd->nz+1), pml_double*);
    PARMS_malloc(phu_hz_y_min, (_grd->nz+1), pml_double*);
    PARMS_malloc(phu_hz_y_max, (_grd->nz+1), pml_double*);
    for (k = 0; k < (_grd->nz+1); ++k) {
        PARMS_malloc(phu_hz_x_min[k], (_grd->ny), pml_double);
        PARMS_malloc(phu_hz_x_max[k], (_grd->ny), pml_double);

        memset(phu_hz_x_min[k], 0, (_grd->ny)*sizeof(pml_double));
        memset(phu_hz_x_max[k], 0, (_grd->ny)*sizeof(pml_double));

        PARMS_malloc(phu_hz_y_min[k], (_grd->nx), pml_double);
        PARMS_malloc(phu_hz_y_max[k], (_grd->nx), pml_double);

        memset(phu_hz_y_min[k], 0, (_grd->nx)*sizeof(pml_double));
        memset(phu_hz_y_max[k], 0, (_grd->nx)*sizeof(pml_double));
    }
}

void set_PML_profile()
{
    int i;
    grid _grd = &fdtd_grd;
    double sigma_max;
    double m = 4.0;//порядок полимиального роста проводиомсти PML слоя
    double d;//толщина слоя в см.
    double lnR = 6.0;//-логарифм коэффициента ослабления
    double eta = 1.0;//волновое сопротивление PML слоя (импенданс)
    double kmax = 1.01;//параметр геометрического роста разностного шага
    double factor;

    for(i = 0; i < PML_WIDTH; ++i) {
        //xmin
        sig_e_x_min[i] = 0.0;
        sig_h_x_min[i] = 0.0;
        k_x_min[i] = 1.0;
        alpha_x_min[i] = 0.0;
        //xmax
        sig_e_x_max[i] = 0.0;
        sig_h_x_max[i] = 0.0;
        k_x_max[i] = 1.0;
        alpha_x_max[i] = 0.0;

        //ymin
        sig_e_y_min[i] = 0.0;
        sig_h_y_min[i] = 0.0;
        k_y_min[i] = 1.0;
        alpha_y_min[i] = 0.0;
        //ymax
        sig_e_y_max[i] = 0.0;
        sig_h_y_max[i] = 0.0;
        k_y_max[i] = 1.0;
        alpha_y_max[i] = 0.0;

        //zmin
        sig_e_z_min[i] = 0.0;
        sig_h_z_min[i] = 0.0;
        k_z_min[i] = 1.0;
        alpha_z_min[i] = 0.0;
        //zmax
        sig_e_z_max[i] = 0.0;
        sig_h_z_max[i] = 0.0;
        k_z_max[i] = 1.0;
        alpha_z_max[i] = 0.0;
    }

    //эти формулы выведены в СИ!

    //профиль для xmin
    d = _grd->xi[PML_WIDTH] - _grd->xi[0];

    sigma_max = (_cl*(m + 1.0)*lnR)/(2.0*eta*d);

    for (i = 0; i < PML_WIDTH; ++i) {
        factor = pow(((_grd->xi[PML_WIDTH] - _grd->xi[i]) / d), m);
        sig_e_x_min[i] = factor*sigma_max;
        k_x_min[i] = 1.0 + (kmax - 1.0)*factor;
    }

    for (i = 0; i < PML_WIDTH; ++i) {
        factor = pow(((_grd->xi[PML_WIDTH] - _grd->xi05[i + 1]) / d), m);
        sig_h_x_min[i] = factor*sigma_max;
    }

    //профиль для xmax
    d = _grd->xi[_grd->nx] - _grd->xi[_grd->nx - PML_WIDTH];

    sigma_max = (_cl*(m + 1.0)*lnR)/(2.0*eta*d);

    for (i = 1; i <= PML_WIDTH; ++i) {
        factor = pow(((_grd->xi[i + _grd->nx - PML_WIDTH] - _grd->xi[_grd->nx - PML_WIDTH]) / d), m);
        sig_e_x_max[i-1] = factor*sigma_max;
        k_x_max[i-1] = 1.0 + (kmax - 1.0)*factor;
    }

    for (i = 0; i < PML_WIDTH; ++i) {
        factor = pow(((_grd->xi05[i + 1 + _grd->nx - PML_WIDTH] - _grd->xi[_grd->nx - PML_WIDTH]) / d), m);
        sig_h_x_max[i] = factor*sigma_max;
    }

    //профиль для ymin
    d = _grd->yi[PML_WIDTH] - _grd->yi[0];

    sigma_max = (_cl*(m + 1.0)*lnR)/(2.0*eta*d);

    for (i = 0; i < PML_WIDTH; ++i) {
        factor = pow(((_grd->yi[PML_WIDTH] - _grd->yi[i]) / d), m);
        sig_e_y_min[i] = factor*sigma_max;
        k_y_min[i] = 1.0 + (kmax - 1.0)*factor;
    }

    for (i = 0; i < PML_WIDTH; ++i) {
        factor = pow(((_grd->yi[PML_WIDTH] - _grd->yi05[i + 1]) / d), m);
        sig_h_y_min[i] = factor*sigma_max;
    }

    //профиль для ymax
    d = _grd->yi[_grd->ny] - _grd->yi[_grd->ny - PML_WIDTH];

    sigma_max = (_cl*(m + 1.0)*lnR)/(2.0*eta*d);

    for (i = 1; i <= PML_WIDTH; ++i) {
        factor = pow(((_grd->yi[i + _grd->ny - PML_WIDTH] - _grd->yi[_grd->ny - PML_WIDTH]) / d), m);
        sig_e_y_max[i-1] = factor*sigma_max;
        k_y_max[i-1] = 1.0 + (kmax - 1.0)*factor;
    }

    for (i = 0; i < PML_WIDTH; ++i) {
        factor = pow(((_grd->yi05[i + 1 + _grd->ny - PML_WIDTH] - _grd->yi[_grd->ny - PML_WIDTH]) / d), m);
        sig_h_y_max[i] = factor*sigma_max;
    }

    //профиль для zmin
    d = _grd->zi[PML_WIDTH] - _grd->zi[0];

    sigma_max = (_cl*(m + 1.0)*lnR)/(2.0*eta*d);

    for (i = 0; i < PML_WIDTH; ++i) {
        factor = pow(((_grd->zi[PML_WIDTH] - _grd->zi[i]) / d), m);
        sig_e_z_min[i] = factor*sigma_max;
        k_z_min[i] = 1.0 + (kmax - 1.0)*factor;
    }

    for (i = 0; i < PML_WIDTH; ++i) {
        factor = pow(((_grd->zi[PML_WIDTH] - _grd->zi05[i + 1]) / d), m);
        sig_h_z_min[i] = factor*sigma_max;
    }

    //профиль для zmax
    d = _grd->zi[_grd->nz] - _grd->zi[_grd->nz - PML_WIDTH];

    sigma_max = (_cl*(m + 1.0)*lnR)/(2.0*eta*d);

    for (i = 1; i <= PML_WIDTH; ++i) {
        factor = pow(((_grd->zi[i + _grd->nz - PML_WIDTH] - _grd->zi[_grd->nz - PML_WIDTH]) / d), m);
        sig_e_z_max[i-1] = factor*sigma_max;
        k_z_max[i-1] = 1.0 + (kmax - 1.0)*factor;
    }

    for (i = 0; i < PML_WIDTH; ++i) {
        factor = pow(((_grd->zi05[i + 1 + _grd->nz - PML_WIDTH] - _grd->zi[_grd->nz - PML_WIDTH]) / d), m);
        sig_h_z_max[i] = factor*sigma_max;
    }
}

double alpha_i_e_min(int i, direction dir)
{
    double val = 0.0;
    double ht = fdtd_ht;
    double *a_i = NULL;
    double *k_i = NULL;
    double *sig_i = NULL;

    switch (dir) {
    case X_DIR :
        a_i = alpha_x_min;
        k_i = k_x_min;
        sig_i = sig_e_x_min;
        break;

    case Y_DIR :
        a_i = alpha_y_min;
        k_i = k_y_min;
        sig_i = sig_e_y_min;
        break;

    case Z_DIR :
        a_i = alpha_z_min;
        k_i = k_z_min;
        sig_i = sig_e_z_min;
        break;
    }

    CHECKERR(a_i == NULL);
    CHECKERR(k_i == NULL);
    CHECKERR(sig_i == NULL);

    if (sig_i[i] == 0.0)
        return 0.0;

    val = (sig_i[i]/(sig_i[i]*k_i[i] + k_i[i]*k_i[i]*a_i[i]))*(exp(-(sig_i[i]/k_i[i] + a_i[i])*ht) - 1.0);

    return val;
}

double alpha_i_e_max(int i, direction dir)
{
    double val = 0.0;
    double ht = fdtd_ht;
    double *a_i = NULL;
    double *k_i = NULL;
    double *sig_i = NULL;

    switch (dir) {
    case X_DIR :
        a_i = alpha_x_max;
        k_i = k_x_max;
        sig_i = sig_e_x_max;
        break;

    case Y_DIR :
        a_i = alpha_y_max;
        k_i = k_y_max;
        sig_i = sig_e_y_max;
        break;

    case Z_DIR :
        a_i = alpha_z_max;
        k_i = k_z_max;
        sig_i = sig_e_z_max;
        break;
    }

    CHECKERR(a_i == NULL);
    CHECKERR(k_i == NULL);
    CHECKERR(sig_i == NULL);

    if (sig_i[i] == 0.0)
        return 0.0;

    val = (sig_i[i]/(sig_i[i]*k_i[i] + k_i[i]*k_i[i]*a_i[i]))*(exp(-(sig_i[i]/k_i[i] + a_i[i])*ht) - 1.0);

    return val;
}

double alpha_i_h_min(int i, direction dir)
{
    double val = 0.0;
    double ht = fdtd_ht05;
    double *a_i = NULL;
    double *k_i = NULL;
    double *sig_i = NULL;

    switch (dir) {
    case X_DIR :
        a_i = alpha_x_min;
        k_i = k_x_min;
        sig_i = sig_h_x_min;
        break;

    case Y_DIR :
        a_i = alpha_y_min;
        k_i = k_y_min;
        sig_i = sig_h_y_min;
        break;

    case Z_DIR :
        a_i = alpha_z_min;
        k_i = k_z_min;
        sig_i = sig_h_z_min;
        break;
    }

    CHECKERR(a_i == NULL);
    CHECKERR(k_i == NULL);
    CHECKERR(sig_i == NULL);

    if (sig_i[i] == 0.0)
        return 0.0;

    val = (sig_i[i]/(sig_i[i]*k_i[i] + k_i[i]*k_i[i]*a_i[i]))*(exp(-(sig_i[i]/k_i[i] + a_i[i])*ht) - 1.0);

    return val;
}

double alpha_i_h_max(int i, direction dir)
{
    double val = 0.0;
    double ht = fdtd_ht05;
    double *a_i = NULL;
    double *k_i = NULL;
    double *sig_i = NULL;

    switch (dir) {
    case X_DIR :
        a_i = alpha_x_max;
        k_i = k_x_max;
        sig_i = sig_h_x_max;
        break;

    case Y_DIR :
        a_i = alpha_y_max;
        k_i = k_y_max;
        sig_i = sig_h_y_max;
        break;

    case Z_DIR :
        a_i = alpha_z_max;
        k_i = k_z_max;
        sig_i = sig_h_z_max;
        break;
    }

    CHECKERR(a_i == NULL);
    CHECKERR(k_i == NULL);
    CHECKERR(sig_i == NULL);

    if (sig_i[i] == 0.0)
        return 0.0;

    val = (sig_i[i]/(sig_i[i]*k_i[i] + k_i[i]*k_i[i]*a_i[i]))*(exp(-(sig_i[i]/k_i[i] + a_i[i])*ht) - 1.0);

    return val;
}

double betta_i_e_min(int i, direction dir)
{
    double val = 0.0;
    double ht = fdtd_ht;
    double *a_i = NULL;
    double *k_i = NULL;
    double *sig_i = NULL;

    switch (dir) {
    case X_DIR :
        a_i = alpha_x_min;
        k_i = k_x_min;
        sig_i = sig_e_x_min;
        break;

    case Y_DIR :
        a_i = alpha_y_min;
        k_i = k_y_min;
        sig_i = sig_e_y_min;
        break;

    case Z_DIR :
        a_i = alpha_z_min;
        k_i = k_z_min;
        sig_i = sig_e_z_min;
        break;
    }

    CHECKERR(a_i == NULL);
    CHECKERR(k_i == NULL);
    CHECKERR(sig_i == NULL);

    val = exp(-(sig_i[i]/k_i[i] + a_i[i])*ht);

    return val;
}

double betta_i_e_max(int i, direction dir)
{
    double val = 0.0;
    double ht = fdtd_ht;
    double *a_i = NULL;
    double *k_i = NULL;
    double *sig_i = NULL;

    switch (dir) {
    case X_DIR :
        a_i = alpha_x_max;
        k_i = k_x_max;
        sig_i = sig_e_x_max;
        break;

    case Y_DIR :
        a_i = alpha_y_max;
        k_i = k_y_max;
        sig_i = sig_e_y_max;
        break;

    case Z_DIR :
        a_i = alpha_z_max;
        k_i = k_z_max;
        sig_i = sig_e_z_max;
        break;
    }

    CHECKERR(a_i == NULL);
    CHECKERR(k_i == NULL);
    CHECKERR(sig_i == NULL);

    val = exp(-(sig_i[i]/k_i[i] + a_i[i])*ht);

    return val;
}

double betta_i_h_min(int i, direction dir)
{
    double val = 0.0;
    double ht = fdtd_ht05;
    double *a_i = NULL;
    double *k_i = NULL;
    double *sig_i = NULL;

    switch (dir) {
    case X_DIR :
        a_i = alpha_x_min;
        k_i = k_x_min;
        sig_i = sig_h_x_min;
        break;

    case Y_DIR :
        a_i = alpha_y_min;
        k_i = k_y_min;
        sig_i = sig_h_y_min;
        break;

    case Z_DIR :
        a_i = alpha_z_min;
        k_i = k_z_min;
        sig_i = sig_h_z_min;
        break;
    }

    CHECKERR(a_i == NULL);
    CHECKERR(k_i == NULL);
    CHECKERR(sig_i == NULL);

    val = exp(-(sig_i[i]/k_i[i] + a_i[i])*ht);

    return val;
}

double betta_i_h_max(int i, direction dir)
{
    double val = 0.0;
    double ht = fdtd_ht05;
    double *a_i = NULL;
    double *k_i = NULL;
    double *sig_i = NULL;

    switch (dir) {
    case X_DIR :
        a_i = alpha_x_max;
        k_i = k_x_max;
        sig_i = sig_h_x_max;
        break;

    case Y_DIR :
        a_i = alpha_y_max;
        k_i = k_y_max;
        sig_i = sig_h_y_max;
        break;

    case Z_DIR :
        a_i = alpha_z_max;
        k_i = k_z_max;
        sig_i = sig_h_z_max;
        break;
    }

    CHECKERR(a_i == NULL);
    CHECKERR(k_i == NULL);
    CHECKERR(sig_i == NULL);

    val = exp(-(sig_i[i]/k_i[i] + a_i[i])*ht);

    return val;
}
