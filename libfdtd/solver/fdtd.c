
/*! \file 
* \brief Функции расчета полей.
* \author A. Kryukov. anton.krv@gmail.com
*/

#include <math.h>
#include <omp.h>
#include "globals.h"
#include "constants.h"
#include "functions.h"
#include "fdtd.h"

//вычисление Ex в точке
static void calc_Ex_point(int i, int j, int k)
{
    double ht = fdtd_ht;
    grid _grd = &fdtd_grd;
    double dhy, dhz;
    double sig_tmp = _pi4*fun_sig(i, j, k, EX_ID);
    double eps_tmp = fun_eps(i, j, k, EX_ID);
    double ky, ay, by, phu_ex_y;
    double kz, az, bz, phu_ex_z;
    int ymin, ymax, j_pml;
    int zmin, zmax, k_pml;

    dhy = (fdtd_Hy[j][k][i-1] - fdtd_Hy[j][k-1][i-1])/_grd->dzi[k];
    dhz = (fdtd_Hz[k][i-1][j] - fdtd_Hz[k][i-1][j-1])/_grd->dyi[j];

    if (fdtd_comp.l_gran == BND_CPML) {
        ymin = (j < PML_WIDTH);
        ymax = (j > _grd->ny - PML_WIDTH);

        zmin = (k < PML_WIDTH);
        zmax = (k > _grd->nz - PML_WIDTH);

        if (!ymin && !ymax &&
            !zmin && !zmax) {
            if (fdtd_comp.l_disp && (i >= 10 && i <= 11) && (j >= 9 && j <= 11) && (k >= 9 && k <= 190)) {
                calc_Ex_point_debye(i, j, k);
            } else {
                fdtd_Ex[i-1][j][k] = ((eps_tmp - 0.5*ht*sig_tmp)*fdtd_Ex[i-1][j][k]
                            + ht*_cl*(dhz - dhy) - ht*_pi4*fdtd_Jx[i-1][j][k])/(eps_tmp + 0.5*ht*sig_tmp);
                }
        } else {
            if (ymin) {
                ky = k_y_min[j];
                ay = alpha_i_e_min(j, Y_DIR);
                by = betta_i_e_min(j, Y_DIR);
                phu_ex_y = phu_ex_y_min[i-1][k][j] = by*phu_ex_y_min[i-1][k][j] + ay*dhz;
            } else if (ymax) {
                j_pml = j - _grd->ny + PML_WIDTH;
                ky = k_y_max[j_pml-1];
                ay = alpha_i_e_max(j_pml-1, Y_DIR);
                by = betta_i_e_max(j_pml-1, Y_DIR);
                phu_ex_y = phu_ex_y_max[i-1][k][j_pml] = by*phu_ex_y_max[i-1][k][j_pml] + ay*dhz;
            } else {
                ky = 1.0;
                ay = 0.0;
                by = 1.0;
                phu_ex_y = 0.0;
            }

            if (zmin) {
                kz = k_z_min[k];
                az = alpha_i_e_min(k, Z_DIR);
                bz = betta_i_e_min(k, Z_DIR);
                phu_ex_z = phu_ex_z_min[i-1][j][k] = bz*phu_ex_z_min[i-1][j][k] + az*dhy;
            } else if (zmax) {
                k_pml = k - _grd->nz + PML_WIDTH;
                kz = k_z_max[k_pml-1];
                az = alpha_i_e_max(k_pml-1, Z_DIR);
                bz = betta_i_e_max(k_pml-1, Z_DIR);
                phu_ex_z = phu_ex_z_max[i-1][j][k_pml] = bz*phu_ex_z_max[i-1][j][k_pml] + az*dhy;
            } else {
                kz = 1.0;
                az = 0.0;
                bz = 1.0;
                phu_ex_z = 0.0;
            }

            fdtd_Ex[i-1][j][k] = ((eps_tmp - 0.5*ht*sig_tmp)*fdtd_Ex[i-1][j][k]
                            + ht*_cl*((dhz/ky - dhy/kz) + (phu_ex_y - phu_ex_z))
                            - ht*_pi4*fdtd_Jx[i-1][j][k])/(eps_tmp + 0.5*ht*sig_tmp);
        }
    } else {
        fdtd_Ex[i-1][j][k] = ((eps_tmp - 0.5*ht*sig_tmp)*fdtd_Ex[i-1][j][k]
                    + ht*_cl*(dhz - dhy) - ht*_pi4*fdtd_Jx[i-1][j][k])/(eps_tmp + 0.5*ht*sig_tmp);
    }
}

//вычисление Ey в точке
static void calc_Ey_point(int i, int j, int k)
{
    double ht = fdtd_ht;
    grid _grd = &fdtd_grd;
    double dhx, dhz;
    double sig_tmp = _pi4*fun_sig(i, j, k, EY_ID);
    double eps_tmp = fun_eps(i, j, k, EY_ID);
    double kx, ax, bx, phu_ey_x;
    double kz, az, bz, phu_ey_z;
    int xmin, xmax, i_pml;
    int zmin, zmax, k_pml;

    dhx = (fdtd_Hx[i][j-1][k] - fdtd_Hx[i][j-1][k-1])/_grd->dzi[k];
    dhz = (fdtd_Hz[k][i][j-1] - fdtd_Hz[k][i-1][j-1])/_grd->dxi[i];

    if (fdtd_comp.l_gran == BND_CPML) {
        xmin = (i < PML_WIDTH);
        xmax = (i > _grd->nx - PML_WIDTH);

        zmin = (k < PML_WIDTH);
        zmax = (k > _grd->nz - PML_WIDTH);

        if (!xmin && !xmax &&
            !zmin && !zmax) {
            if (fdtd_comp.l_disp && (i >= 9 && i <= 11) && (j >= 10 && j <= 11) && (k >= 9 && k <= 190)) {
                calc_Ey_point_debye(i, j, k);
            } else {
                fdtd_Ey[j-1][k][i] = ((eps_tmp - 0.5*ht*sig_tmp)*fdtd_Ey[j-1][k][i]
                            + ht*_cl*(dhx - dhz) - ht*_pi4*fdtd_Jy[j-1][k][i])/(eps_tmp + 0.5*ht*sig_tmp);
            }
        } else {
            if (xmin) {
                kx = k_x_min[i];
                ax = alpha_i_e_min(i, X_DIR);
                bx = betta_i_e_min(i, X_DIR);
                phu_ey_x = phu_ey_x_min[j-1][k][i] = bx*phu_ey_x_min[j-1][k][i] + ax*dhz;
            } else if (xmax) {
                i_pml = i - _grd->nx + PML_WIDTH;
                kx = k_x_max[i_pml-1];
                ax = alpha_i_e_max(i_pml-1, X_DIR);
                bx = betta_i_e_max(i_pml-1, X_DIR);
                phu_ey_x = phu_ey_x_max[j-1][k][i_pml] = bx*phu_ey_x_max[j-1][k][i_pml] + ax*dhz;
            } else {
                kx = 1.0;
                ax = 0.0;
                bx = 1.0;
                phu_ey_x = 0.0;
            }

            if (zmin) {
                kz = k_z_min[k];
                az = alpha_i_e_min(k, Z_DIR);
                bz = betta_i_e_min(k, Z_DIR);
                phu_ey_z = phu_ey_z_min[j-1][i][k] = bz*phu_ey_z_min[j-1][i][k] + az*dhx;
            } else if(zmax) {
                k_pml = k - _grd->nz + PML_WIDTH;
                kz = k_z_max[k_pml-1];
                az = alpha_i_e_max(k_pml-1, Z_DIR);
                bz = betta_i_e_max(k_pml-1, Z_DIR);
                phu_ey_z = phu_ey_z_max[j-1][i][k_pml] = bz*phu_ey_z_max[j-1][i][k_pml] + az*dhx;
            } else {
                kz = 1.0;
                az = 0.0;
                bz = 1.0;
                phu_ey_z = 0.0;
            }

            fdtd_Ey[j-1][k][i] = ((eps_tmp - 0.5*ht*sig_tmp)*fdtd_Ey[j-1][k][i]
                            + ht*_cl*((dhx/kz - dhz/kx) + (phu_ey_z - phu_ey_x))
                            - ht*_pi4*fdtd_Jy[j-1][k][i])/(eps_tmp + 0.5*ht*sig_tmp);
        }
    } else {
        fdtd_Ey[j-1][k][i] = ((eps_tmp - 0.5*ht*sig_tmp)*fdtd_Ey[j-1][k][i]
                        + ht*_cl*(dhx - dhz) - ht*_pi4*fdtd_Jy[j-1][k][i])/(eps_tmp + 0.5*ht*sig_tmp);
    }
}

//вычисление Ez в точке
static void calc_Ez_point(int i, int j, int k)
{
    double ht = fdtd_ht;
    grid _grd = &fdtd_grd;
    double dhx, dhy;
    double sig_tmp = _pi4*fun_sig(i, j, k, EZ_ID);
    double eps_tmp = fun_eps(i, j, k, EZ_ID);
    double kx, ax, bx, phu_ez_x;
    double ky, ay, by, phu_ez_y;
    int xmin, xmax, i_pml;
    int ymin, ymax, j_pml;

    dhx = (fdtd_Hx[i][j][k-1] - fdtd_Hx[i][j-1][k-1])/_grd->dyi[j];
    dhy = (fdtd_Hy[j][k-1][i] - fdtd_Hy[j][k-1][i-1])/_grd->dxi[i];

    if (fdtd_comp.l_gran == BND_CPML) {
        xmin = (i < PML_WIDTH);
        xmax = (i > _grd->nx - PML_WIDTH);

        ymin = (j < PML_WIDTH);
        ymax = (j > _grd->ny - PML_WIDTH);

        if (!xmin && !xmax &&
            !ymin && !ymax) {
            if (fdtd_comp.l_disp && (i >= 9 && i <= 11) && (j >= 9 && j <= 11) && (k >= 10 && k <= 190)) {
                calc_Ez_point_debye(i, j, k);
            } else {
                fdtd_Ez[k-1][i][j] = ((eps_tmp - 0.5*ht*sig_tmp)*fdtd_Ez[k-1][i][j]
                            + ht*_cl*(dhy - dhx) - ht*_pi4*fdtd_Jz[k-1][i][j])/(eps_tmp + 0.5*ht*sig_tmp);
            }
        } else {
            if (xmin) {
                kx = k_x_min[i];
                ax = alpha_i_e_min(i, X_DIR);
                bx = betta_i_e_min(i, X_DIR);
                phu_ez_x = phu_ez_x_min[k-1][j][i] = bx*phu_ez_x_min[k-1][j][i] + ax*dhy;
            } else if (xmax) {
                i_pml = i - _grd->nx + PML_WIDTH;
                kx = k_x_max[i_pml-1];
                ax = alpha_i_e_max(i_pml-1, X_DIR);
                bx = betta_i_e_max(i_pml-1, X_DIR);
                phu_ez_x = phu_ez_x_max[k-1][j][i_pml] = bx*phu_ez_x_max[k-1][j][i_pml] + ax*dhy;
            } else {
                kx = 1.0;
                ax = 0.0;
                bx = 1.0;
                phu_ez_x = 0.0;
            }

            if (ymin) {
                ky = k_y_min[j];
                ay = alpha_i_e_min(j, Y_DIR);
                by = betta_i_e_min(j, Y_DIR);
                phu_ez_y = phu_ez_y_min[k-1][i][j] = by*phu_ez_y_min[k-1][i][j] + ay*dhx;
            } else if (ymax) {
                j_pml = j - _grd->ny + PML_WIDTH;
                ky = k_y_max[j_pml-1];
                ay = alpha_i_e_max(j_pml-1, Y_DIR);
                by = betta_i_e_max(j_pml-1, Y_DIR);
                phu_ez_y = phu_ez_y_max[k-1][i][j_pml] = by*phu_ez_y_max[k-1][i][j_pml] + ay*dhx;
            } else {
                ky = 1.0;
                ay = 0.0;
                by = 1.0;
                phu_ez_y = 0.0;
            }

            fdtd_Ez[k-1][i][j] = ((eps_tmp - 0.5*ht*sig_tmp)*fdtd_Ez[k-1][i][j]
                            + ht*_cl*((dhy/kx - dhx/ky) + (phu_ez_x - phu_ez_y))
                            - ht*_pi4*fdtd_Jz[k-1][i][j])/(eps_tmp + 0.5*ht*sig_tmp);
        }
    } else {
        fdtd_Ez[k-1][i][j] = ((eps_tmp - 0.5*ht*sig_tmp)*fdtd_Ez[k-1][i][j]
                        + ht*_cl*(dhy - dhx) - ht*_pi4*fdtd_Jz[k-1][i][j])/(eps_tmp + 0.5*ht*sig_tmp);
    }
}

//вычисление Hx в точке
static void calc_Hx_point(int i, int j, int k)
{
    double ht = fdtd_ht05;
    grid _grd = &fdtd_grd;
    double dey, dez;
    double mu_tmp = fun_mu(i, j, k, HX_ID);
    double ky, ay, by, phu_hx_y;
    double kz, az, bz, phu_hx_z;
    int ymin, ymax, j_pml;
    int zmin, zmax, k_pml;

    dey = (fdtd_Ey[j-1][k][i] - fdtd_Ey[j-1][k-1][i])/_grd->dzi05[k];
    dez = (fdtd_Ez[k-1][i][j] - fdtd_Ez[k-1][i][j-1])/_grd->dyi05[j];

    if (fdtd_comp.l_gran == BND_CPML) {
        ymin = (j <= PML_WIDTH);
        ymax = (j > _grd->ny - PML_WIDTH);

        zmin = (k <= PML_WIDTH);
        zmax = (k > _grd->nz - PML_WIDTH);

        if (!ymin && !ymax &&
            !zmin && !zmax) {
            fdtd_Hx[i][j-1][k-1] = fdtd_Hx[i][j-1][k-1] - ht*_cl*(dez - dey)/mu_tmp;
        } else {
            if (ymin) {
                ky = k_y_min[j-1];
                ay = alpha_i_h_min(j-1, Y_DIR);
                by = betta_i_h_min(j-1, Y_DIR);
                phu_hx_y = phu_hx_y_min[i][k-1][j-1] = by*phu_hx_y_min[i][k-1][j-1] + ay*dez;
            } else if (ymax) {
                j_pml = j - _grd->ny + PML_WIDTH;
                ky = k_y_max[j_pml-1];
                ay = alpha_i_h_max(j_pml-1, Y_DIR);
                by = betta_i_h_max(j_pml-1, Y_DIR);
                phu_hx_y = phu_hx_y_max[i][k-1][j_pml-1] = by*phu_hx_y_max[i][k-1][j_pml-1] + ay*dez;
            } else {
                ky = 1.0;
                ay = 0.0;
                by = 1.0;
                phu_hx_y = 0.0;
            }

            if (zmin) {
                kz = k_z_min[k-1];
                az = alpha_i_h_min(k-1, Z_DIR);
                bz = betta_i_h_min(k-1, Z_DIR);
                phu_hx_z = phu_hx_z_min[i][j-1][k-1] = bz*phu_hx_z_min[i][j-1][k-1] + az*dey;
            } else if (zmax) {
                k_pml = k - _grd->nz + PML_WIDTH;
                kz = k_z_max[k_pml-1];
                az = alpha_i_h_max(k_pml-1, Z_DIR);
                bz = betta_i_h_max(k_pml-1, Z_DIR);
                phu_hx_z = phu_hx_z_max[i][j-1][k_pml-1] = bz*phu_hx_z_max[i][j-1][k_pml-1] + az*dey;
            } else {
                kz = 1.0;
                az = 0.0;
                bz = 1.0;
                phu_hx_z = 0.0;
            }

            fdtd_Hx[i][j-1][k-1] = fdtd_Hx[i][j-1][k-1] - ht*_cl*((dez/ky - dey/kz) + (phu_hx_y - phu_hx_z))/mu_tmp;
        }
    } else {
        fdtd_Hx[i][j-1][k-1] = fdtd_Hx[i][j-1][k-1] - ht*_cl*(dez - dey)/mu_tmp;
    }
}

//вычисление Hy в точке
static void calc_Hy_point(int i, int j, int k)
{
    double ht = fdtd_ht05;
    grid _grd = &fdtd_grd;
    double dex, dez;
    double mu_tmp = fun_mu(i, j, k, HY_ID);
    double kx, ax, bx, phu_hy_x;
    double kz, az, bz, phu_hy_z;
    int xmin, xmax, i_pml;
    int zmin, zmax, k_pml;

    dex = (fdtd_Ex[i-1][j][k] - fdtd_Ex[i-1][j][k-1])/_grd->dzi05[k];
    dez = (fdtd_Ez[k-1][i][j] - fdtd_Ez[k-1][i-1][j])/_grd->dxi05[i];

    if (fdtd_comp.l_gran == BND_CPML) {
        xmin = (i <= PML_WIDTH);
        xmax = (i > _grd->nx - PML_WIDTH);

        zmin = (k <= PML_WIDTH);
        zmax = (k > _grd->nz - PML_WIDTH);

        if (!xmin && !xmax &&
            !zmin && !zmax) {
            fdtd_Hy[j][k-1][i-1] = fdtd_Hy[j][k-1][i-1] - ht*_cl*(dex - dez)/mu_tmp;
        } else {
            if (xmin) {
                kx = k_x_min[i-1];
                ax = alpha_i_h_min(i-1, X_DIR);
                bx = betta_i_h_min(i-1, X_DIR);
                phu_hy_x = phu_hy_x_min[j][k-1][i-1] = bx*phu_hy_x_min[j][k-1][i-1] + ax*dez;
            } else if(xmax) {
                i_pml = i - _grd->nx + PML_WIDTH;
                kx = k_x_max[i_pml-1];
                ax = alpha_i_h_max(i_pml-1, X_DIR);
                bx = betta_i_h_max(i_pml-1, X_DIR);
                phu_hy_x = phu_hy_x_max[j][k-1][i_pml-1] = bx*phu_hy_x_max[j][k-1][i_pml-1] + ax*dez;
            } else {
                kx = 1.0;
                ax = 0.0;
                bx = 1.0;
                phu_hy_x = 0.0;
            }

            if (zmin) {
                kz = k_z_min[k-1];
                az = alpha_i_h_min(k-1, Z_DIR);
                bz = betta_i_h_min(k-1, Z_DIR);
                phu_hy_z = phu_hy_z_min[j][i-1][k-1] = bz*phu_hy_z_min[j][i-1][k-1] + az*dex;
            } else if (zmax) {
                k_pml = k - _grd->nz + PML_WIDTH;
                kz = k_z_max[k_pml-1];
                az = alpha_i_h_max(k_pml-1, Z_DIR);
                bz = betta_i_h_max(k_pml-1, Z_DIR);
                phu_hy_z = phu_hy_z_max[j][i-1][k_pml-1] = bz*phu_hy_z_max[j][i-1][k_pml-1] + az*dex;
            } else {
                kz = 1.0;
                az = 0.0;
                bz = 1.0;
                phu_hy_z = 0.0;
            }

            fdtd_Hy[j][k-1][i-1] = fdtd_Hy[j][k-1][i-1] - ht*_cl*((dex/kz - dez/kx) + (phu_hy_z - phu_hy_x))/mu_tmp;
        }
    } else {
        fdtd_Hy[j][k-1][i-1] = fdtd_Hy[j][k-1][i-1] - ht*_cl*(dex - dez)/mu_tmp;
    }
}

//вычисление Hz в точке
static void calc_Hz_point(int i, int j, int k)
{
    double ht = fdtd_ht05;
    grid _grd = &fdtd_grd;
    double dex, dey;
    double mu_tmp = fun_mu(i, j, k, HZ_ID);
    double kx, ax, bx, phu_hz_x;
    double ky, ay, by, phu_hz_y;
    int xmin, xmax, i_pml;
    int ymin, ymax, j_pml;

    dex = (fdtd_Ex[i-1][j][k] - fdtd_Ex[i-1][j-1][k])/_grd->dyi05[j];
    dey = (fdtd_Ey[j-1][k][i] - fdtd_Ey[j-1][k][i-1])/_grd->dxi05[i];

    if (fdtd_comp.l_gran == BND_CPML) {
        xmin = (i <= PML_WIDTH);
        xmax = (i > _grd->nx - PML_WIDTH);

        ymin = (j <= PML_WIDTH);
        ymax = (j > _grd->ny - PML_WIDTH);

        if (!xmin && !xmax &&
            !ymin && !ymax) {
            fdtd_Hz[k][i-1][j-1] = fdtd_Hz[k][i-1][j-1] - ht*_cl*(dey - dex)/mu_tmp;
        } else {
            if (xmin) {
                kx = k_x_min[i-1];
                ax = alpha_i_h_min(i-1, X_DIR);
                bx = betta_i_h_min(i-1, X_DIR);
                phu_hz_x = phu_hz_x_min[k][j-1][i-1] = bx*phu_hz_x_min[k][j-1][i-1] + ax*dey;
            } else if (xmax) {
                i_pml = i - _grd->nx + PML_WIDTH;
                kx = k_x_max[i_pml-1];
                ax = alpha_i_h_max(i_pml-1, X_DIR);
                bx = betta_i_h_max(i_pml-1, X_DIR);
                phu_hz_x = phu_hz_x_max[k][j-1][i_pml-1] = bx*phu_hz_x_max[k][j-1][i_pml-1] + ax*dey;
            } else {
                kx = 1.0;
                ax = 0.0;
                bx = 1.0;
                phu_hz_x = 0.0;
            }

            if (ymin) {
                ky = k_y_min[j-1];
                ay = alpha_i_h_min(j-1, Y_DIR);
                by = betta_i_h_min(j-1, Y_DIR);
                phu_hz_y = phu_hz_y_min[k][i-1][j-1] = by*phu_hz_y_min[k][i-1][j-1] + ay*dex;
            }else if (ymax) {
                j_pml = j - _grd->ny + PML_WIDTH;
                ky = k_y_max[j_pml-1];
                ay = alpha_i_h_max(j_pml-1, Y_DIR);
                by = betta_i_h_max(j_pml-1, Y_DIR);
                phu_hz_y = phu_hz_y_max[k][i-1][j_pml-1] = by*phu_hz_y_max[k][i-1][j_pml-1] + ay*dex;
            } else {
                ky = 1.0;
                ay = 0.0;
                by = 1.0;
                phu_hz_y = 0.0;
            }

            fdtd_Hz[k][i-1][j-1] = fdtd_Hz[k][i-1][j-1] - ht*_cl*((dey/kx - dex/ky) + (phu_hz_x - phu_hz_y))/mu_tmp;
        }
    } else {
        fdtd_Hz[k][i-1][j-1] = fdtd_Hz[k][i-1][j-1] - ht*_cl*(dey - dex)/mu_tmp;
    }
}

void solverFDTD(int n, double time)
{
    int i, j, k;
    int jk, ki, ij;
    grid _grd = &fdtd_grd;
    comp _comp = &fdtd_comp;

    if (_comp->l_disp) {
        update_debye_data(fdtd_ht);
    }

    //вычисление Ex внутри расчетной области
    for (i = 1; i <= _grd->nx; ++i) {
#pragma omp parallel for private(j, k, jk)
        for (jk = 0; jk < (_grd->ny-1)*(_grd->nz-1); ++jk) {
            j = jk/(_grd->nz-1);
            k = jk%(_grd->nz-1);

            calc_Ex_point(i, j+1, k+1);
        }
    }

    if (_comp->l_tfsf) {
        //вычисление волновой границы Ex внутри расчетной области
        for (i = TFSF_WIDTH + 1; i <= _grd->nx - TFSF_WIDTH; ++i) {
            apply_Ex_tfsf(i);
        }
    }

    //вычисление Ex на границе Zmin, Zmax
    for (i = 1; i <= _grd->nx; ++i) {
        for (j = 0; j < (_grd->ny+1); ++j) {
            switch (_comp->zmin_bnd) {
            case BND_ZERO_E ://идеальный электрический проводник
            case BND_CPML ://CPML
                //Zmin
                fdtd_Ex[i-1][j][0] = 0.0;
                break;
            }

            switch (_comp->zmax_bnd) {
            case BND_ZERO_E ://идеальный электрический проводник
            case BND_CPML ://CPML
                //Zmax
                fdtd_Ex[i-1][j][_grd->nz] = 0.0;
                break;
            }
        }
    }

    //вычисление Ex на границе Ymin, Ymax
    for (i = 1; i <= _grd->nx; ++i) {
        for (k = 0; k < (_grd->nz+1); ++k) {
            switch (_comp->ymin_bnd) {
            case BND_ZERO_E ://идеальный электрический проводник
            case BND_CPML ://CPML
                //Ymin
                fdtd_Ex[i-1][0][k] = 0.0;
                break;
            }

            switch (_comp->ymax_bnd) {
            case BND_ZERO_E ://идеальный электрический проводник
            case BND_CPML ://CPML
                //Ymax
                fdtd_Ex[i-1][_grd->ny][k] = 0.0;
                break;
            }
        }
    }

    //вычисление Ey внутри расчетной области
    for (j = 1; j <= _grd->ny; ++j) {
#pragma omp parallel for private(k, i, ki)
        for (ki = 0; ki < (_grd->nz-1)*(_grd->nx-1); ++ki) {
            k = ki/(_grd->nx-1);
            i = ki%(_grd->nx-1);

            calc_Ey_point(i+1, j, k+1);
        }
    }

    if (_comp->l_tfsf) {
        //вычисление волновой границы Ey внутри расчетной области
        for (j = TFSF_WIDTH + 1; j <= _grd->ny - TFSF_WIDTH; ++j) {
            apply_Ey_tfsf(j);
        }
    }

    //вычисление Ey на границе Xmin, Xmax
    for (j = 1; j <= _grd->ny; ++j) {
        for (k = 0; k < (_grd->nz+1); ++k) {
            switch (_comp->xmin_bnd) {
            case BND_ZERO_E ://идеальный электрический проводник
            case BND_CPML ://CPML
                //Xmin
                fdtd_Ey[j-1][k][0] = 0.0;
                break;
            }

            switch (_comp->xmax_bnd) {
            case BND_ZERO_E ://идеальный электрический проводник
            case BND_CPML ://CPML
                //Xmax
                fdtd_Ey[j-1][k][_grd->nx] = 0.0;
                break;
            }
        }
    }

    //вычисление Ey на границе Zmin, Zmax
    for (j = 1; j <= _grd->ny; ++j) {
        for (i = 0; i < (_grd->nx+1); ++i) {
            switch (_comp->zmin_bnd) {
            case BND_ZERO_E ://идеальный электрический проводник
            case BND_CPML ://CPML
                //Zmin
                fdtd_Ey[j-1][0][i] = 0.0;
                break;
            }

            switch (_comp->zmax_bnd) {
            case BND_ZERO_E ://идеальный электрический проводник
            case BND_CPML ://CPML
                //Zmax
                fdtd_Ey[j-1][_grd->nz][i] = 0.0;
                break;
            }
        }
    }

    //вычисление Ez внутри расчетной области
    for (k = 1; k <= _grd->nz; ++k) {
#pragma omp parallel for private(i, j, ij)
        for (ij = 0; ij < (_grd->nx-1)*(_grd->ny-1); ++ij) {
            i = ij/(_grd->ny-1);
            j = ij%(_grd->ny-1);

            calc_Ez_point(i+1, j+1, k);
        }
    }

    if (_comp->l_tfsf) {
        //вычисление волновой границы Ez внутри расчетной области
        for (k = TFSF_WIDTH + 1; k <= _grd->nz - TFSF_WIDTH; ++k) {
            apply_Ez_tfsf(k);
        }
    }

    //вычисление Ez на границе Ymin, Ymax
    for (k = 1; k <= _grd->nz; ++k) {
        for (i = 0; i < (_grd->nx+1); ++i) {
            switch (_comp->ymin_bnd) {
            case BND_ZERO_E ://идеальный электрический проводник
            case BND_CPML ://CPML
                //Ymin
                fdtd_Ez[k-1][i][0] = 0.0;
                break;
            }

            switch (_comp->ymax_bnd) {
            case BND_ZERO_E ://идеальный электрический проводник
            case BND_CPML ://CPML
                //Ymax
                fdtd_Ez[k-1][i][_grd->ny] = 0.0;
                break;
            }
        }
    }

    //вычисление Ez на границе Xmin, Xmax
    for (k = 1; k <= _grd->nz; ++k) {
        for (j = 0; j < (_grd->ny+1); ++j) {
            switch (_comp->xmin_bnd) {
            case BND_ZERO_E ://идеальный электрический проводник
            case BND_CPML ://CPML
                //Xmin
                fdtd_Ez[k-1][0][j] = 0.0;
                break;
            }

            switch (_comp->xmax_bnd) {
            case BND_ZERO_E ://идеальный электрический проводник
            case BND_CPML ://CPML
                //Xmax
                fdtd_Ez[k-1][_grd->nx][j] = 0.0;
                break;
            }
        }
    }

    if (_comp->l_tfsf) {
        calc_plane_wave(time);
    }

    if (_comp->l_nfff) {
        load_nfff_domain();
    }

    for (i = 0; i <= _grd->nx; ++i) {
#pragma omp parallel for private(j, k, jk)
        for (jk = 0; jk < (_grd->ny)*(_grd->nz); ++jk) {
            j = jk/(_grd->nz);
            k = jk%(_grd->nz);
            calc_Hx_point(i, j+1, k+1);
        }
    }

    if (_comp->l_tfsf) {
        //вычисление волновой границы Hx внутри расчетной области
        for (i = TFSF_WIDTH; i <= _grd->nx - TFSF_WIDTH; ++i) {
            apply_Hx_tfsf(i);
        }
    }

    for (j = 0; j <= _grd->ny; ++j) {
#pragma omp parallel for private(k, i, ki)
        for (ki = 0; ki < (_grd->nz)*(_grd->nx); ++ki) {
            k = ki/(_grd->nx);
            i = ki%(_grd->nx);

            calc_Hy_point(i+1, j, k+1);
        }
    }

    if (_comp->l_tfsf) {
        //вычисление волновой границы Hy внутри расчетной области
        for (j = TFSF_WIDTH; j <= _grd->ny - TFSF_WIDTH; ++j) {
            apply_Hy_tfsf(j);
        }
    }

    for (k = 0; k <= _grd->nz; ++k) {
#pragma omp parallel for private(i, j, ij)
        for (ij = 0; ij < (_grd->nx)*(_grd->ny); ++ij) {
            i = ij/(_grd->ny);
            j = ij%(_grd->ny);

            calc_Hz_point(i+1, j+1, k);
        }
    }

    if (_comp->l_tfsf) {
        //вычисление волновой границы Hz внутри расчетной области
        for (k = TFSF_WIDTH; k <= _grd->nz - TFSF_WIDTH; ++k) {
            apply_Hz_tfsf(k);
        }
    }

    if (_comp->l_nfff) {
        write_nfff_domain(n, time);
    }
}
