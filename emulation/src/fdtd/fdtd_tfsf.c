
/*! \file 
* \brief Функции расчета плоской падающей волны.
* \author Крюков А.А. anton.krv@gmail.com
*/

#include <math.h>
#include "../../include/globals.h"
#include "../../include/constants.h"
#include "../../include/functions.h"
#include "../../include/fdtd.h"

static double *E;
static double *H;

//вычисление E в точке для падающей волны
static double calc_E_wave(double time)
{
    return sin(time*_pi*5.0e08);
}

void init_plane_wave()
{
    grid _grd = &fdtd_grd;

    PARMS_malloc(E, _grd->nz+1, double);
    memset(E, 0, (_grd->nz+1)*sizeof(double));

    PARMS_malloc(H, _grd->nz, double);
    memset(H, 0, (_grd->nz)*sizeof(double));
}

void calc_plane_wave(double time)
{
    double ht = fdtd_ht;
    grid _grd = &fdtd_grd;
    double dh, de;
    int k;
    double E1 = 0.0;
    double Enz1 = 0.0;

    E1 = E[1];
    Enz1 = E[_grd->nz-1];

    for (k = 1; k < _grd->nz; ++k) {
        dh = (H[k] - H[k-1])/_grd->dzi[k];
        E[k] = E[k] - ht*_cl*dh;

        if (k == 2)
            E[k] = calc_E_wave(time);
    }

    E[0] = E1 + ((_cl*ht - _grd->dzi05[1])/(_cl*ht + _grd->dzi05[1]))*(E[1] - E[0]);
    E[_grd->nz] = Enz1 + ((_cl*ht - _grd->dzi05[_grd->nz])/(_cl*ht + _grd->dzi05[_grd->nz]))*(E[_grd->nz-1] - E[_grd->nz]);

    for (k = 1; k <= _grd->nz; ++k) {
        de = (E[k] - E[k-1])/_grd->dzi05[k];
        H[k-1] = H[k-1] - ht*_cl*de;
    }
}

static double Ex_wave(int i, int j, int k)
{
    double val = E[k];

    return val;
}

//static double Ey_wave(int i, int j, int k)
//{
//    double val = 0.0;
//
//    return val;
//}

//static double Ez_wave(int i, int j, int k)
//{
//    double val = 0.0;
//
//    return val;
//}

//static double Hx_wave(int i, int j, int k)
//{
//    double val = 0.0;
//
//    return val;
//}

static double Hy_wave(int i, int j, int k)
{
    double val = H[k];

    return val;
}

//static double Hz_wave(int i, int j, int k)
//{
//    double val = 0.0;
//
//    return val;
//}

//учет периметра волновой границы Ex
void apply_Ex_tfsf(int i)
{
    int j, k;
    double ht = fdtd_ht;
    grid _grd = &fdtd_grd;
    double sig_tmp;
    double eps_tmp;

    for (j = TFSF_WIDTH; j <= _grd->ny - TFSF_WIDTH; ++j) {
        k = TFSF_WIDTH;//zmin
        //вводим поправку для производной Hy по z, т.к. это граница волновой области
        sig_tmp = _pi4*fun_sig(i, j, k, EX_ID);
        eps_tmp = fun_eps(i, j, k, EX_ID);

        fdtd_Ex[i-1][j][k] += (ht*_cl*(Hy_wave(i, j, k)/_grd->dzi[k]))/(eps_tmp + 0.5*ht*sig_tmp);

        k = _grd->nz - TFSF_WIDTH;//zmax
        //вводим поправку для производной Hy по z, т.к. это граница волновой области
        sig_tmp = _pi4*fun_sig(i, j, k, EX_ID);
        eps_tmp = fun_eps(i, j, k, EX_ID);

        fdtd_Ex[i-1][j][k] -= (ht*_cl*(Hy_wave(i, j, k+1)/_grd->dzi[k]))/(eps_tmp + 0.5*ht*sig_tmp);
    }

    //for (k = TFSF_WIDTH; k <= _grd->nz - TFSF_WIDTH; ++k) {
    //    j = TFSF_WIDTH;//ymin
    //    //вводим поправку для производной Hz по y, т.к. это граница волновой области
    //    sig_tmp = _pi4*fun_sig(i, j, k, EX_ID);
    //    eps_tmp = fun_eps(i, j, k, EX_ID);

    //    fdtd_Ex[i-1][j][k] -= (ht*_cl*(Hz_wave(i, j, k)/_grd->dyi[j]))/(eps_tmp + 0.5*ht*sig_tmp);

    //    j = _grd->ny - TFSF_WIDTH;//ymax
    //    //вводим поправку для производной Hz по y, т.к. это граница волновой области
    //    sig_tmp = _pi4*fun_sig(i, j, k, EX_ID);
    //    eps_tmp = fun_eps(i, j, k, EX_ID);

    //    fdtd_Ex[i-1][j][k] += (ht*_cl*(Hz_wave(i, j+1, k)/_grd->dyi[j]))/(eps_tmp + 0.5*ht*sig_tmp);
    //}
}

//учет периметра волновой границы Ey
void apply_Ey_tfsf(int j)
{
    //int i, k;
    //double ht = fdtd_ht;
    //grid _grd = &fdtd_grd;
    //double sig_tmp;
    //double eps_tmp;

    //for (i = TFSF_WIDTH; i <= _grd->nx - TFSF_WIDTH; ++i) {
    //    k = TFSF_WIDTH;//zmin
    //    //вводим поправку для производной Hx по z, т.к. это граница волновой области
    //    sig_tmp = _pi4*fun_sig(i, j, k, EY_ID);
    //    eps_tmp = fun_eps(i, j, k, EY_ID);

    //    fdtd_Ey[j-1][k][i] -= (ht*_cl*(Hx_wave(i, j, k)/_grd->dzi[k]))/(eps_tmp + 0.5*ht*sig_tmp);

    //    k = _grd->nz - TFSF_WIDTH;//zmax
    //    //вводим поправку для производной Hx по z, т.к. это граница волновой области
    //    sig_tmp = _pi4*fun_sig(i, j, k, EY_ID);
    //    eps_tmp = fun_eps(i, j, k, EY_ID);

    //    fdtd_Ey[j-1][k][i] += (ht*_cl*(Hx_wave(i, j, k+1)/_grd->dzi[k]))/(eps_tmp + 0.5*ht*sig_tmp);
    //}

    //for (k = TFSF_WIDTH; k <= _grd->nz - TFSF_WIDTH; ++k) {
    //    i = TFSF_WIDTH;//xmin
    //    //вводим поправку для производной Hz по x, т.к. это граница волновой области
    //    sig_tmp = _pi4*fun_sig(i, j, k, EY_ID);
    //    eps_tmp = fun_eps(i, j, k, EY_ID);

    //    fdtd_Ey[j-1][k][i] += (ht*_cl*(Hz_wave(i, j, k)/_grd->dxi[i]))/(eps_tmp + 0.5*ht*sig_tmp);

    //    i = _grd->nx - TFSF_WIDTH;//xmax
    //    //вводим поправку для производной Hz по x, т.к. это граница волновой области
    //    sig_tmp = _pi4*fun_sig(i, j, k, EY_ID);
    //    eps_tmp = fun_eps(i, j, k, EY_ID);

    //    fdtd_Ey[j-1][k][i] -= (ht*_cl*(Hz_wave(i+1, j, k)/_grd->dxi[i]))/(eps_tmp + 0.5*ht*sig_tmp);
    //}
}

//учет периметра волновой границы Ez
void apply_Ez_tfsf(int k)
{
    int i, j;
    double ht = fdtd_ht;
    grid _grd = &fdtd_grd;
    double sig_tmp;
    double eps_tmp;

    //for (i = TFSF_WIDTH; i <= _grd->nx - TFSF_WIDTH; ++i) {
    //    j = TFSF_WIDTH;//ymin
    //    //вводим поправку для производной Hx по y, т.к. это граница волновой области
    //    sig_tmp = _pi4*fun_sig(i, j, k, EZ_ID);
    //    eps_tmp = fun_eps(i, j, k, EZ_ID);

    //    fdtd_Ez[k-1][i][j] += (ht*_cl*(Hx_wave(i, j, k)/_grd->dyi[j]))/(eps_tmp + 0.5*ht*sig_tmp);


    //    j = _grd->ny - TFSF_WIDTH;//ymax
    //    //вводим поправку для производной Hx по y, т.к. это граница волновой области
    //    sig_tmp = _pi4*fun_sig(i, j, k, EZ_ID);
    //    eps_tmp = fun_eps(i, j, k, EZ_ID);

    //    fdtd_Ez[k-1][i][j] -= (ht*_cl*(Hx_wave(i, j+1, k)/_grd->dyi[j]))/(eps_tmp + 0.5*ht*sig_tmp);
    //}

    for (j = TFSF_WIDTH; j <= _grd->ny - TFSF_WIDTH; ++j) {
        i = TFSF_WIDTH;//xmin
        //вводим поправку для производной Hy по x, т.к. это граница волновой области
        sig_tmp = _pi4*fun_sig(i, j, k, EZ_ID);
        eps_tmp = fun_eps(i, j, k, EZ_ID);

        fdtd_Ez[k-1][i][j] -= (ht*_cl*(Hy_wave(i, j, k)/_grd->dxi[i]))/(eps_tmp + 0.5*ht*sig_tmp);

        i = _grd->nx - TFSF_WIDTH;//xmax
        //вводим поправку для производной Hy по x, т.к. это граница волновой области
        sig_tmp = _pi4*fun_sig(i, j, k, EZ_ID);
        eps_tmp = fun_eps(i, j, k, EZ_ID);

        fdtd_Ez[k-1][i][j] += (ht*_cl*(Hy_wave(i+1, j, k)/_grd->dxi[i]))/(eps_tmp + 0.5*ht*sig_tmp);
    }
}

//учет периметра волновой границы Hx
void apply_Hx_tfsf(int i)
{
    //int j, k;
    //double ht = fdtd_ht05;
    //grid _grd = &fdtd_grd;
    //double mu_tmp;

    //for (j = TFSF_WIDTH + 1; j <= _grd->ny - TFSF_WIDTH; ++j) {
    //    k = TFSF_WIDTH;//zmin
    //    //вводим поправку для производной Ey по z, т.к. это граница волновой области
    //    mu_tmp = fun_mu(i, j, k, HX_ID);

    //    fdtd_Hx[i][j-1][k-1] -= ht*_cl*(Ey_wave(i, j, k)/_grd->dzi05[k])/mu_tmp;

    //    k = _grd->nz - TFSF_WIDTH + 1;//zmax
    //    //вводим поправку для производной Ey по z, т.к. это граница волновой области
    //    mu_tmp = fun_mu(i, j, k, HX_ID);

    //    fdtd_Hx[i][j-1][k-1] += ht*_cl*(Ey_wave(i, j, k-1)/_grd->dzi05[k])/mu_tmp;
    //}

    //for (k = TFSF_WIDTH + 1; k <= _grd->nz - TFSF_WIDTH; ++k) {
    //    j = TFSF_WIDTH;//ymin
    //    //вводим поправку для производной Ez по y, т.к. это граница волновой области
    //    mu_tmp = fun_mu(i, j, k, HX_ID);

    //    fdtd_Hx[i][j-1][k-1] += ht*_cl*(Ez_wave(i, j, k)/_grd->dyi05[j])/mu_tmp;

    //    j = _grd->ny - TFSF_WIDTH + 1;//ymax
    //    //вводим поправку для производной Ez по y, т.к. это граница волновой области
    //    mu_tmp = fun_mu(i, j, k, HX_ID);

    //    fdtd_Hx[i][j-1][k-1] -= ht*_cl*(Ez_wave(i, j-1, k)/_grd->dyi05[j])/mu_tmp;
    //}
}

//учет периметра волновой границы Hy
void apply_Hy_tfsf(int j)
{
    int i, k;
    double ht = fdtd_ht05;
    grid _grd = &fdtd_grd;
    double mu_tmp;

    for (i = TFSF_WIDTH + 1; i <= _grd->nx - TFSF_WIDTH; ++i) {
        k = TFSF_WIDTH;//zmin
        //вводим поправку для производной Ex по z, т.к. это граница волновой области
        mu_tmp = fun_mu(i, j, k, HY_ID);

        fdtd_Hy[j][k-1][i-1] += ht*_cl*(Ex_wave(i, j, k)/_grd->dzi05[k])/mu_tmp;

        k = _grd->nz - TFSF_WIDTH + 1;//zmax
        //вводим поправку для производной Ex по z, т.к. это граница волновой области
        mu_tmp = fun_mu(i, j, k, HY_ID);

        fdtd_Hy[j][k-1][i-1] -= ht*_cl*(Ex_wave(i, j, k-1)/_grd->dzi05[k])/mu_tmp;
    }

    //for (k = TFSF_WIDTH + 1; k <= _grd->nz - TFSF_WIDTH; ++k) {
    //    i = TFSF_WIDTH;//xmin
    //    //вводим поправку для производной Ez по x, т.к. это граница волновой области
    //    mu_tmp = fun_mu(i, j, k, HY_ID);

    //    fdtd_Hy[j][k-1][i-1] -= ht*_cl*(Ez_wave(i, j, k)/_grd->dxi05[i])/mu_tmp;

    //    i = _grd->nx - TFSF_WIDTH + 1;//xmax
    //    //вводим поправку для производной Ez по x, т.к. это граница волновой области
    //    mu_tmp = fun_mu(i, j, k, HY_ID);

    //    fdtd_Hy[j][k-1][i-1] += ht*_cl*(Ez_wave(i-1, j, k)/_grd->dxi05[i])/mu_tmp;
    //}
}

//учет периметра волновой границы Hz
void apply_Hz_tfsf(int k)
{
    int i, j;
    double ht = fdtd_ht05;
    grid _grd = &fdtd_grd;
    double mu_tmp;

    for (i = TFSF_WIDTH + 1; i <= _grd->nx - TFSF_WIDTH; ++i) {
        j = TFSF_WIDTH;//ymin
        //вводим поправку для производной Ex по y, т.к. это граница волновой области
        mu_tmp = fun_mu(i, j, k, HZ_ID);

        fdtd_Hz[k][i-1][j-1] -= ht*_cl*(Ex_wave(i, j, k)/_grd->dyi05[j])/mu_tmp;

        j = _grd->ny - TFSF_WIDTH + 1;//ymax
        //вводим поправку для производной Ex по y, т.к. это граница волновой области
        mu_tmp = fun_mu(i, j, k, HZ_ID);

        fdtd_Hz[k][i-1][j-1] += ht*_cl*(Ex_wave(i, j-1, k)/_grd->dyi05[j])/mu_tmp;
    }

    //for (j = TFSF_WIDTH + 1; j <= _grd->ny - TFSF_WIDTH; ++j) {
    //    i = TFSF_WIDTH;//xmin
    //    //вводим поправку для производной Ey по x, т.к. это граница волновой области
    //    mu_tmp = fun_mu(i, j, k, HZ_ID);

    //    fdtd_Hz[k][i-1][j-1] += ht*_cl*(Ey_wave(i, j, k)/_grd->dxi05[i])/mu_tmp;

    //    i = _grd->nx - TFSF_WIDTH + 1;//xmax
    //    //вводим поправку для производной Ey по x, т.к. это граница волновой области
    //    mu_tmp = fun_mu(i, j, k, HZ_ID);

    //    fdtd_Hz[k][i-1][j-1] -= ht*_cl*(Ey_wave(i-1, j, k)/_grd->dxi05[i])/mu_tmp;
    //}
}
