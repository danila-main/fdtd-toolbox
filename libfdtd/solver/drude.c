
/*! \file 
* \brief Drude dispersion model.
* \author А.Kryukov anton.krv@gmail.com
*/

#include "globals.h"
#include "constants.h"
#include "functions.h"
#include "fdtd.h"

/*!
* \brief Структура для описания параметров материала.
*/
typedef struct
{
    int nPoles;//число полюсов
    double *res_freq;//резонансная частота на данном полюсе
    double *col_freq;//обратное время релаксации на данном полюсе
    // время релаксации не может быть нулевым, а вот частота столкновений в плазме может

    double *am;//вспомогательные массивы
    double *bm;
    double sig;//эффективная проводимость
} _drudeMaterial;

static _drudeMaterial matData;//материал

/*!
* \brief Структура для хранения параметров расчета в точке.
*/
typedef struct
{
    _drudeMaterial *material;//материал
    double *Jpol;//компоненты поляризационного тока для каждого полюса
} _drude;

_drude ****drude_Jx;/*!< массив x-компонеты поляризационного тока для электрического поля. */
_drude ****drude_Jy;/*!< массив y-компонеты поляризационного тока для электрического поля. */
_drude ****drude_Jz;/*!< массив z-компонеты поляризационного тока для электрического поля. */

static _drude *_init_data(_drudeMaterial *material)
{
    _drude *data;
    PARMS_malloc(data, 1, _drude);

    data->material = material;
    PARMS_malloc(data->Jpol, material->nPoles, double);

    memset(data->Jpol, 0, material->nPoles*sizeof(double));

    return (data);
}

//инициализация дисперсионного материала
void init_drude_data()
{
    int i, j, k;
    grid _grd = &fdtd_grd;

    matData.nPoles = 1;
    PARMS_malloc(matData.res_freq, matData.nPoles, double);
    PARMS_malloc(matData.col_freq, matData.nPoles, double);
    PARMS_malloc(matData.am, matData.nPoles, double);
    PARMS_malloc(matData.bm, matData.nPoles, double);

    for (i = 0; i < matData.nPoles; ++i) {
        matData.res_freq[i] = 1.0e+09;
        matData.col_freq[i] = 1.0e+02;
        matData.am[i] = 0.0;
        matData.bm[i] = 0.0;
    }

    matData.sig = 0.0;

    PARMS_malloc(drude_Jx, _grd->nx, _drude***);
    for (i = 0; i < _grd->nx; ++i) {
        PARMS_malloc(drude_Jx[i], _grd->ny+1, _drude**);
        for (j = 0; j < _grd->ny+1; ++j) {
            PARMS_malloc(drude_Jx[i][j], _grd->nz+1, _drude*);
            for (k = 0; k < _grd->nz+1; ++k) {
                drude_Jx[i][j][k] = _init_data(&matData);
            }
        }
    }

    PARMS_malloc(drude_Jy, _grd->ny, _drude***);
    for (j = 0; j < _grd->ny; ++j) {
        PARMS_malloc(drude_Jy[j], _grd->nz+1, _drude**);
        for (k = 0; k < _grd->nz+1; ++k) {
            PARMS_malloc(drude_Jy[j][k], _grd->nx+1, _drude*);
            for (i = 0; i < _grd->nx+1; ++i) {
                drude_Jy[j][k][i] = _init_data(&matData);
            }
        }
    }

    PARMS_malloc(drude_Jz, _grd->nz, _drude***);
    for (k = 0; k < _grd->nz; ++k) {
        PARMS_malloc(drude_Jz[k], _grd->nx+1, _drude**);
        for (i = 0; i < _grd->nx+1; ++i) {
            PARMS_malloc(drude_Jz[k][i], _grd->ny+1, _drude*);
            for (j = 0; j < _grd->ny+1; ++j) {
                drude_Jz[k][i][j] = _init_data(&matData);
            }
        }
    }
}

//расчет дисперсионного материала на каждом шаге
void update_drude_data(double ht)
{
    int i;

    matData.sig = 0.0;

    for (i = 0; i < matData.nPoles; ++i) {
        matData.am[i] = (2.0 - ht*matData.col_freq[i])/(2.0 + ht*matData.col_freq[i]);
        matData.bm[i] = _pi4*((matData.res_freq[i]*matData.res_freq[i]*ht)/(2.0 + ht*matData.col_freq[i]));

        matData.sig += matData.bm[i];
    }

    matData.sig *= 0.5;//получили поправку для проводимости
}

//никогда не вызывать эти функции в CPML слое!!!!!

//вычисление Ex в точке
void calc_Ex_point_drude(int i, int j, int k)
{
    double ht = fdtd_ht;
    grid _grd = &fdtd_grd;
    double dhy, dhz;
    double sig_tmp = _pi4*fun_sig(i, j, k, EX_ID);
    double eps_tmp = fun_eps(i, j, k, EX_ID);// считаем, что диэлетрическая проницаемость для бесконечной частоты задана базовым материалом
    int l;
    double polarJx = 0.0;
    double Ex;
    _drudeMaterial *_mat = drude_Jx[i-1][j][k]->material;

    //вычисляем поляризационный ток
    for (l = 0; l < _mat->nPoles; ++l) {
        polarJx += drude_Jx[i-1][j][k]->Jpol[l]*(1.0 + _mat->am[l]);
    }

    polarJx *= 0.5;

    sig_tmp += _mat->sig;//добавили дисперсионную поправку

    dhy = (fdtd_Hy[j][k][i-1] - fdtd_Hy[j][k-1][i-1])/_grd->dzi[k];
    dhz = (fdtd_Hz[k][i-1][j] - fdtd_Hz[k][i-1][j-1])/_grd->dyi[j];

    Ex = fdtd_Ex[i-1][j][k];//Ex[n]

    fdtd_Ex[i-1][j][k] = ((eps_tmp - 0.5*ht*sig_tmp)*fdtd_Ex[i-1][j][k]
                    + ht*_cl*(dhz - dhy) - ht*(_pi4*fdtd_Jx[i-1][j][k] + polarJx))/(eps_tmp + 0.5*ht*sig_tmp);

    Ex += fdtd_Ex[i-1][j][k];//Ex[n+1] + Ex[n]
    Ex *= 0.5;//(Ex[n+1] + Ex[n])/2

    for (l = 0; l < _mat->nPoles; ++l) {//расчет поляризационного тока на след шаге
        drude_Jx[i-1][j][k]->Jpol[l] = _mat->am[l]*drude_Jx[i-1][j][k]->Jpol[l] + _mat->bm[l]*Ex;
    }
}

//вычисление Ey в точке
void calc_Ey_point_drude(int i, int j, int k)
{
    double ht = fdtd_ht;
    grid _grd = &fdtd_grd;
    double dhx, dhz;
    double sig_tmp = _pi4*fun_sig(i, j, k, EY_ID);
    double eps_tmp = fun_eps(i, j, k, EY_ID);// считаем, что диэлетрическая проницаемость для бесконечной частоты задана базовым материалом
    int l;
    double polarJy = 0.0;
    double Ey;
    _drudeMaterial *_mat = drude_Jy[j-1][k][i]->material;

    //вычисляем поляризационный ток
    for (l = 0; l < _mat->nPoles; ++l) {
        polarJy += drude_Jy[j-1][k][i]->Jpol[l]*(1.0 + _mat->am[l]);
    }

    polarJy *= 0.5;

    sig_tmp += _mat->sig;//добавили дисперсионную поправку

    dhx = (fdtd_Hx[i][j-1][k] - fdtd_Hx[i][j-1][k-1])/_grd->dzi[k];
    dhz = (fdtd_Hz[k][i][j-1] - fdtd_Hz[k][i-1][j-1])/_grd->dxi[i];

    Ey = fdtd_Ey[j-1][k][i];//Ey[n]

    fdtd_Ey[j-1][k][i] = ((eps_tmp - 0.5*ht*sig_tmp)*fdtd_Ey[j-1][k][i]
                    + ht*_cl*(dhx - dhz) - ht*(_pi4*fdtd_Jy[j-1][k][i] + polarJy))/(eps_tmp + 0.5*ht*sig_tmp);

    Ey += fdtd_Ey[j-1][k][i];//Ey[n+1] + Ey[n]
    Ey *= 0.5;//(Ey[n+1] + Ey[n])/2

    for (l = 0; l < _mat->nPoles; ++l) {//расчет поляризационного тока на след шаге
        drude_Jy[j-1][k][i]->Jpol[l] = _mat->am[l]*drude_Jy[j-1][k][i]->Jpol[l] + _mat->bm[l]*Ey;
    }
}

//вычисление Ez в точке
void calc_Ez_point_drude(int i, int j, int k)
{
    double ht = fdtd_ht;
    grid _grd = &fdtd_grd;
    double dhx, dhy;
    double sig_tmp = _pi4*fun_sig(i, j, k, EZ_ID);
    double eps_tmp = fun_eps(i, j, k, EZ_ID);// считаем, что диэлетрическая проницаемость для бесконечной частоты задана базовым материалом
    int l;
    double polarJz = 0.0;
    double Ez;
    _drudeMaterial *_mat = drude_Jz[k-1][i][j]->material;

    //вычисляем поляризационный ток
    for (l = 0; l < _mat->nPoles; ++l) {
        polarJz += drude_Jz[k-1][i][j]->Jpol[l]*(1.0 + _mat->am[l]);
    }

    polarJz *= 0.5;

    sig_tmp += _mat->sig;//добавили дисперсионную поправку

    dhx = (fdtd_Hx[i][j][k-1] - fdtd_Hx[i][j-1][k-1])/_grd->dyi[j];
    dhy = (fdtd_Hy[j][k-1][i] - fdtd_Hy[j][k-1][i-1])/_grd->dxi[i];

    Ez = fdtd_Ez[k-1][i][j];//Ez[n]

    fdtd_Ez[k-1][i][j] = ((eps_tmp - 0.5*ht*sig_tmp)*fdtd_Ez[k-1][i][j]
                    + ht*_cl*(dhy - dhx) - ht*(_pi4*fdtd_Jz[k-1][i][j] + polarJz))/(eps_tmp + 0.5*ht*sig_tmp);

    Ez += fdtd_Ez[k-1][i][j];//Ez[n+1] + Ez[n]
    Ez *= 0.5;//(Ez[n+1] + Ez[n])/2

    for (l = 0; l < _mat->nPoles; ++l) {
        drude_Jz[k-1][i][j]->Jpol[l] = _mat->am[l]*drude_Jz[k-1][i][j]->Jpol[l] + _mat->bm[l]*Ez;
    }
}
