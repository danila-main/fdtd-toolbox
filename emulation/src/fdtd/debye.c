
/*! \file 
* \brief функции расчета дисперсионной модели Дебая.
* \author Крюков А.А. anton.krv@gmail.com
*/

#include "../../include/globals.h"
#include "../../include/constants.h"
#include "../../include/functions.h"
#include "../../include/fdtd.h"

/*!
* \brief Структура для описания параметров материала.
*/
typedef struct
{
    int nPoles;//число полюсов
    double *deltaEps;//изменение диэлектрической проницаемости на данном полюсе
    double *tauEps;//время релаксации на данном полюсе

    double *km;//вспомогательные массивы
    double *bm;
    double eps;//эффективная диэлектрическая проницаемость
} _debyeMaterial;

static _debyeMaterial matData;//материал

/*!
* \brief Структура для хранения параметров расчета в точке.
*/
typedef struct
{
    _debyeMaterial *material;//материал
    double *Jpol;//компоненты поляризационного тока для каждого полюса
} _debye;

_debye ****debye_Jx;/*!< массив x-компонеты поляризационного тока для электрического поля. */
_debye ****debye_Jy;/*!< массив y-компонеты поляризационного тока для электрического поля. */
_debye ****debye_Jz;/*!< массив z-компонеты поляризационного тока для электрического поля. */

static _debye *_init_data(_debyeMaterial *material)
{
    _debye *data;
    PARMS_malloc(data, 1, _debye);

    data->material = material;
    PARMS_malloc(data->Jpol, material->nPoles, double);

    memset(data->Jpol, 0, material->nPoles*sizeof(double));

    return (data);
}

//инициализация дисперсионного материала
void init_debye_data()
{
    int i, j, k;
    grid _grd = &fdtd_grd;

    matData.nPoles = 1;
    PARMS_malloc(matData.deltaEps, matData.nPoles, double);
    PARMS_malloc(matData.tauEps, matData.nPoles, double);
    PARMS_malloc(matData.km, matData.nPoles, double);
    PARMS_malloc(matData.bm, matData.nPoles, double);

    for (i = 0; i < matData.nPoles; ++i) {
        matData.deltaEps[i] = 1.0;
        matData.tauEps[i] = 1.0/(_pi*5.0e08);
        matData.km[i] = 0.0;
        matData.bm[i] = 0.0;
    }

    matData.eps = 0.0;

    PARMS_malloc(debye_Jx, _grd->nx, _debye***);
    for (i = 0; i < _grd->nx; ++i) {
        PARMS_malloc(debye_Jx[i], _grd->ny+1, _debye**);
        for (j = 0; j < _grd->ny+1; ++j) {
            PARMS_malloc(debye_Jx[i][j], _grd->nz+1, _debye*);
            for (k = 0; k < _grd->nz+1; ++k) {
                debye_Jx[i][j][k] = _init_data(&matData);
            }
        }
    }

    PARMS_malloc(debye_Jy, _grd->ny, _debye***);
    for (j = 0; j < _grd->ny; ++j) {
        PARMS_malloc(debye_Jy[j], _grd->nz+1, _debye**);
        for (k = 0; k < _grd->nz+1; ++k) {
            PARMS_malloc(debye_Jy[j][k], _grd->nx+1, _debye*);
            for (i = 0; i < _grd->nx+1; ++i) {
                debye_Jy[j][k][i] = _init_data(&matData);
            }
        }
    }

    PARMS_malloc(debye_Jz, _grd->nz, _debye***);
    for (k = 0; k < _grd->nz; ++k) {
        PARMS_malloc(debye_Jz[k], _grd->nx+1, _debye**);
        for (i = 0; i < _grd->nx+1; ++i) {
            PARMS_malloc(debye_Jz[k][i], _grd->ny+1, _debye*);
            for (j = 0; j < _grd->ny+1; ++j) {
                debye_Jz[k][i][j] = _init_data(&matData);
            }
        }
    }
}

//расчет дисперсионного материала на каждом шаге
void update_debye_data(double ht)
{
    int i;

    matData.eps = 0.0;

    for (i = 0; i < matData.nPoles; ++i) {
        matData.km[i] = (2.0*matData.tauEps[i] - ht)/(2.0*matData.tauEps[i] + ht);
        matData.bm[i] = _pi4*((2.0*matData.deltaEps[i]*ht)/(2.0*matData.tauEps[i] + ht));

        matData.eps += matData.bm[i];
    }

    matData.eps *= 0.5;//получили поправку для проницаемости
}

//никогда не вызывать эти функции в CPML слое!!!!!

//вычисление Ex в точке
void calc_Ex_point_debye(int i, int j, int k)
{
    double ht = fdtd_ht;
    grid _grd = &fdtd_grd;
    double dhy, dhz;
    double sig_tmp = _pi4*fun_sig(i, j, k, EX_ID);
    double eps_tmp = fun_eps(i, j, k, EX_ID);// считаем, что диэлетрическая проницаемость для бесконечной частоты задана базовым материалом
    int l;
    double polarJx = 0.0;
    double dEx;
    _debyeMaterial *_mat = debye_Jx[i-1][j][k]->material;

    //вычисляем поляризационный ток
    for (l = 0; l < _mat->nPoles; ++l) {
        polarJx += debye_Jx[i-1][j][k]->Jpol[l]*(1.0 + _mat->km[l]);
    }

    polarJx += 0.5;

    eps_tmp += _mat->eps;//добавили дисперсионную поправку

    dhy = (fdtd_Hy[j][k][i-1] - fdtd_Hy[j][k-1][i-1])/_grd->dzi[k];
    dhz = (fdtd_Hz[k][i-1][j] - fdtd_Hz[k][i-1][j-1])/_grd->dyi[j];

    dEx = -fdtd_Ex[i-1][j][k];//-Ex[n]

    fdtd_Ex[i-1][j][k] = ((eps_tmp - 0.5*ht*sig_tmp)*fdtd_Ex[i-1][j][k]
                    + ht*_cl*(dhz - dhy) - ht*(_pi4*fdtd_Jx[i-1][j][k] + polarJx))/(eps_tmp + 0.5*ht*sig_tmp);

    dEx += fdtd_Ex[i-1][j][k];//Ex[n+1] - Ex[n]
    dEx /= ht;//(Ex[n+1] - Ex[n])/ht

    for (l = 0; l < _mat->nPoles; ++l) {//расчет поляризационного тока на след шаге
        debye_Jx[i-1][j][k]->Jpol[l] = _mat->km[l]*debye_Jx[i-1][j][k]->Jpol[l] + _mat->bm[l]*dEx;
    }
}

//вычисление Ey в точке
void calc_Ey_point_debye(int i, int j, int k)
{
    double ht = fdtd_ht;
    grid _grd = &fdtd_grd;
    double dhx, dhz;
    double sig_tmp = _pi4*fun_sig(i, j, k, EY_ID);
    double eps_tmp = fun_eps(i, j, k, EY_ID);// считаем, что диэлетрическая проницаемость для бесконечной частоты задана базовым материалом
    int l;
    double polarJy = 0.0;
    double dEy;
    _debyeMaterial *_mat = debye_Jy[j-1][k][i]->material;

    //вычисляем поляризационный ток
    for (l = 0; l < _mat->nPoles; ++l) {
        polarJy += debye_Jy[j-1][k][i]->Jpol[l]*(1.0 + _mat->km[l]);
    }

    polarJy += 0.5;

    eps_tmp += _mat->eps;//добавили дисперсионную поправку

    dhx = (fdtd_Hx[i][j-1][k] - fdtd_Hx[i][j-1][k-1])/_grd->dzi[k];
    dhz = (fdtd_Hz[k][i][j-1] - fdtd_Hz[k][i-1][j-1])/_grd->dxi[i];

    dEy = -fdtd_Ey[j-1][k][i];//-Ey[n]

    fdtd_Ey[j-1][k][i] = ((eps_tmp - 0.5*ht*sig_tmp)*fdtd_Ey[j-1][k][i]
                    + ht*_cl*(dhx - dhz) - ht*(_pi4*fdtd_Jy[j-1][k][i] + polarJy))/(eps_tmp + 0.5*ht*sig_tmp);

    dEy += fdtd_Ey[j-1][k][i];//Ey[n+1] - Ey[n]
    dEy /= ht;//(Ey[n+1] - Ey[n])/ht

    for (l = 0; l < _mat->nPoles; ++l) {//расчет поляризационного тока на след шаге
        debye_Jy[j-1][k][i]->Jpol[l] = _mat->km[l]*debye_Jy[j-1][k][i]->Jpol[l] + _mat->bm[l]*dEy;
    }
}

//вычисление Ez в точке
void calc_Ez_point_debye(int i, int j, int k)
{
    double ht = fdtd_ht;
    grid _grd = &fdtd_grd;
    double dhx, dhy;
    double sig_tmp = _pi4*fun_sig(i, j, k, EZ_ID);
    double eps_tmp = fun_eps(i, j, k, EZ_ID);// считаем, что диэлетрическая проницаемость для бесконечной частоты задана базовым материалом
    int l;
    double polarJz = 0.0;
    double dEz;
    _debyeMaterial *_mat = debye_Jz[k-1][i][j]->material;

    //вычисляем поляризационный ток
    for (l = 0; l < _mat->nPoles; ++l) {
        polarJz += debye_Jz[k-1][i][j]->Jpol[l]*(1.0 + _mat->km[l]);
    }

    polarJz += 0.5;

    eps_tmp += _mat->eps;//добавили дисперсионную поправку

    dhx = (fdtd_Hx[i][j][k-1] - fdtd_Hx[i][j-1][k-1])/_grd->dyi[j];
    dhy = (fdtd_Hy[j][k-1][i] - fdtd_Hy[j][k-1][i-1])/_grd->dxi[i];

    dEz = -fdtd_Ez[k-1][i][j];//-Ez[n]

    fdtd_Ez[k-1][i][j] = ((eps_tmp - 0.5*ht*sig_tmp)*fdtd_Ez[k-1][i][j]
                    + ht*_cl*(dhy - dhx) - ht*(_pi4*fdtd_Jz[k-1][i][j] + polarJz))/(eps_tmp + 0.5*ht*sig_tmp);

    dEz += fdtd_Ez[k-1][i][j];//Ez[n+1] - Ez[n]
    dEz /= ht;//(Ez[n+1] - Ez[n])/ht

    for (l = 0; l < _mat->nPoles; ++l) {
        debye_Jz[k-1][i][j]->Jpol[l] = _mat->km[l]*debye_Jz[k-1][i][j]->Jpol[l] + _mat->bm[l]*dEz;
    }
}
