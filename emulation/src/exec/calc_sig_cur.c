

/*! \file 
* \brief Функции расчета внешних токов, проводимости и т.п.
*
* Функции заполняют массивы токов, проводимости, вторичных електронов, концентрации ионов.
* Для расчета используются физические модели из файла src/user/user.c
* \author Крюков А.А. anton.krv@gmail.com
*/

#include <stdlib.h>
#include <math.h>

#include "../../include/globals.h"
#include "../../include/constants.h"
#include "../../include/functions.h"
#include "../../include/fdtd.h"

/*!
* \brief Заполнение массива проводимости.
*
* Проводимость вычисляется только в материале где разрешена соотвествующая модель проводимости.
* Кроме этого заполняются массивы плотности вторичных электронов, концентрации ионов.
* \param time - время.
* \return void
*/
void initSig()
{
    int i, j, k;
    grid _grd = &fdtd_grd;
    layers _lay = &fdtd_lay;
    int i_lay; //номер слоя
    double x, y, z;

    for (i = 0; i <= _grd->nx; ++i) {
        for (j = 0; j <= _grd->ny; ++j) {
            for (k = 0; k <= _grd->nz; ++k) {
                i_lay = fdtd_gi[i][j][k];
                x = _grd->xi05[i];
                y = _grd->yi05[j];
                z = _grd->zi05[k];
                fdtd_Sig[i][j][k] = _lay->data_lay[i_lay].lay_sig;
            }
        }
    }
}

#define TEST_FREQ 5.0e+08 //частота в герцах

void setCurr(double time)
{
    int i, j, k;
    double x, y, z;
    grid _grd = &fdtd_grd;

    double wave_num = (2.0*_pi*TEST_FREQ)/_cl;
    double wave_len = _cl/TEST_FREQ;
    double l_antenna = wave_len;
    double I0 = 1.0e+10;

    return;

    for (k = 0; k < _grd->nz; ++k) {
        for (i = 0; i < _grd->nx+1; ++i) {
            for (j = 0; j < _grd->ny+1; ++j) {
                x = _grd->xi[i];
                y = _grd->yi[j];
                z = _grd->zi05[k+1];

                if (fabs(x) <= 1.0e-10 && fabs(y) <= 1.0e-10) {
                    if (z >= 0.0 && z <= 0.5*l_antenna) {
                        fdtd_Jz[k][i][j] = I0*sin(wave_num*(0.5*l_antenna - z));
                        //fdtd_Jz[k][i][j] = I0*(1.0 - 2.0*z / l_antenna);


                        fdtd_Jz[k][i][j] *= sin((2.0*_pi*TEST_FREQ)*time);
                    } else if (z <= 0.0 && z >= -0.5*l_antenna) {
                        fdtd_Jz[k][i][j] = I0*sin(wave_num*(0.5*l_antenna + z));
                        //fdtd_Jz[k][i][j] = I0*(1.0 + 2.0*z / l_antenna);


                        fdtd_Jz[k][i][j] *= sin((2.0*_pi*TEST_FREQ)*time);
                    } else {
                        fdtd_Jz[k][i][j] = 0.0;
                    }
                }
            }
        }
    }
}