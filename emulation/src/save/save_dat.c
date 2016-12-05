
/*! \file 
* \brief Функции для сохранения результатов.
*
* Функции открытия потока вывода результатов расчета, записи результатов и закрытия потока вывода.
* \author Крюков А.А. anton.krv@gmail.com
*/

/*! Отменяет предупреждения об опастности функций вывода. */
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../../include/globals.h"
#include "../../include/constants.h"
#include "../../include/functions.h"

extern void save_data_vtr(const char *file_name, const double time,
                          const int nx, const double *xi,
                          const int ny, const double *yi,
                          const int nz, const double *zi);

double getEx(int i, int j, int k)
{
    grid _grd = &fdtd_grd;
    double sum = 0.0;

    sum += fdtd_Ex[i-1][j-1][k-1]*_grd->dyi[j-1]*_grd->dzi[k-1];
    sum += fdtd_Ex[i-1][j-1][k]*_grd->dyi[j-1]*_grd->dzi[k];
    sum += fdtd_Ex[i-1][j][k-1]*_grd->dyi[j]*_grd->dzi[k-1];
    sum += fdtd_Ex[i-1][j][k]*_grd->dyi[j]*_grd->dzi[k];

    return (0.25*sum/(_grd->dyi05[j]*_grd->dzi05[k]));
}

double getEy(int i, int j, int k)
{
    grid _grd = &fdtd_grd;
    double sum = 0.0;

    sum += fdtd_Ey[j-1][k-1][i-1]*_grd->dzi[k-1]*_grd->dxi[i-1];
    sum += fdtd_Ey[j-1][k-1][i]*_grd->dzi[k-1]*_grd->dxi[i];
    sum += fdtd_Ey[j-1][k][i-1]*_grd->dzi[k]*_grd->dxi[i-1];
    sum += fdtd_Ey[j-1][k][i]*_grd->dzi[k]*_grd->dxi[i];

    return (0.25*sum/(_grd->dzi05[k]*_grd->dxi05[i]));
}

double getEz(int i, int j, int k)
{
    grid _grd = &fdtd_grd;
    double sum = 0.0;

    sum += fdtd_Ez[k-1][i-1][j-1]*_grd->dxi[i-1]*_grd->dyi[j-1];
    sum += fdtd_Ez[k-1][i-1][j]*_grd->dxi[i-1]*_grd->dyi[j];
    sum += fdtd_Ez[k-1][i][j-1]*_grd->dxi[i]*_grd->dyi[j-1];
    sum += fdtd_Ez[k-1][i][j]*_grd->dxi[i]*_grd->dyi[j];

    return (0.25*sum/(_grd->dxi05[i]*_grd->dyi05[j]));
}

double getJx(int i, int j, int k)
{
    grid _grd = &fdtd_grd;
    double sum = 0.0;

    sum += fdtd_Jx[i-1][j-1][k-1]*_grd->dyi[j-1]*_grd->dzi[k-1];
    sum += fdtd_Jx[i-1][j-1][k]*_grd->dyi[j-1]*_grd->dzi[k];
    sum += fdtd_Jx[i-1][j][k-1]*_grd->dyi[j]*_grd->dzi[k-1];
    sum += fdtd_Jx[i-1][j][k]*_grd->dyi[j]*_grd->dzi[k];

    return (0.25*sum/(_grd->dyi05[j]*_grd->dzi05[k]));
}

double getJy(int i, int j, int k)
{
    grid _grd = &fdtd_grd;
    double sum = 0.0;

    sum += fdtd_Jy[j-1][k-1][i-1]*_grd->dzi[k-1]*_grd->dxi[i-1];
    sum += fdtd_Jy[j-1][k-1][i]*_grd->dzi[k-1]*_grd->dxi[i];
    sum += fdtd_Jy[j-1][k][i-1]*_grd->dzi[k]*_grd->dxi[i-1];
    sum += fdtd_Jy[j-1][k][i]*_grd->dzi[k]*_grd->dxi[i];

    return (0.25*sum/(_grd->dzi05[k]*_grd->dxi05[i]));
}

double getJz(int i, int j, int k)
{
    grid _grd = &fdtd_grd;
    double sum = 0.0;

    sum += fdtd_Jz[k-1][i-1][j-1]*_grd->dxi[i-1]*_grd->dyi[j-1];
    sum += fdtd_Jz[k-1][i-1][j]*_grd->dxi[i-1]*_grd->dyi[j];
    sum += fdtd_Jz[k-1][i][j-1]*_grd->dxi[i]*_grd->dyi[j-1];
    sum += fdtd_Jz[k-1][i][j]*_grd->dxi[i]*_grd->dyi[j];

    return (0.25*sum/(_grd->dxi05[i]*_grd->dyi05[j]));
}

double getHx(int i, int j, int k)
{
    grid _grd = &fdtd_grd;
    double sum = 0.0;

    sum += fdtd_Hx[i-1][j-1][k-1]*_grd->dxi[i-1];
    sum += fdtd_Hx[i][j-1][k-1]*_grd->dxi[i];

    return (0.5*sum/(_grd->dxi05[i]));
}

double getHy(int i, int j, int k)
{
    grid _grd = &fdtd_grd;
    double sum = 0.0;

    sum += fdtd_Hy[j-1][k-1][i-1]*_grd->dyi[j-1];
    sum += fdtd_Hy[j][k-1][i-1]*_grd->dyi[j];

    return (0.5*sum/(_grd->dyi05[j]));
}

double getHz(int i, int j, int k)
{
    grid _grd = &fdtd_grd;
    double sum = 0.0;

    sum += fdtd_Hz[k-1][i-1][j-1]*_grd->dzi[k-1];
    sum += fdtd_Hz[k][i-1][j-1]*_grd->dzi[k];

    return (0.5*sum/(_grd->dzi05[k]));
}

double getE(int i, int j, int k)
{
    double val;
    double val_x, val_y, val_z;

    val_x = getEx(i, j, k);

    val_y = getEy(i, j, k);

    val_z = getEz(i, j, k);

    val = sqrt(val_x*val_x + val_y*val_y + val_z*val_z);

    return val;
}

double getJ(int i, int j, int k)
{
    double val;
    double val_x, val_y, val_z;

    val_x = getJx(i, j, k);

    val_y = getJy(i, j, k);

    val_z = getJz(i, j, k);

    val = sqrt(val_x*val_x + val_y*val_y + val_z*val_z);

    return val;
}

double getH(int i, int j, int k)
{
    double val;
    double val_x, val_y, val_z;

    val_x = getHx(i, j, k);

    val_y = getHy(i, j, k);

    val_z = getHz(i, j, k);

    val = sqrt(val_x*val_x + val_y*val_y + val_z*val_z);

    return val;
}

void save_dat(int it)
{
    char file_name[256];
    grid _grd = &fdtd_grd;

    sprintf(file_name, "%sstep_%d.vtr", fdtd_start.result_dir, it);

    save_data_vtr(file_name, _grd->ti[it],
                  _grd->nx, _grd->xi,
                  _grd->ny, _grd->yi,
                  _grd->nz, _grd->zi);
}

/*!
* Создание имени с числовым идентификатором (расширения вида .xxx игнорируются) name.xxx -> name_id.xxx.
* \param name - имя.
* \param new_name - новое имя.
* \param id - идентификатор.
* \return void
*/
void create_name_with_id(char *name, char *new_name, int id)
{
    char extension[16];
    char str_id[16];
    char *dot_pos = 0;

    sprintf(str_id, "%d", id);

    strcpy(new_name, name);
    dot_pos = strrchr(new_name, '.');

    if(dot_pos)
    {
        strcpy(extension, dot_pos);
        *dot_pos = '\0';
        strcat(new_name, "_");
        strcat(new_name, str_id);
        strcat(new_name, extension);
    }
    else
    {
        strcat(new_name, "_");
        strcat(new_name, str_id);
    }
}
