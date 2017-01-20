

/*! \file 
* \brief Определение информации о разбиении области.
*
* Функция формирует список соседних процессов, тип граничных условий, выделяет память для синхронизации между процесссами.
* \author Крюков А.А. anton.krv@gmail.com
*/

#include <stdlib.h>
#include <math.h>

#include "../../include/globals.h"
#include "../../include/data.h"

/*!
* Функция получения информации о разбиении области.
*/
void get_local()
{
    comp _comp = &spars_comp;
    grid _grd = &spars_grd;

    //вычисляем диапазон ячеек
    _pos->nz_l = 1;
    _pos->nz_r = _grd->nz;
    _pos->num_z = _pos->nz_r - _pos->nz_l + 1;

    _pos->ny_l = 1;
    _pos->ny_r = _grd->ny;
    _pos->num_y = _pos->ny_r - _pos->ny_l + 1;

    _pos->nx_l = 1;
    _pos->nx_r = _grd->nx;
    _pos->num_x = _pos->nx_r - _pos->nx_l + 1;
}
