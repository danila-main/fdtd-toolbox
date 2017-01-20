
/*! \file 
* \brief ќбъ¤влени¤ глобальных структур данных, используемых внутри программы.
*
* ƒл¤ использовани¤ глобальных данных внутри единицы компил¤ции необходимо подключить это файл.
* \author  рюков ј.ј. anton.krv@gmail.com
*/

#ifndef _GLOBALS_HEADER_
#define _GLOBALS_HEADER_

#include "data_types.h"

extern double fdtd_ht;/*!< шаг по времени. */
extern double fdtd_ht05;/*!< шаг по времени между полуцелыми шагами. */
extern double fdtd_max_ht;/*!< максимальный шаг по времени. */
extern double fdtd_time;/*!< время. */

extern _start fdtd_start;/*!< параметры запуска программы. */
extern _project fdtd_prj;/*!< параметры проекта. */
extern _comp fdtd_comp;/*!< параметры расчета. */
extern _grid fdtd_grd;/*!< параметры разбиения и разностной сетки. */
extern _layers fdtd_lay;/*!< параметры материалов. */

extern double ***fdtd_Ex;/*!< массив x-компоненты электрического пол¤. */
extern double ***fdtd_Ey;/*!< массив y- компоненты электрического пол¤. */
extern double ***fdtd_Ez;/*!< массив z- компоненты электрического пол¤. */
extern double ***fdtd_Hx;/*!< массив x- компоненты магнитного пол¤. */
extern double ***fdtd_Hy;/*!< массив y- компоненты магнитного пол¤. */
extern double ***fdtd_Hz;/*!< массив z- компоненты магнитного пол¤. */
extern double ***fdtd_Jx;/*!< массив x- компоненты стороннего тока. */
extern double ***fdtd_Jy;/*!< массив y- компоненты стороннего тока. */
extern double ***fdtd_Jz;/*!< массив z- компоненты стороннего тока. */

extern double ***fdtd_Sig;/*!< массив проводимости. */
extern int ***fdtd_gi;/*!< массив номеров материалов. */

#endif
