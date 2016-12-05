
/*! \file 
* \brief Описание структур и макросов.
*
* Данные структуры предназначены для хранения матрицы, векторов, предобусловливателя,
* параметров итерационного решателя и предобусловливателя,
* данных о взаимодействии между процессами.
* \author Крюков А.А. anton.krv@gmail.com
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <memory.h>

#ifndef _DATA_HEADER_
#define _DATA_HEADER_

/*! Сообщение об ошибке. */
#define CHECKERR(ierr) assert(!(ierr))

/*! Оператор выделения памяти */
#define PARMS_malloc(base, nmem, type) {\
    (base) = (type *)malloc((nmem)*sizeof(type)); \
    CHECKERR((base) == NULL); \
}

#define TRUE      1 /*!< логическая истина. */
#define FALSE     0 /*!< логическая ложь. */

#endif

