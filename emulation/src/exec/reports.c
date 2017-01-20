
/*! \file 
* \brief Функции вывода протокола работы программы.
*
* Данные функции используются для вывода протокола работы программы.
* Протокол предназначен для диагностики процесса выполнения программы.
* Все процессы выполняют одинаковую работу.
* Все функции требующие обменов неявно синхронизируются.
* Нет причины делать протокол, основанный на широковещательном опросе всех процессов.
* Протокол по умолчанию выводится для нулевого.
* \author Крюков А.А. anton.krv@gmail.com
*/

/*! Отменяет предупреждения об опастности функций вывода. */
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../include/globals.h"

/*! Поток вывода информации. */
static FILE *_fp_log = NULL;

/*!
* Создание потока вывода информации.
* \param log_name - имя протокола.
* \return void
*/
void open_log_file(char *log_name)//открывает поток вывода протокола
{
    if ((_fp_log = fopen(log_name, "w+b")) == NULL) {
        fprintf(stderr, "Cannot open file %s\n", log_name);
        exit(1);
    }
}

/*!
* Закрытие потока вывода информации.
* \return void
*/
void close_log_file()//закрывает потоки вывода протокола
{
    if (_fp_log) {
        fclose(_fp_log);
    }
}

/*!
* Печеть сообщения в протокол.
* \param msg - текст сообщения.
* \return void
*/
void print_log_message(char *msg)//печатает сообщение
{
    fprintf(_fp_log, "%s\n", msg);

    fflush(_fp_log);
}
