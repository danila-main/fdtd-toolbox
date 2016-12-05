
/*! \file 
* \brief Функции для работы с проектом и разностной схемой.
* 
* Функции для работы с проектом осуществляют считывание и обработку данных проекта, вывод результатов и запись/использование контрольных точек.
* Функции для работы с разностной схемой осуществляют запись разностного шаблона в матрицу,
* изменение элементов шаблона, вычисление индекса элемента, вычисление правой части, синхронизацию полей между процессами.
* \author Крюков А.А. anton.krv@gmail.com
*/

#ifndef _SETPAR_HEADER_
#define _SETPAR_HEADER_

#include "data.h"
#include "data_types.h"

/*! Считывание параметров запуска программы. */
extern void readStart(char *filename);
/*! Считывание параметров проекта. */
extern void readProject();
/*! Считывание параметров расчета. */
extern void readTok();
/*! Считывание сетки. */
extern void readGrd();
/*! Считывание параметров материалов. */
extern void readLay();
/*! Считывание привязки материалов к сетке. */
extern void readGi();
/*! Считывание параметров вывода результатов. */
extern void readVid();
/*! Сохранение результатов. */
extern void save_dat(int it);
/*! Выделение памяти для хранения полей, массивов проводимости и т.п. */
extern void allocateCells();
/*! Проводимость в точке, усредненная по направлению. */
extern double fun_sig(int i, int j, int k, fieldType idf);
/*! Диэлектрическая проницаемость в точке, усредненная по направлению. */
extern double fun_eps(int i, int j, int k, fieldType idf);
/*! Магнитная проницаемость в точке, усредненная по направлению. */
extern double fun_mu(int i, int j, int k, fieldType idf);
/*! Создание имени с числовым идентификатором. */
extern void create_name_with_id(char *name, char *new_name, int id);

/*! Создание точки сохранения. */
extern void make_save_point(char *filename, int nt);
/*! Создание точки автоматического сохранения. */
extern void make_auto_save_point(int nt);
/*! Загрузка точки сохранения. */
extern void use_save_point(char *filename, int *nt);
/*! Загрузка точки автоматического сохранения. */
extern int use_auto_save_point(int *nt);

/*! Открытие файла протокола. */
extern void open_log_file(char *log_name);
/*! Закрытие файла протокола. */
extern void close_log_file();
/*! Запись сообщения в протокол. */
extern void print_log_message(char *msg);

/*! Системное время в секундах. */
extern double dwalltime();

extern double getEx(int i, int j, int k);
extern double getEy(int i, int j, int k);
extern double getEz(int i, int j, int k);
extern double getJx(int i, int j, int k);
extern double getJy(int i, int j, int k);
extern double getJz(int i, int j, int k);
extern double getHx(int i, int j, int k);
extern double getHy(int i, int j, int k);
extern double getHz(int i, int j, int k);

extern double getE(int i, int j, int k);
extern double getJ(int i, int j, int k);
extern double getH(int i, int j, int k);

#endif
