
/*! \file 
* \brief Основная программа.
*
* Программа предназначена для запуска в среде с поддержкой MPI.
* Программа производит инициализацию многопроцессорной MPI среды, считывает параметры расчета,
* на основании этих данных формирует геометрическое разбиение расчетной области.
* В зависимости от параметров осуществляет сохранение и(или) автосохранение, для возможности возобновлениея расчета.
* После инициализации и обработки данных проета запускается основной расчетный цикл.
* По его завершению производится освобождение выделенной памяти
* \author Крюков А.А. anton.krv@gmail.com
*/

/*! Отменяет предупреждения об опастности функций вывода. */
#define _CRT_SECURE_NO_WARNINGS

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../../include/globals.h"
#include "../../include/functions.h"
#include "../../include/fdtd.h"
#include "../../include/constants.h"

/*! Имя файла содержащего имя проекта и параметры запуска программы. */
#define START_NAME "START_N"

/*! Вычисление решения для верхнего временного слоя. */
extern void solverFDTD(int n, double time);

extern void zeroCurrents();

/*! Вычисление проводимости. */
extern void initSig();

extern double dwalltime();

extern void init_data_pvd();

extern void setCurr(double time);

extern void init_data_nfff(const double time);

/*!
 * Главная функция программы.
 * \param argc - число аргументов строки запуска, (необходимо для инициализации среды MPI).
 * \param argv - строка запуска, (необходима для инициализации среды MPI).
 * \return 0, успешное завершение.
 */
int main(int argc, char* argv[])
{
    int n;
    int n0;
    char log_msg[256];
    double time;//время исполнения
    project _prj = &fdtd_prj;

    n0 = 0;

    open_log_file("protocol.log");

    readStart(START_NAME);//считываем параметры запуска

    readProject();//считываем параметры проекта
    print_log_message("Project DONE");

    readTok();//считываем параметры расчета
    print_log_message("Calculation DONE");

    readGrd();//считываем сетку
    print_log_message("Grid  DONE");

    allocateCells();//выделение памяти для полей
    print_log_message("Fields memory has allocated");

    readLay();//считываем слои
    print_log_message("Layers data has read");

    readGi();//считываем привязку к модели
    print_log_message("Medium data has read");

    initSig();//задание проводимости из материалов
    print_log_message("Conductivity has synced");

    if (strcmp(_prj->start_name, "NONE")) {//используем точку сохранения
        use_save_point(_prj->start_name, &n0);
        print_log_message("We have a save point!");
        sprintf(log_msg, "Start time step is %d", n0);
        print_log_message(log_msg);
    } else {
        print_log_message("We don't have a save point");
    }

    if (use_auto_save_point(&n0)) {//проверяем автосохранение
        print_log_message("We have an auto save point!");
        print_log_message("Save points are ignored in this case!");
        sprintf(log_msg, "Start time step is %d", n0);
        print_log_message(log_msg);
    } else {
        print_log_message("We don't have an auto save point");
    }

    print_log_message("All output streams are adjusted to the start time point");

    init_data_pvd();

    if (fdtd_comp.l_nfff) {
        init_data_nfff(fdtd_grd.ti[0]);
        write_nfff_grid();
    }

    for (n = n0; n < fdtd_grd.nt; ++n) {//основной цикл по времени
        printf("Step %d\n", n);
        sprintf(log_msg, "Step %d starts", n);
        print_log_message(log_msg);

        time = dwalltime();

        if (n%10 == 0) {
            save_dat(n);
        }

        time = dwalltime() - time;

        sprintf(log_msg, "Save data: step is %d, time %f", n, time);
        print_log_message(log_msg);

        time = dwalltime();

        if (_prj->cur_it_save < _prj->nt_save && _prj->it_save[_prj->cur_it_save] == n) {//сохранение
            make_save_point("SAVE_POINT", n);

            ++(_prj->cur_it_save);

            time = dwalltime() - time;

            sprintf(log_msg, "Made a save point: step is %d, time %f", n, time);
            print_log_message(log_msg);

            time = dwalltime();
        }

        if(fdtd_start.auto_save != 0 && n%fdtd_start.auto_save == 0) {//автосохранение
            make_auto_save_point(n);

            time = dwalltime() - time;

            sprintf(log_msg, "Made an auto save point: step is %d, time %f", n, time);
            print_log_message(log_msg);

            time = dwalltime();
        }

        //вычисление шага по времени
        fdtd_ht = fdtd_grd.ti[n+1] - fdtd_grd.ti[n];
        fdtd_ht05 = fdtd_ht;

        if (n+1 < fdtd_grd.nt && n >= 1)
            fdtd_ht05 = 0.5*(fdtd_grd.ti[n+1] - fdtd_grd.ti[n-1]);

        fdtd_time = fdtd_grd.ti[n+1];

        zeroCurrents();//обнуление токов
        setCurr(fdtd_grd.ti[n+1]);//задание токов

        //запуск решателей
        if (fdtd_comp.l_shem == 1) {
            int loc_n_step;//число локальных шагов с учетом курранта
            int n_loc;
            double ht = fdtd_ht;//сохраняем исходные шаги по времени
            double ht05 = fdtd_ht05;

            loc_n_step = (int)ceil(fdtd_ht/fdtd_max_ht);//делим на максимальный шаг по времени и приводим к наибольшему целому

            fdtd_ht /= loc_n_step;//делим на число локальных шагов
            fdtd_ht05 /= loc_n_step;

            for (n_loc = 0;n_loc < loc_n_step;n_loc++) {
                solverFDTD(n+1, fdtd_grd.ti[n+1]);
            }

            fdtd_ht = ht;//восстанавливаем исходные шаги по времени
            fdtd_ht05 = ht05;
        }

        time = dwalltime() - time;

        sprintf(log_msg, "Solver at step %d has completed, time %f", n, time);
        print_log_message(log_msg);

        sprintf(log_msg, "Step %d ends off", n);
        print_log_message(log_msg);
    }

    close_log_file();//завершили запись протокола

    return 0;
}
