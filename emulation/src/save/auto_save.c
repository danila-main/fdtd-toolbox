

/*! \file 
* \brief Функции созранения.
*
* Функции позволяют сохранить результаты расчета для последующего его возобновления.
* \author Крюков А.А. anton.krv@gmail.com
*/

/*! Отменяет предупреждения об опастности функций вывода. */
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../../include/globals.h"
#include "../../include/functions.h"
#include "../../include/fdtd.h"

static int transaction = 0;/*!< номер транзакции для корректного автосохранения. */

/*!
* Запись контольной точки на диск.
* \param filename - имя контрольной точки.
* \return void
*/
static void _auto_write(char *filename)
{
    int i, j, k;
    FILE *fp;
    grid _grd = &fdtd_grd;

    if ((fp = fopen(filename, "w+b")) == NULL) {
        fprintf(stderr, "Cannot open file %s\n", filename);
        exit(1);
    }

    for (i = 0; i < _grd->nx; ++i) {
        for (j = 0; j < _grd->ny+1; ++j) {
            fwrite(fdtd_Ex[i][j], sizeof(double), _grd->nz+1, fp);
            fwrite(fdtd_Jx[i][j], sizeof(double), _grd->nz+1, fp);
        }
    }

    for (j = 0; j < _grd->ny; ++j) {
        for (k = 0; k < _grd->nz+1; ++k) {
            fwrite(fdtd_Ey[j][k], sizeof(double), _grd->nx+1, fp);
            fwrite(fdtd_Jy[j][k], sizeof(double), _grd->nx+1, fp);
        }
    }

    for (k = 0; k < _grd->nz; ++k) {
        for(i = 0; i < _grd->nx+1; ++i) {
            fwrite(fdtd_Ez[k][i], sizeof(double), _grd->ny+1, fp);
            fwrite(fdtd_Jz[k][i], sizeof(double), _grd->ny+1, fp);
        }
    }

    for (i = 0; i < _grd->nx+1; ++i) {
        for (j = 0; j < _grd->ny; ++j) {
            fwrite(fdtd_Hx[i][j], sizeof(double), _grd->nz, fp);
        }
    }

    for (j = 0; j < _grd->ny+1; ++j) {
        for (k = 0; k < _grd->nz; ++k) {
            fwrite(fdtd_Hy[j][k], sizeof(double), _grd->nx, fp);
        }
    }

    for(k = 0; k < _grd->nz+1; ++k) {
        for(i = 0; i < _grd->nx; ++i) {
            fwrite(fdtd_Hz[k][i], sizeof(double), _grd->ny, fp);
        }
    }

    fclose(fp);
}

/*!
* Считывание контрольной точки с диска.
* \param filename - имя контрольной точки.
* \return void
*/
static void _auto_read(char *filename)
{
    int i, j, k;
    FILE *fp;
    grid _grd = &fdtd_grd;

    if ((fp = fopen(filename, "r+b")) == NULL) {
        fprintf(stderr, "Cannot open file %s\n", filename);
        exit(1);
    }

    for (i = 0; i < _grd->nx; ++i) {
        for (j = 0; j < _grd->ny+1; ++j) {
            fread(fdtd_Ex[i][j], sizeof(double), _grd->nz+1, fp);
            fread(fdtd_Jx[i][j], sizeof(double), _grd->nz+1, fp);
        }
    }

    for (j = 0; j < _grd->ny; ++j) {
        for (k = 0; k < _grd->nz+1; ++k) {
            fread(fdtd_Ey[j][k], sizeof(double), _grd->nx+1, fp);
            fread(fdtd_Jy[j][k], sizeof(double), _grd->nx+1, fp);
        }
    }

    for (k = 0; k < _grd->nz; ++k) {
        for(i = 0; i < _grd->nx+1; ++i) {
            fread(fdtd_Ez[k][i], sizeof(double), _grd->ny+1, fp);
            fread(fdtd_Jz[k][i], sizeof(double), _grd->ny+1, fp);
        }
    }

    for (i = 0; i < _grd->nx+1; ++i) {
        for (j = 0; j < _grd->ny; ++j) {
            fread(fdtd_Hx[i][j], sizeof(double), _grd->nz, fp);
        }
    }

    for (j = 0; j < _grd->ny+1; ++j) {
        for (k = 0; k < _grd->nz; ++k) {
            fread(fdtd_Hy[j][k], sizeof(double), _grd->nx, fp);
        }
    }

    for(k = 0; k < _grd->nz+1; ++k) {
        for(i = 0; i < _grd->nx; ++i) {
            fread(fdtd_Hz[k][i], sizeof(double), _grd->ny, fp);
        }
    }

    fclose(fp);
}

/*!
* Создание контольной точки.
* \param filename - имя контрольной точки.
* \param nt - номер шага по времени.
* \return void
*/
void make_save_point(char *filename, int nt)
{
    char tmp_name[256];
    char save_name[256];
    char tran_name[256];
    FILE *fp;

    create_name_with_id(filename, save_name, nt);

    strcpy(tmp_name, fdtd_start.save_dir);
    strcat(tmp_name, save_name);

    _auto_write(tmp_name);

    //пишем файл транзакции
    strcpy(tran_name, save_name);
    strcat(tran_name, ".SAV");

    strcpy(tmp_name, fdtd_start.save_dir);
    strcat(tmp_name, tran_name);

    if ((fp = fopen(tmp_name, "w+t")) == NULL) {
        fprintf(stderr, "Cannot open file %s\n", tmp_name);
        exit(1);
    }

    fprintf(fp, "%s\n", save_name);
    fprintf(fp, "%d", nt);

    fclose(fp);
}

/*!
* Загрузка контольной точки.
* \param filename - имя контрольной точки.
* \param nt - указатель на номер шага по времени.
* \return void
*/
void use_save_point(char *filename, int *nt)
{
    char tmp_name[256];
    char save_name[256];
    FILE *fp;

    strcpy(tmp_name, fdtd_start.save_dir);
    strcat(tmp_name, filename);

    if ((fp = fopen(tmp_name, "r+t")) == NULL) {
        fprintf(stderr, "Cannot open file %s\n", tmp_name);
        exit(1);
    }

    fscanf(fp, "%s\n", save_name);
    fscanf(fp, "%d", nt);

    fclose(fp);

    strcpy(tmp_name, fdtd_start.save_dir);
    strcat(tmp_name, save_name);

    _auto_read(tmp_name);
}

/*!
* Создание контольной точки автосохранения.
* \param nt - номер шага по времени.
* \return void
*/
void make_auto_save_point(int nt)
{
    char tmp_name[256];
    char save_name[256];
    char tran_name[256];
    FILE *fp;

    if(transaction == 0)
        transaction = 1;
    else
        transaction = 0;

    create_name_with_id("AUTO_SAVE", save_name, transaction);

    strcpy(tmp_name, fdtd_start.auto_save_dir);
    strcat(tmp_name, save_name);

    _auto_write(tmp_name);

    //пишем файл транзакции
    strcpy(tran_name, "AUTO_SAVE");
    strcat(tran_name, ".SAV");

    strcpy(tmp_name, fdtd_start.auto_save_dir);
    strcat(tmp_name, tran_name);

    if ((fp = fopen(tmp_name, "w+t")) == NULL) {
        fprintf(stderr, "Cannot open file %s\n", tran_name);
        exit(1);
    }

    fprintf(fp, "%s\n", save_name);
    fprintf(fp, "%d\n", nt);
    fprintf(fp, "%d", transaction);

    fclose(fp);
}

/*!
* Загрузка контольной точки автосохранения.
* \param nt - указатель на номер шага по времени.
* \return void
*/
int use_auto_save_point(int *nt)
{
    char tmp_name[256];
    char save_name[256];
    FILE *fp;

    strcpy(tmp_name, fdtd_start.auto_save_dir);
    strcat(tmp_name, "AUTO_SAVE.SAV");

    fp = fopen(tmp_name, "r+t");

    if(fp) {
        fscanf(fp, "%s\n", save_name);
        fscanf(fp, "%d", nt);
        fscanf(fp, "%d", &transaction);

        fclose(fp);

        strcpy(tmp_name, fdtd_start.auto_save_dir);
        strcat(tmp_name, save_name);

        _auto_read(tmp_name);

        return 1;
    }

    return 0;
}
