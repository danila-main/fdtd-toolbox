

/*! \file 
* \brief Основная программа.
*
* 
* \author Крюков А.А. anton.krv@gmail.com
*/

/*! Отменяет предупреждения об опастности функций вывода. */
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include "../../include/globals.h"
#include "../../include/data.h"
#include "../../include/constants.h"

/*!
* Чтение файла с именем проекта и параметрами запуска программы.
* \param filename - имя файла с данными.
* \return void
*/
void readStart(char *filename)
{
    FILE *fp;

    if ((fp = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "Cannot open file %s\n", filename);
        exit(1);
    }

    fscanf(fp,"%s", fdtd_start.prj_name);
    fscanf(fp,"%s", fdtd_start.result_dir);
    fscanf(fp,"%s", fdtd_start.save_dir);
    fscanf(fp,"%s", fdtd_start.auto_save_dir);
    fscanf(fp,"%d", &(fdtd_start.auto_save));

    fclose(fp);
}

/*!
* Чтение данных проекта.
*/
void readProject()
{
    FILE *fp;
    char buf[256];
    char key_word[64];
    int i;
    int n_out;
    _name64 vid_name;
    _name64 dat_name;

    if ((fp = fopen(fdtd_start.prj_name, "r")) == NULL) {
        fprintf(stderr, "Cannot open file %s\n", fdtd_start.prj_name);
        exit(1);
    }

    key_word[0] = '\0';

    while(fgets(buf, 256, fp) != NULL) {
       sscanf(buf, "%*[ ]<%[^\n>]", key_word);//считываем ключевое слово

        if (!strcmp(key_word, "Grd name")) {//имя сетки
            fgets(buf, 256, fp);//переход на след. строку
            sscanf(buf,"%s", fdtd_prj.grd_name);
        } else if (!strcmp(key_word, "Cell name")) {//имя файла привязки
            fgets(buf, 256, fp);//переход на след. строку
            sscanf(buf,"%s", fdtd_prj.cel_name);
        } else if (!strcmp(key_word, "Tok name")) {//имя файла параметров расчета
            fgets(buf, 256, fp);//переход на след. строку
            sscanf(buf,"%s", fdtd_prj.tok_name);
        } else if (!strcmp(key_word, "Layers name")) {//имя файла параметров слоев
            fgets(buf, 256, fp);//переход на след. строку
            sscanf(buf,"%s", fdtd_prj.lay_name);
        } else if (!strcmp(key_word, "Output")) {//имя файлов вывода результатов
            fgets(buf, 256, fp);//переход на след. строку
            fdtd_prj.n_output = 0;

            sscanf(buf,"%d", &(n_out));

            if (n_out > 0) {
                PARMS_malloc(fdtd_prj.vid_names, n_out, _name64);
                PARMS_malloc(fdtd_prj.dat_names, n_out, _name64);
                PARMS_malloc(fdtd_prj.out_data, n_out, _output);

                for (i = 0; i < n_out; ++i) {
                    fgets(buf, 256, fp);
                    sscanf(buf,"%s", vid_name.data);

                    fgets(buf, 256, fp);
                    sscanf(buf,"%s", dat_name.data);

                    strcpy(fdtd_prj.vid_names[fdtd_prj.n_output].data, vid_name.data);
                    strcpy(fdtd_prj.dat_names[fdtd_prj.n_output].data, dat_name.data);
                    fdtd_prj.n_output++;
                }
            }
        } else if (!strcmp(key_word, "Recent calculation")) {//имя контрольной точки
            fgets(buf, 256, fp);//переход на след. строку
            sscanf(buf,"%s", fdtd_prj.start_name);
        }
    }

    fclose(fp);
}

/*!
* Чтение данных расчета.
*/
void readTok()
{
    FILE *fp;
    char buf[256];
    char tmp[32];
    int num;

    if ((fp = fopen(fdtd_prj.tok_name, "r")) == NULL) {
        fprintf(stderr, "Cannot open file %s\n", fdtd_prj.tok_name);
        exit(1);
    }

    num = 1;
    while (fgets(buf, 256, fp) != NULL) {
        switch(num) {
            case 3:
                sscanf(buf,"%d", &(fdtd_comp.l_eh0));
            break;
            case 5:
                sscanf(buf,"%d", &(fdtd_comp.l_shem));
            break;

            case 7:
                sscanf(buf,"%d", &(fdtd_comp.l_time));
            break;

            case 9:
                sscanf(buf,"%d", &(fdtd_comp.l_eh_ext));
            break;

            case 11:
                sscanf(buf,"%d", &(fdtd_comp.l_gran));
            break;

            case 13:
                sscanf(buf,"%lf%*[ ]%lf%*[ ]%lf", &(fdtd_comp.front_x), &(fdtd_comp.front_y), &(fdtd_comp.front_z));
            break;

            case 15:
                sscanf(buf,"%s%*[ ]%d%*[ ]%d", tmp, &(fdtd_comp.xmin_bnd), &(fdtd_comp.xmax_bnd));
            break;

            case 16:
                sscanf(buf,"%s%*[ ]%d%*[ ]%d", tmp, &(fdtd_comp.ymin_bnd), &(fdtd_comp.ymax_bnd));
            break;

            case 17:
                sscanf(buf,"%s%*[ ]%d%*[ ]%d", tmp, &(fdtd_comp.zmin_bnd), &(fdtd_comp.zmax_bnd));
            break;

            case 21:
                sscanf(buf,"%d", &(fdtd_comp.l_part));
            break;

            case 23:
                sscanf(buf,"%d", &(fdtd_comp.l_j));
            break;

            case 25:
                sscanf(buf,"%d", &(fdtd_comp.l_qc));
            break;

            case 27:
                sscanf(buf,"%d", &(fdtd_comp.l_sec));
            break;
        }
        num++;
    }

    fclose(fp);

    fdtd_comp.l_tfsf = TRUE;//включена плоская волна
    fdtd_comp.l_nfff = FALSE;//включен расчет диаграммы направленности
    fdtd_comp.l_disp = TRUE;//наличие дисперсионного материала

    //0 Идеальный электрический проводник
    //1 Идеальный магнитный проводник
    //2 Условия Мура 2-го порядка
    //3 CPML

    //обработка граничных условий
    switch (fdtd_comp.l_gran) {
    case BND_ZERO_E ://Идеальный электрический проводник
        fdtd_comp.xmin_bnd = BND_ZERO_E;
        fdtd_comp.xmax_bnd = BND_ZERO_E;
        fdtd_comp.ymin_bnd = BND_ZERO_E;
        fdtd_comp.ymax_bnd = BND_ZERO_E;
        fdtd_comp.zmin_bnd = BND_ZERO_E;
        fdtd_comp.zmax_bnd = BND_ZERO_E;

        break;

    case BND_ZERO_H ://Идеальный магнитный проводник
        fdtd_comp.xmin_bnd = BND_ZERO_H;
        fdtd_comp.xmax_bnd = BND_ZERO_H;
        fdtd_comp.ymin_bnd = BND_ZERO_H;
        fdtd_comp.ymax_bnd = BND_ZERO_H;
        fdtd_comp.zmin_bnd = BND_ZERO_H;
        fdtd_comp.zmax_bnd = BND_ZERO_H;

        break;

    case BND_MUR_2 ://Условия Мура 2-го порядка
        fdtd_comp.xmin_bnd = BND_MUR_2;
        fdtd_comp.xmax_bnd = BND_MUR_2;
        fdtd_comp.ymin_bnd = BND_MUR_2;
        fdtd_comp.ymax_bnd = BND_MUR_2;
        fdtd_comp.zmin_bnd = BND_MUR_2;
        fdtd_comp.zmax_bnd = BND_MUR_2;

        break;

    case BND_CPML ://CPML
        fdtd_comp.xmin_bnd = BND_CPML;
        fdtd_comp.xmax_bnd = BND_CPML;
        fdtd_comp.ymin_bnd = BND_CPML;
        fdtd_comp.ymax_bnd = BND_CPML;
        fdtd_comp.zmin_bnd = BND_CPML;
        fdtd_comp.zmax_bnd = BND_CPML;

        break;
    }
}

/*!
* Считывание регулярной сетки по одному направлению.
* \param _grd - массив сетки.
* \param no_grd - число больших отрезков.
* \param n_grd - число точек.
* \param fp - поток ввода для файла сетки.
* \return void
*/
static void reading(double* _grd, int no_grd, int n_grd, FILE *fp)
{
    int ix, j, i;
    double dx_l;
    double dx_r;
    double c;
    int n;
    int l_dir;
    double q;
    char buf[256];

    fgets(buf, 256, fp);
    fgets(buf, 256, fp);

    ix = 0;
    sscanf(buf,"%lf", &(_grd[ix]));
    ix++;

    for(j=1;j<=no_grd;j++)
    {
        fgets(buf, 256, fp);//шаг слева
        sscanf(buf,"%lf", &dx_l);

        fgets(buf, 256, fp);//шаг справа
        sscanf(buf,"%lf", &dx_r);

        fgets(buf, 256, fp);//коэффициент
        sscanf(buf,"%lf", &c);

        fgets(buf, 256, fp);//количество точек
        sscanf(buf,"%d", &n);
        n--;

        fgets(buf, 256, fp);

        fgets(buf, 256, fp);//направление - 1 если слева направо, 0 если справа налево
        sscanf(buf,"%d", &l_dir);

        fgets(buf, 256, fp);

        fgets(buf, 256, fp);
        sscanf(buf,"%lf", &(_grd[ix+n-1]));

        if(l_dir==1)// Отрезок обрабатывается слева направо  -->
        {
            q = dx_l;
            for(i=1;i<=n-1;i++)
            {
                _grd[ix] = _grd[ix-1] + q;
                ix++;
                q = q*c;
            }
            ix++;
        }
        else//Отрезок обрабатывается справа налево  <--
        {
            q = dx_r;
            ix = ix + n - 2;
            for(i=n-1;i>=1;i--)
            {
                _grd[ix] = _grd[ix+1] - q;
                ix--;
                q = q*c;
            }
            ix = ix + n + 1;
        }
    }
}

/*!
* Считывание разностной сетки.

*/
void readGrd()
{
    FILE *fp;
    char buf[256];
    char tmp[32];
    int no_grd;
    int i;
    double min_hx, min_hy, min_hz;
    grid _grd = &fdtd_grd;

    if ((fp = fopen(fdtd_prj.grd_name, "r")) == NULL) {
        fprintf(stderr, "Cannot open file %s\n", fdtd_prj.grd_name);
        exit(1);
    }

    while (fgets(buf, 256, fp) != NULL) {
        sscanf(buf,"%s", tmp);

        if (!strcmp(tmp, "X")) {
            fgets(buf, 256, fp);
            sscanf(buf,"%s", tmp);
            if(!strcmp(tmp, "X"))
                break;
        }
    }

    fgets(buf, 256, fp);
    sscanf(buf, "%d", &no_grd);

    fgets(buf, 256, fp);
    sscanf(buf, "%d", &(_grd->nx));

    --(_grd->nx);

    PARMS_malloc(_grd->xi, (_grd->nx+1), double);
    PARMS_malloc(_grd->xi05, (_grd->nx+2), double);
    PARMS_malloc(_grd->dxi, (_grd->nx+1), double);
    PARMS_malloc(_grd->dxi05, (_grd->nx+2), double);

    reading(_grd->xi, no_grd, _grd->nx, fp);

    _grd->xi05[0] = 1.5*_grd->xi[0] - 0.5*_grd->xi[1];
    for(i=1;i<=_grd->nx;++i) {
        _grd->xi05[i] = 0.5*(_grd->xi[i] + _grd->xi[i-1]);
    }
    _grd->xi05[_grd->nx+1] = 1.5*_grd->xi[_grd->nx] - 0.5*_grd->xi[_grd->nx-1];


    for(i=0;i<=_grd->nx;++i) {
        _grd->dxi[i] = _grd->xi05[i+1] - _grd->xi05[i];
    }

    min_hx = DBL_MAX;
    _grd->dxi05[0] = 0.0;
    for(i=1;i<=_grd->nx;++i) {
        _grd->dxi05[i] = _grd->xi[i] - _grd->xi[i-1];
        if (_grd->dxi05[i] < min_hx) {
            min_hx = _grd->dxi05[i];
        }
    }
    _grd->dxi05[_grd->nx+1] = 0.0;

    fgets(buf, 256, fp);
    fgets(buf, 256, fp);

    fgets(buf, 256, fp);
    sscanf(buf, "%d", &no_grd);

    fgets(buf, 256, fp);
    sscanf(buf, "%d", &(_grd->ny));

    --(_grd->ny);

    PARMS_malloc(_grd->yi, (_grd->ny+1), double);
    PARMS_malloc(_grd->yi05, (_grd->ny+2), double);
    PARMS_malloc(_grd->dyi, (_grd->ny+1), double);
    PARMS_malloc(_grd->dyi05, (_grd->ny+2), double);

    reading(_grd->yi, no_grd, _grd->ny, fp);

    _grd->yi05[0] = 1.5*_grd->yi[0] - 0.5*_grd->yi[1];
    for(i=1;i<=_grd->ny;++i) {
        _grd->yi05[i] = 0.5*(_grd->yi[i] + _grd->yi[i-1]);
    }
    _grd->yi05[_grd->ny+1] = 1.5*_grd->yi[_grd->ny] - 0.5*_grd->yi[_grd->ny-1];


    for(i=0;i<=_grd->ny;++i) {
        _grd->dyi[i] = _grd->yi05[i+1] - _grd->yi05[i];
    }

    min_hy = DBL_MAX;
    _grd->dyi05[0] = 0.0;
    for(i=1;i<=_grd->ny;++i) {
        _grd->dyi05[i] = _grd->yi[i] - _grd->yi[i-1];
        if (_grd->dyi05[i] < min_hy) {
            min_hy = _grd->dyi05[i];
        }
    }
    _grd->dyi05[_grd->ny+1] = 0.0;

    fgets(buf, 256, fp);
    fgets(buf, 256, fp);

    fgets(buf, 256, fp);
    sscanf(buf, "%d", &no_grd);

    fgets(buf, 256, fp);
    sscanf(buf, "%d", &(_grd->nz));

    --(_grd->nz);

    PARMS_malloc(_grd->zi, (_grd->nz+1), double);
    PARMS_malloc(_grd->zi05, (_grd->nz+2), double);
    PARMS_malloc(_grd->dzi, (_grd->nz+1), double);
    PARMS_malloc(_grd->dzi05, (_grd->nz+2), double);

    reading(_grd->zi, no_grd, _grd->nz, fp);

    _grd->zi05[0] = 1.5*_grd->zi[0] - 0.5*_grd->zi[1];
    for(i=1;i<=_grd->nz;++i) {
        _grd->zi05[i] = 0.5*(_grd->zi[i] + _grd->zi[i-1]);
    }
    _grd->zi05[_grd->nz+1] = 1.5*_grd->zi[_grd->nz] - 0.5*_grd->zi[_grd->nz-1];

    for(i=0;i<=_grd->nz;++i) {
        _grd->dzi[i] = _grd->zi05[i+1] - _grd->zi05[i];
    }

    min_hz = DBL_MAX;
    _grd->dzi05[0] = 0.0;
    for(i=1;i<=_grd->nz;++i) {
        _grd->dzi05[i] = _grd->zi[i] - _grd->zi[i-1];
        if (_grd->dzi05[i] < min_hz)
            min_hz = _grd->dzi05[i];
    }
    _grd->dzi05[_grd->nz+1] = 0.0;

    fgets(buf, 256, fp);
    fgets(buf, 256, fp);

    fgets(buf, 256, fp);
    sscanf(buf, "%d", &no_grd);

    fgets(buf, 256, fp);
    sscanf(buf, "%d", &(_grd->nt));

    --(_grd->nt);

    PARMS_malloc(_grd->ti, (_grd->nt+1), double);

    reading(_grd->ti, no_grd, _grd->nt, fp);

    fdtd_max_ht = (1.0/_cl)*(1.0/sqrt(1.0/(min_hx*min_hx) + 1.0/(min_hy*min_hy) + 1.0/(min_hz*min_hz)));

    fclose(fp);
}

/*!
* Считывание параметров материалов.

*/
void readLay()
{
    FILE *fp;
    char buf[256];
    int i, j;
    layers _lay = &fdtd_lay;
    int n_goods, tmp_ival;

    if((fp = fopen(fdtd_prj.lay_name, "r")) == NULL) {
        fprintf(stderr, "Cannot open file %s\n", fdtd_prj.lay_name);
        exit(1);
    }

    fgets(buf, 256, fp);//Описание слоев
    fgets(buf, 256, fp);//<Количество слоев>
    fgets(buf, 256, fp);//_lay->n_lay

    sscanf(buf, "%d", &(_lay->n_lay));

    PARMS_malloc(_lay->data_lay, (_lay->n_lay+1), _layer);

    //окамление области
    _lay->data_lay[0].lay_eps = 1.0;
    _lay->data_lay[0].lay_mu = 1.0;
    _lay->data_lay[0].lay_sig = 1.0;

    for(i=1;i<=_lay->n_lay;i++) {//цикл по слоям
        fgets(buf, 256, fp);//<Номер, название слоя>
        fgets(buf, 256, fp);//0   Воздух

        sscanf(buf, "%d", &(_lay->data_lay[i].lay_num));

        fgets(buf, 256, fp);//<количество веществ, газ(0)/не газ(1), давление в слое(газ-атм.)/плотн.(не газ-г/см3), темп.[Цельсий], мод. провод. №, расчет полей(0-нет,1-да)>
        fgets(buf, 256, fp);//1     0     0.00132     20     0     1

        sscanf(buf, "%d", &n_goods);

        fgets(buf, 256, fp);//<ЭМ свойства слоя(0-заданы здесь/1-из веществ): диэл.прон.(0-зав. от частоты), магн. прон., проводимость>
        fgets(buf, 256, fp);//0     1     1     0

        sscanf(buf, "%d%*[ ]%lf%*[ ]%lf%*[ ]%lf", &(tmp_ival),
                                                  &(_lay->data_lay[i].lay_eps),
                                                  &(_lay->data_lay[i].lay_mu),
                                                  &(_lay->data_lay[i].lay_sig));

        for(j=0;j<n_goods;j++) {
            fgets(buf, 256, fp);//<Вещество: номер, доля, название>
            fgets(buf, 256, fp);//0   1   Азот
        }
    }

    fclose(fp);
}

/*!
* Считывание привязки материалов к ячейкам.

*/
void readGi()
{
    FILE *fp;
    char buf[256];
    int i, j, k, n;
    int nind;//количество ячеек 
    int ind;//номер слоя
    int ind_lay;//номер слоя в списке
    grid _grd = &fdtd_grd;

    for (i = 0; i <= _grd->nx+1; ++i) {
        for (j = 0; j <= _grd->ny+1; ++j) {
            for (k = 0; k <= _grd->nz+1; ++k) {
                fdtd_gi[i][j][k] = 1;//вся область заполнена воздухом
            }
        }
    }

    if (!strcmp(fdtd_prj.cel_name, "NONE"))
        return;

    if ((fp = fopen(fdtd_prj.cel_name, "r")) == NULL) {
        fprintf(stderr, "Cannot open file %s\n", fdtd_prj.cel_name);
        exit(1);
    }

    nind = 0;

    for (i = 1; i <= _grd->nx; ++i) {
        for (j = 1; j <= _grd->ny; ++j) {
            for (k = 1; k <= _grd->nz; ++k) {
                if (nind == 0) {
                    do {
                        fgets(buf, 256, fp);
                        sscanf(buf, "%d", &ind);//Номер слоя

                        for (n = 1;n <= fdtd_lay.n_lay; ++n) {
                            if (ind == fdtd_lay.data_lay[n].lay_num)
                                break;
                        }

                        ind_lay = n;//Номер слоя в списке

                        fgets(buf, 256, fp);
                        sscanf(buf, "%d", &nind);//Количество ячеек на этом слое
                    } while (nind <= 0);
                }

                fdtd_gi[i][j][k] = ind_lay;

                nind--;
            }
        }
    }

    fclose(fp);
}

/*!
* Считывание данных о выводе результатов расчета.

*/
void readVid()
{
    FILE *fp;
    char buf[256];
    char tmp[64];
    int n, i, j;
    int data_size;
    project _prj = &fdtd_prj;
    grid _grd = &fdtd_grd;

    for (n = 0; n < _prj->n_output; ++n) {//цикл по файлам выдачи
        if ((fp = fopen(_prj->vid_names[n].data, "r")) == NULL) {
            fprintf(stderr, "Cannot open file %s\n", _prj->vid_names[n].data);
            exit(1);
        }

        do {
            fgets(buf, 256, fp);
            sscanf(buf, "%s", tmp);//ищем "DATA"
        } while (strcmp(tmp, "DATA"));

        fgets(buf, 256, fp);
        fgets(buf, 256, fp);
        sscanf(buf, "%d", &(_prj->out_data[n].out_narr));

        PARMS_malloc(_prj->out_data[n].out_arr, _prj->out_data[n].out_narr, _output_arr);

        fgets(buf, 256, fp);
        fgets(buf, 256, fp);
        sscanf(buf, "%d", &(_prj->out_data[n].out_nt));

        PARMS_malloc(_prj->out_data[n].out_it, _prj->out_data[n].out_nt, int);
        _prj->out_data[n].cur_out_it = 0;
        _prj->out_data[n].max_out_num_data = 0;//начальный максимальный размер

        fgets(buf, 256, fp);
        for(i=0;i<_prj->out_data[n].out_nt;i++)
        {
            fscanf(fp, "%d ", &(_prj->out_data[n].out_it[i]));
        }
        fscanf(fp, "\n");

        for(i=0;i<_prj->out_data[n].out_narr;i++)
        {
            fgets(buf, 256, fp);
            fgets(buf, 256, fp);
            sscanf(buf, "%d", &(_prj->out_data[n].out_arr[i].out_index));

            fgets(buf, 256, fp);
            fgets(buf, 256, fp);
            sscanf(buf, "%d", &(_prj->out_data[n].out_arr[i].out_sgs));

            fgets(buf, 256, fp);
            fgets(buf, 256, fp);
            fgets(buf, 256, fp);
            fgets(buf, 256, fp);

            fgets(buf, 256, fp);
            fgets(buf, 256, fp);
            sscanf(buf, "%d", &(_prj->out_data[n].out_arr[i].out_nx));

            PARMS_malloc(_prj->out_data[n].out_arr[i].out_ix, _prj->out_data[n].out_arr[i].out_nx, int);

            fgets(buf, 256, fp);

            for(j=0;j<_prj->out_data[n].out_arr[i].out_nx;j++)
            {
                fscanf(fp, "%d", &(_prj->out_data[n].out_arr[i].out_ix[j]));
            }

            fscanf(fp, "\n");

            fgets(buf, 256, fp);
            fgets(buf, 256, fp);
            sscanf(buf, "%d", &(_prj->out_data[n].out_arr[i].out_ny));

            PARMS_malloc(_prj->out_data[n].out_arr[i].out_iy, _prj->out_data[n].out_arr[i].out_ny, int);

            fgets(buf, 256, fp);

            for(j=0;j<_prj->out_data[n].out_arr[i].out_ny;j++)
            {
                fscanf(fp, "%d", &(_prj->out_data[n].out_arr[i].out_iy[j]));
            }

            fscanf(fp, "\n");

            fgets(buf, 256, fp);
            fgets(buf, 256, fp);
            sscanf(buf, "%d", &(_prj->out_data[n].out_arr[i].out_nz));

            PARMS_malloc(_prj->out_data[n].out_arr[i].out_iz, _prj->out_data[n].out_arr[i].out_nz, int);

            fgets(buf, 256, fp);

            for(j=0;j<_prj->out_data[n].out_arr[i].out_nz;j++)
            {
                fscanf(fp, "%d", &(_prj->out_data[n].out_arr[i].out_iz[j]));
            }

            fscanf(fp, "\n");

            data_size = _prj->out_data[n].out_arr[i].out_nx *
                        _prj->out_data[n].out_arr[i].out_ny *
                        _prj->out_data[n].out_arr[i].out_nz;

            _prj->out_data[n].max_out_num_data = (data_size > _prj->out_data[n].max_out_num_data) ? data_size : _prj->out_data[n].max_out_num_data;
        }

        if(n==0)//сохранение
        {
            fgets(buf, 256, fp);
            fgets(buf, 256, fp);
            sscanf(buf, "%d", &(_prj->nt_save));

            _prj->cur_it_save = 0;

            if(_prj->nt_save)
            {
                PARMS_malloc(_prj->it_save, _prj->nt_save, int);

                fgets(buf, 256, fp);
                for(i=0;i<_prj->nt_save;i++)
                {
                    fscanf(fp, "%d", &(_prj->it_save[i]));
                }
                fscanf(fp, "\n");
            }
        }
        fclose(fp);
    }
}
