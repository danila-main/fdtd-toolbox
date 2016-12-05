
#include <vtkQuad.h>
#include <vtkPointData.h>
#include <vtkVersion.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkSmartPointer.h>
#include <vtkREctilinearGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLRectilinearGridReader.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <complex>
#include <omp.h>

#define NUM_THETA 60
#define NUM_PHI 60

extern double *nfff_ex_y_min;/*!< массив z- компоненты поля Ex. */
extern double *nfff_ex_y_max;/*!< массив z- компоненты поля Ex. */
extern double *nfff_ex_z_min;/*!< массив y- компоненты поля Ex. */
extern double *nfff_ex_z_max;/*!< массив y- компоненты поля Ex. */

extern double *nfff_ey_x_min;/*!< массив z- компоненты поля Ey. */
extern double *nfff_ey_x_max;/*!< массив z- компоненты поля Ey. */
extern double *nfff_ey_z_min;/*!< массив x- компоненты поля Ey. */
extern double *nfff_ey_z_max;/*!< массив x- компоненты поля Ey. */

extern double *nfff_ez_x_min;/*!< массив y- компоненты поля Ez. */
extern double *nfff_ez_x_max;/*!< массив y- компоненты поля Ez. */
extern double *nfff_ez_y_min;/*!< массив x- компоненты поля Ez. */
extern double *nfff_ez_y_max;/*!< массив x- компоненты поля Ez. */

extern double *nfff_hx_y_min;/*!< массив z- компоненты поля Hx. */
extern double *nfff_hx_y_max;/*!< массив z- компоненты поля Hx. */
extern double *nfff_hx_z_min;/*!< массив y- компоненты поля Hx. */
extern double *nfff_hx_z_max;/*!< массив y- компоненты поля Hx. */

extern double *nfff_hy_x_min;/*!< массив z- компоненты поля Hy. */
extern double *nfff_hy_x_max;/*!< массив z- компоненты поля Hy. */
extern double *nfff_hy_z_min;/*!< массив x- компоненты поля Hy. */
extern double *nfff_hy_z_max;/*!< массив x- компоненты поля Hy. */

extern double *nfff_hz_x_min;/*!< массив y- компоненты поля Hz. */
extern double *nfff_hz_x_max;/*!< массив y- компоненты поля Hz. */
extern double *nfff_hz_y_min;/*!< массив x- компоненты поля Hz. */
extern double *nfff_hz_y_max;/*!< массив x- компоненты поля Hz. */

extern int nfff_nx;
extern int nfff_ny;
extern int nfff_nz;

extern double *grid_x;
extern double *grid_dx;
extern double *grid_y;
extern double *grid_dy;
extern double *grid_z;
extern double *grid_dz;

extern double xMin, xMax;
extern double yMin, yMax;
extern double zMin, zMax;

using namespace std;

complex<double> *dft_ex_y_min;/*!< образ z- компоненты поля Ex. */
complex<double> *dft_ex_y_max;/*!< образ z- компоненты поля Ex. */
complex<double> *dft_ex_z_min;/*!< образ y- компоненты поля Ex. */
complex<double> *dft_ex_z_max;/*!< образ y- компоненты поля Ex. */

complex<double> *dft_ey_x_min;/*!< образ z- компоненты поля Ey. */
complex<double> *dft_ey_x_max;/*!< образ z- компоненты поля Ey. */
complex<double> *dft_ey_z_min;/*!< образ x- компоненты поля Ey. */
complex<double> *dft_ey_z_max;/*!< образ x- компоненты поля Ey. */

complex<double> *dft_ez_x_min;/*!< образ y- компоненты поля Ez. */
complex<double> *dft_ez_x_max;/*!< образ y- компоненты поля Ez. */
complex<double> *dft_ez_y_min;/*!< образ x- компоненты поля Ez. */
complex<double> *dft_ez_y_max;/*!< образ x- компоненты поля Ez. */

complex<double> *dft_hx_y_min;/*!< образ z- компоненты поля Hx. */
complex<double> *dft_hx_y_max;/*!< образ z- компоненты поля Hx. */
complex<double> *dft_hx_z_min;/*!< образ y- компоненты поля Hx. */
complex<double> *dft_hx_z_max;/*!< образ y- компоненты поля Hx. */

complex<double> *dft_hy_x_min;/*!< образ z- компоненты поля Hy. */
complex<double> *dft_hy_x_max;/*!< образ z- компоненты поля Hy. */
complex<double> *dft_hy_z_min;/*!< образ x- компоненты поля Hy. */
complex<double> *dft_hy_z_max;/*!< образ x- компоненты поля Hy. */

complex<double> *dft_hz_x_min;/*!< образ y- компоненты поля Hz. */
complex<double> *dft_hz_x_max;/*!< образ y- компоненты поля Hz. */
complex<double> *dft_hz_y_min;/*!< образ x- компоненты поля Hz. */
complex<double> *dft_hz_y_max;/*!< образ x- компоненты поля Hz. */

complex<double> *dft_JSy_x_min;
complex<double> *dft_JSz_x_min;
complex<double> *dft_JSy_x_max;
complex<double> *dft_JSz_x_max;

complex<double> *dft_JSx_y_min;
complex<double> *dft_JSz_y_min;
complex<double> *dft_JSx_y_max;
complex<double> *dft_JSz_y_max;

complex<double> *dft_JSx_z_min;
complex<double> *dft_JSy_z_min;
complex<double> *dft_JSx_z_max;
complex<double> *dft_JSy_z_max;

complex<double> *dft_MSy_x_min;
complex<double> *dft_MSz_x_min;
complex<double> *dft_MSy_x_max;
complex<double> *dft_MSz_x_max;

complex<double> *dft_MSx_y_min;
complex<double> *dft_MSz_y_min;
complex<double> *dft_MSx_y_max;
complex<double> *dft_MSz_y_max;

complex<double> *dft_MSx_z_min;
complex<double> *dft_MSy_z_min;
complex<double> *dft_MSx_z_max;
complex<double> *dft_MSy_z_max;

int *dft_ex_y_min_nstep;
int *dft_ex_y_max_nstep;
int *dft_ex_z_min_nstep;
int *dft_ex_z_max_nstep;

int *dft_ey_x_min_nstep;
int *dft_ey_x_max_nstep;
int *dft_ey_z_min_nstep;
int *dft_ey_z_max_nstep;

int *dft_ez_x_min_nstep;
int *dft_ez_x_max_nstep;
int *dft_ez_y_min_nstep;
int *dft_ez_y_max_nstep;

int *dft_hx_y_min_nstep;
int *dft_hx_y_max_nstep;
int *dft_hx_z_min_nstep;
int *dft_hx_z_max_nstep;

int *dft_hy_x_min_nstep;
int *dft_hy_x_max_nstep;
int *dft_hy_z_min_nstep;
int *dft_hy_z_max_nstep;

int *dft_hz_x_min_nstep;
int *dft_hz_x_max_nstep;
int *dft_hz_y_min_nstep;
int *dft_hz_y_max_nstep;

double freq = 2.0e+09;//частота в герцах (не круговая)
double total_nfff_power;//сумарная энергия

void allocate_dft_data()
{
    //Xmin
    dft_ey_x_min = new complex<double>[nfff_ny*nfff_nz];
    dft_ez_x_min = new complex<double>[nfff_ny*nfff_nz];
    dft_hy_x_min = new complex<double>[nfff_ny*nfff_nz];
    dft_hz_x_min = new complex<double>[nfff_ny*nfff_nz];

    dft_ey_x_min_nstep = new int[nfff_ny*nfff_nz];
    dft_ez_x_min_nstep = new int[nfff_ny*nfff_nz];
    dft_hy_x_min_nstep = new int[nfff_ny*nfff_nz];
    dft_hz_x_min_nstep = new int[nfff_ny*nfff_nz];

    //Xmax
    dft_ey_x_max = new complex<double>[nfff_ny*nfff_nz];
    dft_ez_x_max = new complex<double>[nfff_ny*nfff_nz];
    dft_hy_x_max = new complex<double>[nfff_ny*nfff_nz];
    dft_hz_x_max = new complex<double>[nfff_ny*nfff_nz];

    dft_ey_x_max_nstep = new int[nfff_ny*nfff_nz];
    dft_ez_x_max_nstep = new int[nfff_ny*nfff_nz];
    dft_hy_x_max_nstep = new int[nfff_ny*nfff_nz];
    dft_hz_x_max_nstep = new int[nfff_ny*nfff_nz];

    //Ymin
    dft_ex_y_min = new complex<double>[nfff_nx*nfff_nz];
    dft_ez_y_min = new complex<double>[nfff_nx*nfff_nz];
    dft_hx_y_min = new complex<double>[nfff_nx*nfff_nz];
    dft_hz_y_min = new complex<double>[nfff_nx*nfff_nz];

    dft_ex_y_min_nstep = new int[nfff_nx*nfff_nz];
    dft_ez_y_min_nstep = new int[nfff_nx*nfff_nz];
    dft_hx_y_min_nstep = new int[nfff_nx*nfff_nz];
    dft_hz_y_min_nstep = new int[nfff_nx*nfff_nz];

    //Ymax
    dft_ex_y_max = new complex<double>[nfff_nx*nfff_nz];
    dft_ez_y_max = new complex<double>[nfff_nx*nfff_nz];
    dft_hx_y_max = new complex<double>[nfff_nx*nfff_nz];
    dft_hz_y_max = new complex<double>[nfff_nx*nfff_nz];

    dft_ex_y_max_nstep = new int[nfff_nx*nfff_nz];
    dft_ez_y_max_nstep = new int[nfff_nx*nfff_nz];
    dft_hx_y_max_nstep = new int[nfff_nx*nfff_nz];
    dft_hz_y_max_nstep = new int[nfff_nx*nfff_nz];

    //Zmin
    dft_ex_z_min = new complex<double>[nfff_nx*nfff_ny];
    dft_ey_z_min = new complex<double>[nfff_nx*nfff_ny];
    dft_hx_z_min = new complex<double>[nfff_nx*nfff_ny];
    dft_hy_z_min = new complex<double>[nfff_nx*nfff_ny];

    dft_ex_z_min_nstep = new int[nfff_nx*nfff_ny];
    dft_ey_z_min_nstep = new int[nfff_nx*nfff_ny];
    dft_hx_z_min_nstep = new int[nfff_nx*nfff_ny];
    dft_hy_z_min_nstep = new int[nfff_nx*nfff_ny];

    //Zmax
    dft_ex_z_max = new complex<double>[nfff_nx*nfff_ny];
    dft_ey_z_max = new complex<double>[nfff_nx*nfff_ny];
    dft_hx_z_max = new complex<double>[nfff_nx*nfff_ny];
    dft_hy_z_max = new complex<double>[nfff_nx*nfff_ny];

    dft_ex_z_max_nstep = new int[nfff_nx*nfff_ny];
    dft_ey_z_max_nstep = new int[nfff_nx*nfff_ny];
    dft_hx_z_max_nstep = new int[nfff_nx*nfff_ny];
    dft_hy_z_max_nstep = new int[nfff_nx*nfff_ny];

    for (int i = 0; i < nfff_ny*nfff_nz; ++i) {
        //Xmin
        dft_ey_x_min[i] = 0.0;
        dft_ez_x_min[i] = 0.0;
        dft_hy_x_min[i] = 0.0;
        dft_hz_x_min[i] = 0.0;

        dft_ey_x_min_nstep[i] = 1;
        dft_ez_x_min_nstep[i] = 1;
        dft_hy_x_min_nstep[i] = 1;
        dft_hz_x_min_nstep[i] = 1;

        //Xmax
        dft_ey_x_max[i] = 0.0;
        dft_ez_x_max[i] = 0.0;
        dft_hy_x_max[i] = 0.0;
        dft_hz_x_max[i] = 0.0;

        dft_ey_x_max_nstep[i] = 1;
        dft_ez_x_max_nstep[i] = 1;
        dft_hy_x_max_nstep[i] = 1;
        dft_hz_x_max_nstep[i] = 1;
    }

    for (int i = 0; i < nfff_nx*nfff_nz; ++i) {
        //Ymin
        dft_ex_y_min[i] = 0.0;
        dft_ez_y_min[i] = 0.0;
        dft_hx_y_min[i] = 0.0;
        dft_hz_y_min[i] = 0.0;

        dft_ex_y_min_nstep[i] = 1;
        dft_ez_y_min_nstep[i] = 1;
        dft_hx_y_min_nstep[i] = 1;
        dft_hz_y_min_nstep[i] = 1;

        //Ymax
        dft_ex_y_max[i] = 0.0;
        dft_ez_y_max[i] = 0.0;
        dft_hx_y_max[i] = 0.0;
        dft_hz_y_max[i] = 0.0;

        dft_ex_y_max_nstep[i] = 1;
        dft_ez_y_max_nstep[i] = 1;
        dft_hx_y_max_nstep[i] = 1;
        dft_hz_y_max_nstep[i] = 1;
    }

    for (int i = 0; i < nfff_nx*nfff_ny; ++i) {
        //Zmin
        dft_ex_z_min[i] = 0.0;
        dft_ey_z_min[i] = 0.0;
        dft_hx_z_min[i] = 0.0;
        dft_hy_z_min[i] = 0.0;

        dft_ex_z_min_nstep[i] = 1;
        dft_ey_z_min_nstep[i] = 1;
        dft_hx_z_min_nstep[i] = 1;
        dft_hy_z_min_nstep[i] = 1;

        //Zmax
        dft_ex_z_max[i] = 0.0;
        dft_ey_z_max[i] = 0.0;
        dft_hx_z_max[i] = 0.0;
        dft_hy_z_max[i] = 0.0;

        dft_ex_z_max_nstep[i] = 1;
        dft_ey_z_max_nstep[i] = 1;
        dft_hx_z_max_nstep[i] = 1;
        dft_hy_z_max_nstep[i] = 1;
    }
}

static void _update_dff_component(complex<double> &dft_val, const double val, const complex<double> &power, int &nstep)
{
    if (nstep > 1 || val != 0.0) {
        dft_val += val*power;
        ++nstep;
    }
}

void make_dft_step(const double t2, const double t1)
{
    complex<double> power = polar(1.0, -2.0*M_PI*freq*t2);

    for (int i = 0; i < nfff_ny*nfff_nz; ++i) {
        //Xmin
        _update_dff_component(dft_ey_x_min[i], nfff_ey_x_min[i], power, dft_ey_x_min_nstep[i]);
        _update_dff_component(dft_ez_x_min[i], nfff_ez_x_min[i], power, dft_ez_x_min_nstep[i]);
        _update_dff_component(dft_hy_x_min[i], nfff_hy_x_min[i], power, dft_hy_x_min_nstep[i]);
        _update_dff_component(dft_hz_x_min[i], nfff_hz_x_min[i], power, dft_hz_x_min_nstep[i]);

        //Xmax
        _update_dff_component(dft_ey_x_max[i], nfff_ey_x_max[i], power, dft_ey_x_max_nstep[i]);
        _update_dff_component(dft_ez_x_max[i], nfff_ez_x_max[i], power, dft_ez_x_max_nstep[i]);
        _update_dff_component(dft_hy_x_max[i], nfff_hy_x_max[i], power, dft_hy_x_max_nstep[i]);
        _update_dff_component(dft_hz_x_max[i], nfff_hz_x_max[i], power, dft_hz_x_max_nstep[i]);
    }

    for (int i = 0; i < nfff_nx*nfff_nz; ++i) {
        //Ymin
        _update_dff_component(dft_ex_y_min[i], nfff_ex_y_min[i], power, dft_ex_y_min_nstep[i]);
        _update_dff_component(dft_ez_y_min[i], nfff_ez_y_min[i], power, dft_ez_y_min_nstep[i]);
        _update_dff_component(dft_hx_y_min[i], nfff_hx_y_min[i], power, dft_hx_y_min_nstep[i]);
        _update_dff_component(dft_hz_y_min[i], nfff_hz_y_min[i], power, dft_hz_y_min_nstep[i]);

        //Ymax
        _update_dff_component(dft_ex_y_max[i], nfff_ex_y_max[i], power, dft_ex_y_max_nstep[i]);
        _update_dff_component(dft_ez_y_max[i], nfff_ez_y_max[i], power, dft_ez_y_max_nstep[i]);
        _update_dff_component(dft_hx_y_max[i], nfff_hx_y_max[i], power, dft_hx_y_max_nstep[i]);
        _update_dff_component(dft_hz_y_max[i], nfff_hz_y_max[i], power, dft_hz_y_max_nstep[i]);
    }

    for (int i = 0; i < nfff_nx*nfff_ny; ++i) {
        //Zmin
        _update_dff_component(dft_ex_z_min[i], nfff_ex_z_min[i], power, dft_ex_z_min_nstep[i]);
        _update_dff_component(dft_ey_z_min[i], nfff_ey_z_min[i], power, dft_ey_z_min_nstep[i]);
        _update_dff_component(dft_hx_z_min[i], nfff_hx_z_min[i], power, dft_hx_z_min_nstep[i]);
        _update_dff_component(dft_hy_z_min[i], nfff_hy_z_min[i], power, dft_hy_z_min_nstep[i]);

        //Zmax
        _update_dff_component(dft_ex_z_max[i], nfff_ex_z_max[i], power, dft_ex_z_max_nstep[i]);
        _update_dff_component(dft_ey_z_max[i], nfff_ey_z_max[i], power, dft_ey_z_max_nstep[i]);
        _update_dff_component(dft_hx_z_max[i], nfff_hx_z_max[i], power, dft_hx_z_max_nstep[i]);
        _update_dff_component(dft_hy_z_max[i], nfff_hy_z_max[i], power, dft_hy_z_max_nstep[i]);
    }
}

void norm_dft_data()
{
    for (int i = 0; i < nfff_ny*nfff_nz; ++i) {
        //Xmin
        dft_ey_x_min[i] /= ((double)dft_ey_x_min_nstep[i]);
        dft_ez_x_min[i] /= ((double)dft_ez_x_min_nstep[i]);
        dft_hy_x_min[i] /= ((double)dft_hy_x_min_nstep[i]);
        dft_hz_x_min[i] /= ((double)dft_hz_x_min_nstep[i]);

        //Xmax
        dft_ey_x_max[i] /= ((double)dft_ey_x_max_nstep[i]);
        dft_ez_x_max[i] /= ((double)dft_ez_x_max_nstep[i]);
        dft_hy_x_max[i] /= ((double)dft_hy_x_max_nstep[i]);
        dft_hz_x_max[i] /= ((double)dft_hz_x_max_nstep[i]);
    }

    for (int i = 0; i < nfff_nx*nfff_nz; ++i) {
        //Ymin
        dft_ex_y_min[i] /= ((double)dft_ex_y_min_nstep[i]);
        dft_ez_y_min[i] /= ((double)dft_ez_y_min_nstep[i]);
        dft_hx_y_min[i] /= ((double)dft_hx_y_min_nstep[i]);
        dft_hz_y_min[i] /= ((double)dft_hz_y_min_nstep[i]);

        //Ymax
        dft_ex_y_max[i] /= ((double)dft_ex_y_max_nstep[i]);
        dft_ez_y_max[i] /= ((double)dft_ez_y_max_nstep[i]);
        dft_hx_y_max[i] /= ((double)dft_hx_y_max_nstep[i]);
        dft_hz_y_max[i] /= ((double)dft_hz_y_max_nstep[i]);
    }

    for (int i = 0; i < nfff_nx*nfff_ny; ++i) {
        //Zmin
        dft_ex_z_min[i] /= ((double)dft_ex_z_min_nstep[i]);
        dft_ey_z_min[i] /= ((double)dft_ey_z_min_nstep[i]);
        dft_hx_z_min[i] /= ((double)dft_hx_z_min_nstep[i]);
        dft_hy_z_min[i] /= ((double)dft_hy_z_min_nstep[i]);

        //Zmax
        dft_ex_z_max[i] /= ((double)dft_ex_z_max_nstep[i]);
        dft_ey_z_max[i] /= ((double)dft_ey_z_max_nstep[i]);
        dft_hx_z_max[i] /= ((double)dft_hx_z_max_nstep[i]);
        dft_hy_z_max[i] /= ((double)dft_hy_z_max_nstep[i]);
    }
}

static double _get_total_power()
{
    //<S> = Re(ExH*])/2
    double val = 0.0;
    complex<double> Ex;
    complex<double> Ey;
    complex<double> Ez;
    complex<double> Hx;
    complex<double> Hy;
    complex<double> Hz;
    complex<double> poyting;

    double powerXmin = 0.0;
    double powerXmax = 0.0;
    double powerYmin = 0.0;
    double powerYmax = 0.0;
    double powerZmin = 0.0;
    double powerZmax = 0.0;

    //Xmin, Xmax
    for (int jk = 0; jk < nfff_ny*nfff_nz; ++jk) {
        int j = jk/nfff_nz;
        int k = jk%nfff_nz;

        //Xmin
        Ey = dft_ey_x_min[jk];
        Ez = dft_ez_x_min[jk];

        Hy = dft_hy_x_min[jk];
        Hz = dft_hz_x_min[jk];

        poyting = Ey*std::conj(Hz) - Ez*std::conj(Hy);

        powerXmin += std::real(poyting)*-0.5*grid_dy[j]*grid_dz[k];

        //Xmax
        Ey = dft_ey_x_max[jk];
        Ez = dft_ez_x_max[jk];

        Hy = dft_hy_x_max[jk];
        Hz = dft_hz_x_max[jk];

        poyting = Ey*std::conj(Hz) - Ez*std::conj(Hy);

        powerXmax += std::real(poyting)*0.5*grid_dy[j]*grid_dz[k];
    }

    //Ymin, Ymax
    for (int ik = 0; ik < nfff_nx*nfff_nz; ++ik) {
        int i = ik/nfff_nz;
        int k = ik%nfff_nz;

        //Ymin
        Ex = dft_ex_y_min[ik];
        Ez = dft_ez_y_min[ik];

        Hx = dft_hx_y_min[ik];
        Hz = dft_hz_y_min[ik];

        poyting = Ez*std::conj(Hx) - Ex*std::conj(Hz);

        powerYmin += std::real(poyting)*-0.5*grid_dx[i]*grid_dz[k];

        //Ymax
        Ex = dft_ex_y_max[ik];
        Ez = dft_ez_y_max[ik];

        Hx = dft_hx_y_max[ik];
        Hz = dft_hz_y_max[ik];

        poyting = Ez*std::conj(Hx) - Ex*std::conj(Hz);

        powerYmax += std::real(poyting)*0.5*grid_dx[i]*grid_dz[k];
    }

    //Zmin, Zmax
    for (int ij = 0; ij < nfff_nx*nfff_ny; ++ij) {
        int i = ij/nfff_ny;
        int j = ij%nfff_ny;

        //Zmin
        Ex = dft_ex_z_min[ij];
        Ey = dft_ey_z_min[ij];

        Hx = dft_hx_z_min[ij];
        Hy = dft_hy_z_min[ij];

        poyting = Ex*std::conj(Hy) - Ey*std::conj(Hx);

        powerZmin += std::real(poyting)*-0.5*grid_dx[i]*grid_dy[j];

        //Zmax
        Ex = dft_ex_z_max[ij];
        Ey = dft_ey_z_max[ij];

        Hx = dft_hx_z_max[ij];
        Hy = dft_hy_z_max[ij];

        poyting = Ex*std::conj(Hy) - Ey*std::conj(Hx);

        powerZmax += std::real(poyting)*0.5*grid_dx[i]*grid_dy[j];
    }

    val = powerXmin +
          powerXmax +
          powerYmin +
          powerYmax +
          powerZmin +
          powerZmax;

    return val;
}

static void _convert_fields_to_curr()
{
    //             |i  j  k |   {ny*Uz - nz*Uy}
    //US = [nxU] = |nx ny nz| = {nz*Ux - nx*Uz}
    //             |Ux Uy Uz|   {nx*Uy - ny*Ux}

    //JS = [nxH]
    //MS = -[nxE]

    //Xmin, Xmax
    dft_JSy_x_min = dft_hz_x_min;
    dft_JSz_x_min = dft_hy_x_min;
    dft_MSy_x_min = dft_ez_x_min;
    dft_MSz_x_min = dft_ey_x_min;

    dft_JSy_x_max = dft_hz_x_max;
    dft_JSz_x_max = dft_hy_x_max;
    dft_MSy_x_max = dft_ez_x_max;
    dft_MSz_x_max = dft_ey_x_max;
    for (int jk = 0; jk < nfff_ny*nfff_nz; ++jk) {
        //Xmin n = (-1, 0, 0)
        //JSy = Hz
        //JSz = -Hy
        //MSy = -Ez
        //MSz = Ey
        dft_JSz_x_min[jk] = -dft_JSz_x_min[jk];
        dft_MSy_x_min[jk] = -dft_MSy_x_min[jk];

        //Xmax n = (1, 0, 0)
        //JSy = -Hz
        //JSz = Hy
        //MSy = Ez
        //MSz = -Ey
        dft_JSy_x_max[jk] = -dft_JSy_x_max[jk];
        dft_MSz_x_max[jk] = -dft_MSz_x_max[jk];
    }

    //Ymin, Ymax
    dft_JSx_y_min = dft_hz_y_min;
    dft_JSz_y_min = dft_hx_y_min;
    dft_MSx_y_min = dft_ez_y_min;
    dft_MSz_y_min = dft_ex_y_min;

    dft_JSx_y_max = dft_hz_y_max;
    dft_JSz_y_max = dft_hx_y_max;
    dft_MSx_y_max = dft_ez_y_max;
    dft_MSz_y_max = dft_ex_y_max;
    for (int ik = 0; ik < nfff_nx*nfff_nz; ++ik) {
        //Ymin n = (0, -1, 0)
        //JSx = -Hz
        //JSz = Hx
        //MSy = Ez
        //MSz = -Ex
        dft_JSx_y_min[ik] = -dft_JSx_y_min[ik];
        dft_MSz_y_min[ik] = -dft_MSz_y_min[ik];

        //Ymax n = (0, 1, 0)
        //JSx = Hz
        //JSz = -Hx
        //MSy = -Ez
        //Msz = Ex
        dft_JSz_y_max[ik] = -dft_JSz_y_max[ik];
        dft_MSx_y_max[ik] = -dft_MSx_y_max[ik];
    }

    //Zmin, Zmax
    dft_JSx_z_min = dft_hy_z_min;
    dft_JSy_z_min = dft_hx_z_min;
    dft_MSx_z_min = dft_ey_z_min;
    dft_MSy_z_min = dft_ex_z_min;

    dft_JSx_z_max = dft_hy_z_max;
    dft_JSy_z_max = dft_hx_z_max;
    dft_MSx_z_max = dft_ey_z_max;
    dft_MSy_z_max = dft_ex_z_max;
    for (int ij = 0; ij < nfff_nx*nfff_ny; ++ij) {
        //Zmin n = (0, 0, -1)
        //JSx = Hy
        //JSy = -Hx
        //MSx = -Ey
        //MSy = Ex
        dft_JSy_z_min[ij] = -dft_JSy_z_min[ij];
        dft_MSx_z_min[ij] = -dft_MSx_z_min[ij];

        //Zmax n = (0, 0, 1)
        //JSx = -Hy
        //JSy = Hx
        //MSx = Ey
        //MSy = -Ex
        dft_JSx_z_max[ij] = -dft_JSx_z_max[ij];
        dft_MSy_z_max[ij] = -dft_MSy_z_max[ij];
    }
}

static complex<double> _calc_exponent(const double x,
                                      const double y,
                                      const double z,
                                      const double sinTeta,
                                      const double cosTeta,
                                      const double sinFi,
                                      const double cosFi)
{
    complex<double> val;

    double k = (2.0*M_PI*freq)/(2.99792458E+10);

    val = polar(1.0, k*(x*sinTeta*cosFi + y*sinTeta*sinFi + z*cosTeta));

    return (val);
}

static complex<double> _calc_dot_prod_teta(const complex<double> &ux,
                                           const complex<double> &uy,
                                           const complex<double> &uz,
                                           const double sinTeta,
                                           const double cosTeta,
                                           const double sinFi,
                                           const double cosFi)
{
    complex<double> val(0.0, 0.0);

    val += ux*cosTeta*cosFi;
    val += uy*cosTeta*sinFi;
    val -= uz*sinTeta;

    return (val);
}

static complex<double> _calc_dot_prod_fi(const complex<double> &ux,
                                         const complex<double> &uy,
                                         const double sinFi,
                                         const double cosFi)
{
    complex<double> val(0.0, 0.0);

    val -= ux*sinFi;
    val += uy*cosFi;

    return (val);
}

static complex<double> _calc_integral_function_teta(const double x,
                                                    const double y,
                                                    const double z,
                                                    const complex<double> &ux,
                                                    const complex<double> &uy,
                                                    const complex<double> &uz,
                                                    const double sinTeta,
                                                    const double cosTeta,
                                                    const double sinFi,
                                                    const double cosFi)
{
    complex<double> val;
    complex<double> dot_prod = _calc_dot_prod_teta(ux, uy, uz, sinTeta, cosTeta, sinFi, cosFi);
    complex<double> exponent = _calc_exponent(x, y, z, sinTeta, cosTeta, sinFi, cosFi);

    val = exponent*dot_prod;

    return (val);
}

static complex<double> _calc_integral_function_fi(const double x,
                                                  const double y,
                                                  const double z,
                                                  const complex<double> &ux,
                                                  const complex<double> &uy,
                                                  const complex<double> &uz,
                                                  const double sinTeta,
                                                  const double cosTeta,
                                                  const double sinFi,
                                                  const double cosFi)
{
    complex<double> val;
    complex<double> dot_prod = _calc_dot_prod_fi(ux, uy, sinFi, cosFi);
    complex<double> exponent = _calc_exponent(x, y, z, sinTeta, cosTeta, sinFi, cosFi);

    val = exponent*dot_prod;

    return (val);
}

static double _calc_radiation_pattern(const double teta, const double fi)
{
    complex<double> Nteta(0.0, 0.0);
    complex<double> Nfi(0.0, 0.0);
    complex<double> Lteta(0.0, 0.0);
    complex<double> Lfi(0.0, 0.0);

    double sinTeta = sin(teta);
    double cosTeta = cos(teta);
    double sinFi = sin(fi);
    double cosFi = cos(fi);

    double x, y, z;
    double site;

    //Xmin, Xmax
    for (int jk = 0; jk < nfff_ny*nfff_nz; ++jk) {
        int j = jk/nfff_nz;
        int k = jk%nfff_nz;

        y = grid_y[j];
        z = grid_z[k];
        site = grid_dy[j]*grid_dz[k];

        //Xmin
        x = xMin;
        //JS
        Nteta += site*_calc_integral_function_teta(x, y, z,
                                                   0.0, dft_JSy_x_min[jk], dft_JSz_x_min[jk],
                                                   sinTeta, cosTeta, sinFi, cosFi);
        Nfi += site*_calc_integral_function_fi(x, y, z,
                                               0.0, dft_JSy_x_min[jk], dft_JSz_x_min[jk],
                                               sinTeta, cosTeta, sinFi, cosFi);

        //MS
        Lteta += site*_calc_integral_function_teta(x, y, z,
                                                   0.0, dft_MSy_x_min[jk], dft_MSz_x_min[jk],
                                                   sinTeta, cosTeta, sinFi, cosFi);
        Lfi += site*_calc_integral_function_fi(x, y, z,
                                               0.0, dft_MSy_x_min[jk], dft_MSz_x_min[jk],
                                               sinTeta, cosTeta, sinFi, cosFi);

        //Xmax
        x = xMax;
        //JS
        Nteta += site*_calc_integral_function_teta(x, y, z,
                                                   0.0, dft_JSy_x_max[jk], dft_JSz_x_max[jk],
                                                   sinTeta, cosTeta, sinFi, cosFi);
        Nfi += site*_calc_integral_function_fi(x, y, z,
                                               0.0, dft_JSy_x_max[jk], dft_JSz_x_max[jk],
                                               sinTeta, cosTeta, sinFi, cosFi);

        //MS
        Lteta += site*_calc_integral_function_teta(x, y, z,
                                                   0.0, dft_MSy_x_max[jk], dft_MSz_x_max[jk],
                                                   sinTeta, cosTeta, sinFi, cosFi);
        Lfi += site*_calc_integral_function_fi(x, y, z,
                                               0.0, dft_MSy_x_max[jk], dft_MSz_x_max[jk],
                                               sinTeta, cosTeta, sinFi, cosFi);
    }

    //Ymin, Ymax
    for (int ik = 0; ik < nfff_nx*nfff_nz; ++ik) {
        int i = ik/nfff_nz;
        int k = ik%nfff_nz;

        x = grid_x[i];
        z = grid_z[k];
        site = grid_dx[i]*grid_dz[k];

        //Ymin
        y = yMin;
        //JS
        Nteta += site*_calc_integral_function_teta(x, y, z,
                                                   dft_JSx_y_min[ik], 0.0, dft_JSz_y_min[ik],
                                                   sinTeta, cosTeta, sinFi, cosFi);
        Nfi += site*_calc_integral_function_fi(x, y, z,
                                               dft_JSx_y_min[ik], 0.0, dft_JSz_y_min[ik],
                                               sinTeta, cosTeta, sinFi, cosFi);

        //MS
        Lteta += site*_calc_integral_function_teta(x, y, z,
                                                   dft_MSx_y_min[ik], 0.0, dft_MSz_y_min[ik],
                                                   sinTeta, cosTeta, sinFi, cosFi);
        Lfi += site*_calc_integral_function_fi(x, y, z,
                                               dft_MSx_y_min[ik], 0.0, dft_MSz_y_min[ik],
                                               sinTeta, cosTeta, sinFi, cosFi);

        //Ymax
        y = yMax;
        //JS
        Nteta += site*_calc_integral_function_teta(x, y, z,
                                                   dft_JSx_y_max[ik], 0.0, dft_JSz_y_max[ik],
                                                   sinTeta, cosTeta, sinFi, cosFi);
        Nfi += site*_calc_integral_function_fi(x, y, z,
                                               dft_JSx_y_max[ik], 0.0, dft_JSz_y_max[ik],
                                               sinTeta, cosTeta, sinFi, cosFi);

        //MS
        Lteta += site*_calc_integral_function_teta(x, y, z,
                                                   dft_MSx_y_max[ik], 0.0, dft_MSz_y_max[ik],
                                                   sinTeta, cosTeta, sinFi, cosFi);
        Lfi += site*_calc_integral_function_fi(x, y, z,
                                               dft_MSx_y_max[ik], 0.0, dft_MSz_y_max[ik],
                                               sinTeta, cosTeta, sinFi, cosFi);
    }

    //Zmin, Zmax
    for (int ij = 0; ij < nfff_nx*nfff_ny; ++ij) {
        int i = ij/nfff_ny;
        int j = ij%nfff_ny;

        x = grid_x[i];
        y = grid_y[j];
        site = grid_dx[i]*grid_dy[j];

        //Zmin
        z = zMin;
        //JS
        Nteta += site*_calc_integral_function_teta(x, y, z,
                                                   dft_JSx_z_min[ij], dft_JSy_z_min[ij], 0.0,
                                                   sinTeta, cosTeta, sinFi, cosFi);
        Nfi += site*_calc_integral_function_fi(x, y, z,
                                               dft_JSx_z_min[ij], dft_JSy_z_min[ij], 0.0,
                                               sinTeta, cosTeta, sinFi, cosFi);

        //MS
        Lteta += site*_calc_integral_function_teta(x, y, z,
                                                   dft_MSx_z_min[ij], dft_MSy_z_min[ij], 0.0,
                                                   sinTeta, cosTeta, sinFi, cosFi);
        Lfi += site*_calc_integral_function_fi(x, y, z,
                                               dft_MSx_z_min[ij], dft_MSy_z_min[ij], 0.0,
                                               sinTeta, cosTeta, sinFi, cosFi);

        //Zmax
        z = zMax;
        //JS
        Nteta += site*_calc_integral_function_teta(x, y, z,
                                                   dft_JSx_z_max[ij], dft_JSy_z_max[ij], 0.0,
                                                   sinTeta, cosTeta, sinFi, cosFi);
        Nfi += site*_calc_integral_function_fi(x, y, z,
                                               dft_JSx_z_max[ij], dft_JSy_z_max[ij], 0.0,
                                               sinTeta, cosTeta, sinFi, cosFi);

        //MS
        Lteta += site*_calc_integral_function_teta(x, y, z,
                                                   dft_MSx_z_max[ij], dft_MSy_z_max[ij], 0.0,
                                                   sinTeta, cosTeta, sinFi, cosFi);
        Lfi += site*_calc_integral_function_fi(x, y, z,
                                               dft_MSx_z_max[ij], dft_MSy_z_max[ij], 0.0,
                                               sinTeta, cosTeta, sinFi, cosFi);
    }

    double val = pow(abs(Lfi + Nteta), 2) + pow(abs(Lteta - Nfi), 2);

    double k = (2.0*M_PI*freq)/(2.99792458E+10);

    val *= (k*k)/(8.0*M_PI);

    val /= total_nfff_power;

    return (val);
}

static double _get_dB(double val)
{
    double n_val = 10.0*log10(val);

    if (n_val <= -500.0) {
        n_val = -500.0;
    }

    return (n_val);
}

static double _get_pattern(double val)
{
    return fabs(val);
}

static void _write_radiation_pattern(const char *file_name)
{
    double theta, phi;
    double x,y, z;
    double r = 1.0;
    double val;

    vtkPoints *points = vtkPoints::New();
    vtkIdType id;
    vtkDoubleArray *data = vtkDoubleArray::New();
    data->SetName("Power");

    double pattern_value[NUM_PHI][NUM_THETA - 1];

    {
    #pragma omp parallel for private(theta, phi)
        for (int ij = 0; ij < (NUM_PHI)*(NUM_THETA - 1); ++ij) {
            int i = ij/(NUM_THETA - 1);
            int j = ij%(NUM_THETA - 1);

            theta = (j + 1)*(M_PI/NUM_THETA);

            phi = i*(2*M_PI/NUM_PHI);

            pattern_value[i][j] = _calc_radiation_pattern(theta, phi);

            std::cout << "point: (" << j <<", " << i << ")" << std::endl;
        }
    }

    //for (int j = 0; j < NUM_THETA - 1; ++j) {
    //    theta = (j + 1)*(M_PI/NUM_THETA);
    //    for (int i = 0; i < NUM_PHI; ++i) {
    //        phi = i*(2*M_PI/NUM_PHI);

    //        pattern_value[i][j] = _calc_radiation_pattern(theta, phi);

    //        std::cout << "point: (" << j <<", " << i << ")" << std::endl;
    //    }
    //}

    for (int j = 0; j < NUM_THETA - 1; ++j) {
        theta = (j + 1)*(M_PI/NUM_THETA);
        for (int i = 0; i < NUM_PHI; ++i) {
            phi = i*(2*M_PI/NUM_PHI);

            val = pattern_value[i][j];
            r = _get_pattern(val);
            val = _get_dB(val);//dB

            x = r*sin(theta)*cos(phi);
            y = r*sin(theta)*sin(phi);
            z = r*cos(theta);

            id = points->InsertNextPoint(x, y, z);

            data->InsertNextValue(val);
        }
    }

    vtkUnstructuredGrid *unstructuredGrid = vtkUnstructuredGrid::New();
    unstructuredGrid->Allocate();

    vtkIdType pts[4];

    for (int j = 1; j < NUM_THETA - 1; ++j) {
        pts[1] = (j-1)*NUM_PHI + (NUM_PHI-1);//     3------2
        pts[0] = (j-1)*NUM_PHI;//                   |(0,j) |
        pts[3] = (j)*NUM_PHI;//                     |      |
        pts[2] = (j)*NUM_PHI + (NUM_PHI-1);//       0------1

        unstructuredGrid->InsertNextCell(VTK_QUAD, 4, pts);

        for (int i = 1; i < NUM_PHI; ++i) {
            pts[1] = (j-1)*NUM_PHI + (i-1);//       3------2
            pts[0] = (j-1)*NUM_PHI + (i);//         |(i,j) |
            pts[3] = (j)*NUM_PHI + (i);//           |      |
            pts[2] = (j)*NUM_PHI + (i-1);//         0------1

            unstructuredGrid->InsertNextCell(VTK_QUAD, 4, pts);
        }
    }

    //полюс +
    val = _calc_radiation_pattern(0.0, 0.0);
    r = _get_pattern(val);
    val = _get_dB(val);//dB

    id = points->InsertNextPoint(0.0, 0.0, r);

    pts[0] = id;
    pts[1] = NUM_PHI - 1;
    pts[2] = 0;

    unstructuredGrid->InsertNextCell(VTK_TRIANGLE, 3, pts);

    for (int i = 1; i < NUM_PHI; ++i) {
        pts[0] = id;
        pts[1] = i - 1;
        pts[2] = i;

        unstructuredGrid->InsertNextCell(VTK_TRIANGLE, 3, pts);
    }

    data->InsertNextValue(val);

    //полюс -
    val = _calc_radiation_pattern(M_PI, 0.0);
    r = _get_pattern(val);
    val = _get_dB(val);//dB

    id = points->InsertNextPoint(0.0, 0.0, -r);

    pts[0] = id;
    pts[1] = (NUM_THETA-2)*NUM_PHI;
    pts[2] = (NUM_THETA-2)*NUM_PHI + NUM_PHI - 1;

    unstructuredGrid->InsertNextCell(VTK_TRIANGLE, 3, pts);

    for (int i = 1; i < NUM_PHI; ++i) {
        pts[0] = id;
        pts[1] = (NUM_THETA-2)*NUM_PHI + i;
        pts[2] = (NUM_THETA-2)*NUM_PHI + i - 1;

        unstructuredGrid->InsertNextCell(VTK_TRIANGLE, 3, pts);
    }

    data->InsertNextValue(val);

    unstructuredGrid->SetPoints(points);
    unstructuredGrid->GetPointData()->AddArray(data);

    // Write file
    vtkXMLUnstructuredGridWriter *writer = vtkXMLUnstructuredGridWriter::New();
    writer->SetFileName(file_name);
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(unstructuredGrid);
#else
    writer->SetInputData(unstructuredGrid);
#endif
    writer->Write();

    points->Delete();
    unstructuredGrid->Delete();
    data->Delete();
    writer->Delete();
}

void get_radiation_pattern()
{
    total_nfff_power = _get_total_power();

    //конвертация
    _convert_fields_to_curr();

    _write_radiation_pattern("pattern.vtu");
}
