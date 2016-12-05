
#include <iostream>
#include <string>

#include <pugixml.hpp>
#include <H5Cpp.h>

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

double nfff_time1;
double nfff_time2;

int nfff_nx;
int nfff_ny;
int nfff_nz;

double xMin, xMax;
double yMin, yMax;
double zMin, zMax;

double *grid_x;
double *grid_dx;
double *grid_y;
double *grid_dy;
double *grid_z;
double *grid_dz;

double *nfff_ex_y_min;/*!< массив z- компоненты поля Ex. */
double *nfff_ex_y_max;/*!< массив z- компоненты поля Ex. */
double *nfff_ex_z_min;/*!< массив y- компоненты поля Ex. */
double *nfff_ex_z_max;/*!< массив y- компоненты поля Ex. */

double *nfff_ey_x_min;/*!< массив z- компоненты поля Ey. */
double *nfff_ey_x_max;/*!< массив z- компоненты поля Ey. */
double *nfff_ey_z_min;/*!< массив x- компоненты поля Ey. */
double *nfff_ey_z_max;/*!< массив x- компоненты поля Ey. */

double *nfff_ez_x_min;/*!< массив y- компоненты поля Ez. */
double *nfff_ez_x_max;/*!< массив y- компоненты поля Ez. */
double *nfff_ez_y_min;/*!< массив x- компоненты поля Ez. */
double *nfff_ez_y_max;/*!< массив x- компоненты поля Ez. */

double *nfff_hx_y_min;/*!< массив z- компоненты поля Hx. */
double *nfff_hx_y_max;/*!< массив z- компоненты поля Hx. */
double *nfff_hx_z_min;/*!< массив y- компоненты поля Hx. */
double *nfff_hx_z_max;/*!< массив y- компоненты поля Hx. */

double *nfff_hy_x_min;/*!< массив z- компоненты поля Hy. */
double *nfff_hy_x_max;/*!< массив z- компоненты поля Hy. */
double *nfff_hy_z_min;/*!< массив x- компоненты поля Hy. */
double *nfff_hy_z_max;/*!< массив x- компоненты поля Hy. */

double *nfff_hz_x_min;/*!< массив y- компоненты поля Hz. */
double *nfff_hz_x_max;/*!< массив y- компоненты поля Hz. */
double *nfff_hz_y_min;/*!< массив x- компоненты поля Hz. */
double *nfff_hz_y_max;/*!< массив x- компоненты поля Hz. */

extern void allocate_dft_data();
extern void make_dft_step(const double t2, const double t1);
extern void norm_dft_data();

static void _allocate_nfff_data()
{
    //Xmin
    nfff_ey_x_min = new double[nfff_ny*nfff_nz];
    nfff_ez_x_min = new double[nfff_ny*nfff_nz];
    nfff_hy_x_min = new double[nfff_ny*nfff_nz];
    nfff_hz_x_min = new double[nfff_ny*nfff_nz];

    //Xmax
    nfff_ey_x_max = new double[nfff_ny*nfff_nz];
    nfff_ez_x_max = new double[nfff_ny*nfff_nz];
    nfff_hy_x_max = new double[nfff_ny*nfff_nz];
    nfff_hz_x_max = new double[nfff_ny*nfff_nz];

    //Ymin
    nfff_ex_y_min = new double[nfff_nx*nfff_nz];
    nfff_ez_y_min = new double[nfff_nx*nfff_nz];
    nfff_hx_y_min = new double[nfff_nx*nfff_nz];
    nfff_hz_y_min = new double[nfff_nx*nfff_nz];

    //Ymax
    nfff_ex_y_max = new double[nfff_nx*nfff_nz];
    nfff_ez_y_max = new double[nfff_nx*nfff_nz];
    nfff_hx_y_max = new double[nfff_nx*nfff_nz];
    nfff_hz_y_max = new double[nfff_nx*nfff_nz];

    //Zmin
    nfff_ex_z_min = new double[nfff_nx*nfff_ny];
    nfff_ey_z_min = new double[nfff_nx*nfff_ny];
    nfff_hx_z_min = new double[nfff_nx*nfff_ny];
    nfff_hy_z_min = new double[nfff_nx*nfff_ny];

    //Zmax
    nfff_ex_z_max = new double[nfff_nx*nfff_ny];
    nfff_ey_z_max = new double[nfff_nx*nfff_ny];
    nfff_hx_z_max = new double[nfff_nx*nfff_ny];
    nfff_hy_z_max = new double[nfff_nx*nfff_ny];
}

static void _read_axis(double* &x, double* &dx, int n, const std::string &id, const H5File &file)
{
    std::string name1("/axis_");
    name1 += id;

    DataSet dataset1 = file.openDataSet(name1.c_str());

    x = new double[n];

    // Read the data.
    dataset1.read(x, PredType::NATIVE_DOUBLE);

    std::string name2("/steps_");
    name2 += id;

    DataSet dataset2 = file.openDataSet(name2.c_str());

    dx = new double[n];

    // Read the data.
    dataset2.read(dx, PredType::NATIVE_DOUBLE);
}

void read_nfff_grid()
{
    // Open an existing file and dataset.
    H5File file("nfff_grid.h5", H5F_ACC_RDONLY);

    DataSet dim_dataset = file.openDataSet("/dimensions");

    int dim_data[3];

    // Read the data.
    dim_dataset.read(dim_data, PredType::NATIVE_INT);

    nfff_nx = dim_data[0];
    nfff_ny = dim_data[1];
    nfff_nz = dim_data[2];

    DataSet bnd_dataset = file.openDataSet("/boundaries");

    double bnd_data[2][3];

    // Read the data.
    bnd_dataset.read(bnd_data, PredType::NATIVE_DOUBLE);

    xMin = bnd_data[0][0];
    xMax = bnd_data[1][0];
    yMin = bnd_data[0][1];
    yMax = bnd_data[1][1];
    zMin = bnd_data[0][2];
    zMax = bnd_data[1][2];

    _read_axis(grid_x, grid_dx, nfff_nx, "x", file);
    _read_axis(grid_y, grid_dy, nfff_ny, "y", file);
    _read_axis(grid_z, grid_dz, nfff_nz, "z", file);

    //закрыли файл
    file.close();

    _allocate_nfff_data();

    allocate_dft_data();
}

void _load_nfff_data(const char *file_name)
{
    //открыли файл на чтение
    H5File file(file_name, H5F_ACC_RDONLY);

    DataSet dataset;

    //чтение данных
    //Xmin
    dataset= file.openDataSet("/ey_x_min");
    dataset.read(nfff_ey_x_min, PredType::NATIVE_DOUBLE);

    dataset= file.openDataSet("/ez_x_min");
    dataset.read(nfff_ez_x_min, PredType::NATIVE_DOUBLE);

    dataset= file.openDataSet("/hy_x_min");
    dataset.read(nfff_hy_x_min, PredType::NATIVE_DOUBLE);

    dataset= file.openDataSet("/hz_x_min");
    dataset.read(nfff_hz_x_min, PredType::NATIVE_DOUBLE);

    //Xmax
    dataset= file.openDataSet("/ey_x_max");
    dataset.read(nfff_ey_x_max, PredType::NATIVE_DOUBLE);

    dataset= file.openDataSet("/ez_x_max");
    dataset.read(nfff_ez_x_max, PredType::NATIVE_DOUBLE);

    dataset= file.openDataSet("/hy_x_max");
    dataset.read(nfff_hy_x_max, PredType::NATIVE_DOUBLE);

    dataset= file.openDataSet("/hz_x_max");
    dataset.read(nfff_hz_x_max, PredType::NATIVE_DOUBLE);

    //Ymin
    dataset= file.openDataSet("/ex_y_min");
    dataset.read(nfff_ex_y_min, PredType::NATIVE_DOUBLE);

    dataset= file.openDataSet("/ez_y_min");
    dataset.read(nfff_ez_y_min, PredType::NATIVE_DOUBLE);

    dataset= file.openDataSet("/hx_y_min");
    dataset.read(nfff_hx_y_min, PredType::NATIVE_DOUBLE);

    dataset= file.openDataSet("/hz_y_min");
    dataset.read(nfff_hz_y_min, PredType::NATIVE_DOUBLE);

    //Ymax
    dataset= file.openDataSet("/ex_y_max");
    dataset.read(nfff_ex_y_max, PredType::NATIVE_DOUBLE);

    dataset= file.openDataSet("/ez_y_max");
    dataset.read(nfff_ez_y_max, PredType::NATIVE_DOUBLE);

    dataset= file.openDataSet("/hx_y_max");
    dataset.read(nfff_hx_y_max, PredType::NATIVE_DOUBLE);

    dataset= file.openDataSet("/hz_y_max");
    dataset.read(nfff_hz_y_max, PredType::NATIVE_DOUBLE);

    //Zmin
    dataset= file.openDataSet("/ex_z_min");
    dataset.read(nfff_ex_z_min, PredType::NATIVE_DOUBLE);

    dataset= file.openDataSet("/ey_z_min");
    dataset.read(nfff_ey_z_min, PredType::NATIVE_DOUBLE);

    dataset= file.openDataSet("/hx_z_min");
    dataset.read(nfff_hx_z_min, PredType::NATIVE_DOUBLE);

    dataset= file.openDataSet("/hy_z_min");
    dataset.read(nfff_hy_z_min, PredType::NATIVE_DOUBLE);

    //Zmax
    dataset= file.openDataSet("/ex_z_max");
    dataset.read(nfff_ex_z_max, PredType::NATIVE_DOUBLE);

    dataset= file.openDataSet("/ey_z_max");
    dataset.read(nfff_ey_z_max, PredType::NATIVE_DOUBLE);

    dataset= file.openDataSet("/hx_z_max");
    dataset.read(nfff_hx_z_max, PredType::NATIVE_DOUBLE);

    dataset= file.openDataSet("/hy_z_max");
    dataset.read(nfff_hy_z_max, PredType::NATIVE_DOUBLE);

    //закрыли файл
    file.close();
}

void process_nfff_data()
{
    pugi::xml_document doc;

    pugi::xml_parse_result result = doc.load_file("nfff.xml");

    pugi::xml_node nfff = doc.child("NFFF");

    //начальное время расчета
    nfff_time1 = 0.0;
    nfff_time2 = nfff.attribute("time").as_double();

    pugi::xml_node hdf5 = nfff.child("HDF5");

    for (pugi::xml_node_iterator it = hdf5.begin(); it != hdf5.end(); ++it) {
        nfff_time1 = nfff_time2;
        nfff_time2 = it->attribute("timestep").as_double();
        std::string file(it->attribute("file").as_string());

        _load_nfff_data(file.c_str());

        make_dft_step(nfff_time2, nfff_time1);

        std::cout << nfff_time2 << " " << file << std::endl;
    }

    norm_dft_data();
}
