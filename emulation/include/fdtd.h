/*! \file 
* \brief Функции расчета полей с идеальным поглащающим слоем на границе.
* \author Крюков А.А. anton.krv@gmail.com
*/

#ifndef _FDTD_HEADER_
#define _FDTD_HEADER_

/*! Толщина PML слоя. */
#define PML_WIDTH 6
/*! Толщина TFSF слоя. */
#define TFSF_WIDTH 7 //TFSF слой должен быть толще PML слоя

extern void allocatePML();
extern void set_PML_profile();

extern void init_plane_wave();
extern void calc_plane_wave(double time);

extern void init_debye_data();
extern void update_debye_data(double ht);

extern void calc_Ex_point_debye(int i, int j, int k);
extern void calc_Ey_point_debye(int i, int j, int k);
extern void calc_Ez_point_debye(int i, int j, int k);

extern void init_drude_data();
extern void update_drude_data(double ht);

extern void calc_Ex_point_drude(int i, int j, int k);
extern void calc_Ey_point_drude(int i, int j, int k);
extern void calc_Ez_point_drude(int i, int j, int k);

typedef enum { X_DIR = 0, Y_DIR, Z_DIR } direction;//направление

typedef double pml_double[PML_WIDTH];
typedef double mur_double[2][2];

extern pml_double k_x_min;/*!< массив коэффициентов маштабирования PML слоя по x. */
extern pml_double k_y_min;/*!< массив коэффициентов маштабирования PML слоя по y. */
extern pml_double k_z_min;/*!< массив коэффициентов маштабирования PML слоя по z. */

extern pml_double k_x_max;/*!< массив коэффициентов маштабирования PML слоя по x. */
extern pml_double k_y_max;/*!< массив коэффициентов маштабирования PML слоя по y. */
extern pml_double k_z_max;/*!< массив коэффициентов маштабирования PML слоя по z. */

extern double alpha_i_e_min(int i, direction dir);
extern double alpha_i_h_min(int i, direction dir);
extern double alpha_i_e_max(int i, direction dir);
extern double alpha_i_h_max(int i, direction dir);

extern double betta_i_e_min(int i, direction dir);
extern double betta_i_h_min(int i, direction dir);
extern double betta_i_e_max(int i, direction dir);
extern double betta_i_h_max(int i, direction dir);

extern pml_double **phu_ex_y_min;/*!< массив y- компоненты свертки для поля Ex. */
extern pml_double **phu_ex_z_min;/*!< массив z- компоненты свертки для поля Ex. */
extern pml_double **phu_ey_x_min;/*!< массив x- компоненты свертки для поля Ey. */
extern pml_double **phu_ey_z_min;/*!< массив z- компоненты свертки для поля Ey. */
extern pml_double **phu_ez_x_min;/*!< массив x- компоненты свертки для поля Ez. */
extern pml_double **phu_ez_y_min;/*!< массив y- компоненты свертки для поля Ez. */

extern pml_double **phu_ex_y_max;/*!< массив y- компоненты свертки для поля Ex. */
extern pml_double **phu_ex_z_max;/*!< массив z- компоненты свертки для поля Ex. */
extern pml_double **phu_ey_x_max;/*!< массив x- компоненты свертки для поля Ey. */
extern pml_double **phu_ey_z_max;/*!< массив z- компоненты свертки для поля Ey. */
extern pml_double **phu_ez_x_max;/*!< массив x- компоненты свертки для поля Ez. */
extern pml_double **phu_ez_y_max;/*!< массив y- компоненты свертки для поля Ez. */

extern pml_double **phu_hx_y_min;/*!< массив y- компоненты свертки для поля Hx. */
extern pml_double **phu_hx_z_min;/*!< массив z- компоненты свертки для поля Hx. */
extern pml_double **phu_hy_x_min;/*!< массив x- компоненты свертки для поля Hy. */
extern pml_double **phu_hy_z_min;/*!< массив z- компоненты свертки для поля Hy. */
extern pml_double **phu_hz_x_min;/*!< массив x- компоненты свертки для поля Hz. */
extern pml_double **phu_hz_y_min;/*!< массив y- компоненты свертки для поля Hz. */

extern pml_double **phu_hx_y_max;/*!< массив y- компоненты свертки для поля Hx. */
extern pml_double **phu_hx_z_max;/*!< массив z- компоненты свертки для поля Hx. */
extern pml_double **phu_hy_x_max;/*!< массив x- компоненты свертки для поля Hy. */
extern pml_double **phu_hy_z_max;/*!< массив z- компоненты свертки для поля Hy. */
extern pml_double **phu_hz_x_max;/*!< массив x- компоненты свертки для поля Hz. */
extern pml_double **phu_hz_y_max;/*!< массив y- компоненты свертки для поля Hz. */

extern void apply_Ex_tfsf(int i);
extern void apply_Ey_tfsf(int j);
extern void apply_Ez_tfsf(int k);

extern void apply_Hx_tfsf(int i);
extern void apply_Hy_tfsf(int j);
extern void apply_Hz_tfsf(int k);

extern void init_nfff_domain();
extern void load_nfff_domain();
extern void write_nfff_domain(int n, double time);
extern void write_nfff_grid();

#endif
