
/*! \file 
* \brief ‘ункции вычисл¤ющие усрдненную сеточную проводимость, диэлектрическую,
         магнитную проницаемость, градиент фронта.
* \author  рюков ј.ј. anton.krv@gmail.com
*/

#include <math.h>
#include "../../include/globals.h"
#include "../../include/constants.h"

/*!
* ”средненна¤ сеточна¤ проводимость.
* \param i - указатель индекса ¤чейки по x-направлению.
* \param j - указатель индекса ¤чейки по y-направлению.
* \param k - указатель индекса ¤чейки по z-направлению.
* \param ind - направление усреднени¤ (1-x, 2-y, 3-z).
* \return усредненна¤ проводимость.
*/
double fun_sig(int i, int j, int k, fieldType idf)
{
    double sig_tmp = 0.0;

    int i0, i1;
    int j0, j1;
    int k0, k1;
    double dx_tmp, dy_tmp, dz_tmp;
    grid _grd = &fdtd_grd;

    switch(idf) {
    case EX_ID :
        if(j==0) {
            j0 = 1;
            j1 = 1;
            dy_tmp = _grd->dyi05[1];
        } else if(j==_grd->ny) {
            j0 = 0;
            j1 = 0;
            dy_tmp = _grd->dyi05[_grd->ny];
        } else {
            j0 = 0;
            j1 = 1;
            dy_tmp = _grd->dyi[j];
        }

        if(k==0) {
            k0 = 1;
            k1 = 1;
            dz_tmp = _grd->dzi05[1];
        } else if(k==_grd->nz) {
            k0 = 0;
            k1 = 0;
            dz_tmp = _grd->dzi05[_grd->nz];
        } else {
            k0 = 0;
            k1 = 1;
            dz_tmp = _grd->dzi[k];
        }

        sig_tmp = 0.25*_grd->dyi05[j+j0]*_grd->dzi05[k+k0]*fdtd_Sig[i][j+j0][k+k0] +
                  0.25*_grd->dyi05[j+j1]*_grd->dzi05[k+k0]*fdtd_Sig[i][j+j1][k+k0] +
                  0.25*_grd->dyi05[j+j0]*_grd->dzi05[k+k1]*fdtd_Sig[i][j+j0][k+k1] +
                  0.25*_grd->dyi05[j+j1]*_grd->dzi05[k+k1]*fdtd_Sig[i][j+j1][k+k1];
        sig_tmp /= (dy_tmp*dz_tmp);

        return sig_tmp;

    case EY_ID :
        if(i==0) {
            i0 = 1;
            i1 = 1;
            dx_tmp = _grd->dxi05[1];
        } else if(i==_grd->nx) {
            i0 = 0;
            i1 = 0;
            dx_tmp = _grd->dxi05[_grd->nx];
        } else {
            i0 = 0;
            i1 = 1;
            dx_tmp = _grd->dxi[i];
        }

        if(k==0) {
            k0 = 1;
            k1 = 1;
            dz_tmp = _grd->dzi05[1];
        } else if(k==_grd->nz) {
            k0 = 0;
            k1 = 0;
            dz_tmp = _grd->dzi05[_grd->nz];
        } else {
            k0 = 0;
            k1 = 1;
            dz_tmp = _grd->dzi[k];
        }

        sig_tmp = 0.25*_grd->dxi05[i+i0]*_grd->dzi05[k+k0]*fdtd_Sig[i+i0][j][k+k0] +
                  0.25*_grd->dxi05[i+i1]*_grd->dzi05[k+k0]*fdtd_Sig[i+i1][j][k+k0] +
                  0.25*_grd->dxi05[i+i0]*_grd->dzi05[k+k1]*fdtd_Sig[i+i0][j][k+k1] +
                  0.25*_grd->dxi05[i+i1]*_grd->dzi05[k+k1]*fdtd_Sig[i+i1][j][k+k1];
        sig_tmp /= (dx_tmp*dz_tmp);

        return sig_tmp;

    case EZ_ID :
        if(i==0) {
            i0 = 1;
            i1 = 1;
            dx_tmp = _grd->dxi05[1];
        } else if(i==_grd->nx) {
            i0 = 0;
            i1 = 0;
            dx_tmp = _grd->dxi05[_grd->nx];
        } else {
            i0 = 0;
            i1 = 1;
            dx_tmp = _grd->dxi[i];
        }

        if(j==0) {
            j0 = 1;
            j1 = 1;
            dy_tmp = _grd->dyi05[1];
        } else if(j==_grd->ny) {
            j0 = 0;
            j1 = 0;
            dy_tmp = _grd->dyi05[_grd->ny];
        } else {
            j0 = 0;
            j1 = 1;
            dy_tmp = _grd->dyi[j];
        }

        sig_tmp = 0.25*_grd->dxi05[i+i0]*_grd->dyi05[j+j0]*fdtd_Sig[i+i0][j+j0][k] +
                  0.25*_grd->dxi05[i+i1]*_grd->dyi05[j+j0]*fdtd_Sig[i+i1][j+j0][k] +
                  0.25*_grd->dxi05[i+i0]*_grd->dyi05[j+j1]*fdtd_Sig[i+i0][j+j1][k] +
                  0.25*_grd->dxi05[i+i1]*_grd->dyi05[j+j1]*fdtd_Sig[i+i1][j+j1][k];
        sig_tmp /= (dx_tmp*dy_tmp);

        return sig_tmp;

    default :
        sig_tmp = 0.0;
    }

    return sig_tmp;
}

/*!
* ”средненна¤ диэлектрическа¤ проницаемость.
* \param i - указатель индекса ¤чейки по x-направлению.
* \param j - указатель индекса ¤чейки по y-направлению.
* \param k - указатель индекса ¤чейки по z-направлению.
* \param ind - направление усреднени¤ (1-x, 2-y, 3-z).
* \return усредненна¤ диэлектрическа¤ проницаемость.
*/
double fun_eps(int i, int j, int k, fieldType idf)
{
    double eps_tmp = 1.0;

    int i0, i1;
    int j0, j1;
    int k0, k1;
    double dx_tmp, dy_tmp, dz_tmp;
    grid _grd = &fdtd_grd;
    layers _lay = &fdtd_lay;
    double x, y, z;

    switch(idf) {
    case EX_ID :
        if(j==0) {
            j0 = 1;
            j1 = 1;
            dy_tmp = _grd->dyi05[1];
        } else if(j==_grd->ny) {
            j0 = 0;
            j1 = 0;
            dy_tmp = _grd->dyi05[_grd->ny];
        } else {
            j0 = 0;
            j1 = 1;
            dy_tmp = _grd->dyi[j];
        }

        if(k==0) {
            k0 = 1;
            k1 = 1;
            dz_tmp = _grd->dzi05[1];
        } else if(k==_grd->nz) {
            k0 = 0;
            k1 = 0;
            dz_tmp = _grd->dzi05[_grd->nz];
        } else {
            k0 = 0;
            k1 = 1;
            dz_tmp = _grd->dzi[k];
        }

        eps_tmp = 0.25*_grd->dyi05[j+j0]*_grd->dzi05[k+k0]*_lay->data_lay[fdtd_gi[i][j+j0][k+k0]].lay_eps +
                  0.25*_grd->dyi05[j+j1]*_grd->dzi05[k+k0]*_lay->data_lay[fdtd_gi[i][j+j1][k+k0]].lay_eps +
                  0.25*_grd->dyi05[j+j0]*_grd->dzi05[k+k1]*_lay->data_lay[fdtd_gi[i][j+j0][k+k1]].lay_eps +
                  0.25*_grd->dyi05[j+j1]*_grd->dzi05[k+k1]*_lay->data_lay[fdtd_gi[i][j+j1][k+k1]].lay_eps;
        eps_tmp /= (dy_tmp*dz_tmp);

        x = _grd->xi05[i];
        y = _grd->yi[j];
        z = _grd->zi[k];

        break;

    case EY_ID :
        if(i==0) {
            i0 = 1;
            i1 = 1;
            dx_tmp = _grd->dxi05[1];
        } else if(i==_grd->nx) {
            i0 = 0;
            i1 = 0;
            dx_tmp = _grd->dxi05[_grd->nx];
        } else {
            i0 = 0;
            i1 = 1;
            dx_tmp = _grd->dxi[i];
        }

        if(k==0) {
            k0 = 1;
            k1 = 1;
            dz_tmp = _grd->dzi05[1];
        } else if(k==_grd->nz) {
            k0 = 0;
            k1 = 0;
            dz_tmp = _grd->dzi05[_grd->nz];
        } else {
            k0 = 0;
            k1 = 1;
            dz_tmp = _grd->dzi[k];
        }

        eps_tmp = 0.25*_grd->dxi05[i+i0]*_grd->dzi05[k+k0]*_lay->data_lay[fdtd_gi[i+i0][j][k+k0]].lay_eps +
                  0.25*_grd->dxi05[i+i1]*_grd->dzi05[k+k0]*_lay->data_lay[fdtd_gi[i+i1][j][k+k0]].lay_eps +
                  0.25*_grd->dxi05[i+i0]*_grd->dzi05[k+k1]*_lay->data_lay[fdtd_gi[i+i0][j][k+k1]].lay_eps +
                  0.25*_grd->dxi05[i+i1]*_grd->dzi05[k+k1]*_lay->data_lay[fdtd_gi[i+i1][j][k+k1]].lay_eps;
        eps_tmp /= (dx_tmp*dz_tmp);

        x = _grd->xi[i];
        y = _grd->yi05[j];
        z = _grd->zi[k];

        break;

    case EZ_ID :
        if(i==0) {
            i0 = 1;
            i1 = 1;
            dx_tmp = _grd->dxi05[1];
        } else if(i==_grd->nx) {
            i0 = 0;
            i1 = 0;
            dx_tmp = _grd->dxi05[_grd->nx];
        } else {
            i0 = 0;
            i1 = 1;
            dx_tmp = _grd->dxi[i];
        }

        if(j==0) {
            j0 = 1;
            j1 = 1;
            dy_tmp = _grd->dyi05[1];
        } else if(j==_grd->ny) {
            j0 = 0;
            j1 = 0;
            dy_tmp = _grd->dyi05[_grd->ny];
        } else {
            j0 = 0;
            j1 = 1;
            dy_tmp = _grd->dyi[j];
        }

        eps_tmp = 0.25*_grd->dxi05[i+i0]*_grd->dyi05[j+j0]*_lay->data_lay[fdtd_gi[i+i0][j+j0][k]].lay_eps +
                  0.25*_grd->dxi05[i+i1]*_grd->dyi05[j+j0]*_lay->data_lay[fdtd_gi[i+i1][j+j0][k]].lay_eps +
                  0.25*_grd->dxi05[i+i0]*_grd->dyi05[j+j1]*_lay->data_lay[fdtd_gi[i+i0][j+j1][k]].lay_eps +
                  0.25*_grd->dxi05[i+i1]*_grd->dyi05[j+j1]*_lay->data_lay[fdtd_gi[i+i1][j+j1][k]].lay_eps;
        eps_tmp /= (dx_tmp*dy_tmp);

        x = _grd->xi[i];
        y = _grd->yi[j];
        z = _grd->zi05[k];

        break;

    default :
        eps_tmp = 1.0;
    }

    return eps_tmp;
}

/*!
* ”средненна¤ магнитна¤ проницаемость.
* \param i - указатель индекса ¤чейки по x-направлению.
* \param j - указатель индекса ¤чейки по y-направлению.
* \param k - указатель индекса ¤чейки по z-направлению.
* \param ind - направление усреднени¤ (1-x, 2-y, 3-z).
* \return усредненна¤ диэлектрическа¤ проницаемость.
*/
double fun_mu(int i, int j, int k, fieldType idf)
{
    double mu_tmp = 1.0;

    int i0, i1;
    int j0, j1;
    int k0, k1;
    double dx_tmp, dy_tmp, dz_tmp;
    grid _grd = &fdtd_grd;
    layers _lay = &fdtd_lay;

    switch(idf) {
    case HX_ID :
        if(i==0) {
            i0 = 1;
            i1 = 1;
            dx_tmp = _grd->dxi05[1];
        } else if(i==_grd->nx) {
            i0 = 0;
            i1 = 0;
            dx_tmp = _grd->dxi05[_grd->nx];
        } else {
            i0 = 0;
            i1 = 1;
            dx_tmp = _grd->dxi[i];
        }

        mu_tmp = 0.5*_grd->dxi05[i+i0]*_lay->data_lay[fdtd_gi[i+i0][j][k]].lay_mu +
                 0.5*_grd->dxi05[i+i1]*_lay->data_lay[fdtd_gi[i+i1][j][k]].lay_mu;
        mu_tmp /= dx_tmp;

        return mu_tmp;

    case HY_ID :
        if(j==0) {
            j0 = 1;
            j1 = 1;
            dy_tmp = _grd->dyi05[1];
        } else if(j==_grd->ny) {
            j0 = 0;
            j1 = 0;
            dy_tmp = _grd->dyi05[_grd->ny];
        } else {
            j0 = 0;
            j1 = 1;
            dy_tmp = _grd->dyi[j];
        }

        mu_tmp = 0.5*_grd->dyi05[j+j0]*_lay->data_lay[fdtd_gi[i][j+j0][k]].lay_mu +
                 0.5*_grd->dyi05[j+j1]*_lay->data_lay[fdtd_gi[i][j+j1][k]].lay_mu;
        mu_tmp /= dy_tmp;

        return mu_tmp;

    case HZ_ID :
        if(k==0) {
            k0 = 1;
            k1 = 1;
            dz_tmp = _grd->dzi05[1];
        } else if(k==_grd->nz) {
            k0 = 0;
            k1 = 0;
            dz_tmp = _grd->dzi05[_grd->nz];
        } else {
            k0 = 0;
            k1 = 1;
            dz_tmp = _grd->dzi[k];
        }

        mu_tmp = 0.5*_grd->dzi05[k+k0]*_lay->data_lay[fdtd_gi[i][j][k+k0]].lay_mu +
                 0.5*_grd->dzi05[k+k1]*_lay->data_lay[fdtd_gi[i][j][k+k1]].lay_mu;
        mu_tmp /= dz_tmp;

        return mu_tmp;

    default :
        mu_tmp = 1.0;
    }

    return mu_tmp;
}
