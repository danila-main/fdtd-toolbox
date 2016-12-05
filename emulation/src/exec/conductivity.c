
/*! \file 
* \brief Функции расчета проводимости.
* \author Крюков А.А. anton.krv@gmail.com
*/

#include <math.h>
#include "../../include/globals.h"
#include "../../include/constants.h"

//вычисление тока проводимости в точке 
//(модель 2 - для полупроводников по формулам )
void sub_sigma2(double wt, double me0, double mp0, double s, double e, double qc,
                double *ce_ptr, double *sg_ptr)
{
//wt - температура слоя
//me0, mp0 - подвижность электронов и дырок
//s - скорость звука
//e - модуль электрического поля
//qc - источник вторичных электронов
//ce - концентрация электронов
//sg - проводимость

    double qc1, ce1, ce2, t0dt, t0dt32, me0t, mp0t, x1, x2, y1, y2;
    double me, mp;
    double ht = spars_ht;
    double ce = *ce_ptr;
    double sg = *sg_ptr;

    t0dt = _wt0/wt;
    t0dt32 = t0dt*sqrt(t0dt);

    qc1 = qc;
    ce1 = ce;
    ce2 = ce1 + qc*ht;

    me0t = me0*t0dt32;
    x1 = me0t*e / s;
    x2 = x1*x1;
    me = me0t*sqrt( 2 / ( 1.0+sqrt(1.0+ _pi38*x2) ) ) ;
     
    mp0t = mp0*t0dt32;
    y1 = mp0t*e / s;
    y2 = x1*x1;
    mp = mp0t*sqrt( 2 / ( 1.0+sqrt(1.0+ _pi38*x2) ) ) ;

    sg = ce1 * _qe * ( me + mp );

    ce = ce2;

    *ce_ptr = ce;
    *sg_ptr = sg;
}

//вычисление тока проводимости в точке 
//(модель 3 - для аморфных диэлектриков)
void sub_sigma3(double wt, double wt1, double me0, double nt0, double t, double qc,
                double *ce_ptr, double *ci_ptr, double *sg_ptr)
{
//wt - температура слоя = 1.5*w_bolz*t
//wt1 - температура ловушек = 1.5*w_bolz*t1
//me0 - подвижность электронов
//nt0 - концентрация ловушек
//t - время (c)
//qc - источник вторичных электронов
//ce, ci - концентрации электронов и ионов
//sg - проводимость

    double ce1, ce2;
    double alfa, alfa1, tau, fri;
    double r_trp;
    double d;
    double ht = spars_ht;
    double ce = *ce_ptr;
    double ci = *ci_ptr;
    double sg = *sg_ptr;

    ce1 = ce;//концентрации свободных электронов    

    alfa  = 0.5;//отношение температуры среды к температуре ловушек
    alfa1 = alfa/(1.0 + alfa);
    tau   = 4.0E-11;//время захвата ловушками
    fri   = 3.0E4;//частота освобождения ловушками
    r_trp = 3.0E-8;//радиус ловушки

    d = 2.0*me0*wt/(3.0*_qe);//коэффициент диффузии электрона

    ce2 = ce1 + qc * ht;//интеграл от источника вторичных
     
    sg = _qe*me0*tau*(qc+alfa1*fri*ce1);//проводимость

    ce = ce2;

    *ce_ptr = ce;
    *ci_ptr = ci;
    *sg_ptr = sg;
}

//прилипание к кислороду
static double gam_cfti(double x, double p)
{
//x    - отношение модуля эл поля к давлению[атм]
//p    - давление[атм]
    double absx;
    double val = 0.0;

    absx = fabs(x);

    if (absx == 0.0) {
        val = (0.63E8*p)/sqrt(1.055);
    } else {
        val = (0.63E8*p)/sqrt(1.055+absx) + 1.3E8*exp( (-25.0/(absx + 1.0E-4/p)) );
    }

    return val;
}

//равновесная дрейфовая скорость электронов
static double ue_ter(double x)
{
//x - отношение модуля эл поля к давлению[атм]
    double  ax;
    double val = 0.0;
  
    ax = fabs(x);

    if(ax < 421)
        if(ax < 6.65)
            if(ax < 0.251)
                val = 1.47E6*ax;
            else
                val = 0.737E6*sqrt(ax);
        else
            val = 0.459E6*pow(ax, 0.75);
    else
        val = 2.08E6*sqrt(ax);

    return val;
}

//коэффициент таунсенда
static double al_ter(double x)
{
//x - отношение модуля эл поля к давлению[атм]
    double ax;
    double val = 0.0;
    
    ax = fabs(x);

    if(ax < 330.0)
        if(ax < 97.1)
            if(ax < 1.0E-8)
                val = 0.0;
            else {
                ax = -360.0 / ax;
                val = 5.04E2 * exp(ax);
            }
        else {
            ax =  ax - 66.0;
            val = 1.28E-02*pow(ax, 2.0);
        }
    else
        val = 2.44E2*(sqrt(ax)-14.5);

    return val;
}

//вычисление тока проводимости в точке 
//(модель 4 - стандартная для воздуха)
void sub_sigma4(double p, double rk_i, double alfa_ei, double alfa_ii, double e, double qc,
                double *ce_ptr, double *ci_ptr, double *sg_ptr)
{
//p - давление слоя
//rk_i - подвижность ионов в слое
//alfa_ei - коэффициент рекомбинации для слоя
//alfa_ii - коэффициент рекомбинации для слоя
//e - модуль электрического поля
//qc - источник вторичных электронов
//ce,ci - концентрации электронов и ионов
//sg - проводимость

    double qc1, ce1, ce2, ci1, ci2, gm1, ue1, al1, ry, py, ay, by;

    double ht = spars_ht;
    double ce = *ce_ptr;
    double ci = *ci_ptr;
    double sg = *sg_ptr;

    ce1 = ce;//концентрация вторичных электронов
    ci1 = ci;//концентрация отрицательных ионов 
    qc1 = qc;//источник вторичных
    gm1 = gam_cfti(e/p, p);//прилипание к кислороду
    ue1 = ue_ter(e/p);//модуль скорости вторичных электронов
    al1 = al_ter(e/p);//коэффициент таунсенда

//вычисление концентрации вторичных электронов
    ry = p * (al1*ue1 - gm1) - alfa_ei * (ci1 + ce1);
    py = ry * ht;
    ay = exp(py);
      
    if( py < 1.0E-10 )
        by = ht;
    else
        by = -(1.0 - ay) / ry;

    ce2 = ce1 * ay + qc1 * by;
    ce1=0.5*(ce1+ce2);

//вычисление концентрации отрицательных ионов
    ry = -alfa_ii * (ci1 + ce1);
    py = ry * ht;
    ay = exp(py);
      
    if( py < 1.0E-10 )
        by = ht;
    else
        by = -(1.0 - ay) / ry;

    ci2 = ci1 * ay + (p * gm1 * ce1) * by;

    ce = ce2;
    ci = ci2;

//вычисление проводимости
    if (e == 0.0)
        sg = 0.0;
    else
        sg = _qe*( ce1*ue1/e + rk_i*(ce1+2.0*ci1) );

    *ce_ptr = ce;
    *ci_ptr = ci;
    *sg_ptr = sg;
}

//модуль скорости вторичных электронов в азоте
static double ue_n2(double x)
{
    double ax, x0, x1, x2, x3, x4, x5, x6, val;

    ax = 0.395*fabs(x);//     ! отношение модуля эл поля к давлению[в/(см*торр)]

    x0 = 0.0;
    x1 = 0.003901;
    x2 = 0.01757;
    x3 = 0.03901;
    x4 = 0.09753;
    x5 = 0.7806;
    x6 = 22.08;

    if(ax <= x1)
        val = 46000 * (ax / 0.003901);
    else if(ax <= x2)
        val = 46000 * pow((ax / x1), 0.9);
    else if(ax <= x3)
        val = 171000 * pow((ax / x2), 0.6);
    else if(ax <= x4)
        val = 255000 * pow((ax / x3) ,0.23);
    else if(ax <= x5)
        val = 309000 + 550000 * (ax - x4);
    else if(ax <= x6)
        val = 670000 * pow((ax / x5), 0.73);
    else
        val = 670000 * pow((x6 / x5), 0.73);

    return val;
}

//стандартная для азота
void sub_sigma8(double p, double rk_i, double alfa_ei, double e, double qc, 
                double *ce_ptr, double *sg_ptr)
{
    //p - давление слоя
    //rk_i - подвижность ионов в слое
    //alfa_ei - коэффициент рекомбинации для слоя
    //alfa_ii - коэффициент рекомбинации для слоя
    //e - модуль электрического поля
    //qc - источник вторичных электронов
    //ce,ci - концентрации электронов и ионов
    //sg - проводимость

    double  qc1,ce1,ce2,ue1,al1,ry,py,ay,by;
    double ht = spars_ht;

    ce1 = *ce_ptr;//концентрация вторичных электронов
    qc1 = qc;//источник вторичных

    ue1 = ue_n2(e/p);//модуль скорости вторичных электронов в азоте
    al1 = al_ter(e/p);//коэффициент таунсенда

    //вычисление концентрации вторичных электронов
    ry = p * (al1*ue1) - alfa_ei * ce1;
    py = ry * ht;
    ay = exp(py);
    if( py < 1.0E-10 )
        by = ht;
    else
        by = -(1.0 - ay) / ry;
    
    ce2 = ce1 * ay + qc1 * by;
    ce1 = 0.5*(ce1 + ce2);

    *ce_ptr = ce2;

    //вычисление проводимости
    if (e == 0.0)
        *sg_ptr = 0.0; 
    else
        *sg_ptr = _qe*( ce1*ue1/e + rk_i*ce1 );
}

//(модель 9) для топлива
void sub_sigma9(double lay_evtor, double lay_ro, double qc, double *sg_ptr)
{
    double ev_mev, ev_rad, sg_si;

    ev_mev = lay_evtor * 1.0E-6 * qc/lay_ro;//энерговыделение [мэв/(г*с)]
    ev_rad = ev_mev * _w_1ev*1.0E4;//энерговыделение [рад/с]

    sg_si = 2.0E-16 * ev_rad;//!проводимость [1/(ом*м)]
    *sg_ptr = 9.0E9 * sg_si;//!проводимость [1/c]
}


//равновесная дрейфовая скорость электронов в воздухе
static double u_air(double x)
{
    //x - отношение модуля эл поля к давлению[атм]

    double ax;
    double val = 0.0;

    ax = fabs(x);

    if (ax < 421.0) {
        if (ax < 6.65) {
            if (ax < 0.251) {
                val = 1.47E06*ax;
            } else {
                val = 0.737E06*sqrt(ax);
            }
        } else {
            val = 0.459E06*pow(ax, 0.75);
        }
    } else {
        val = 2.08E06*sqrt(ax);
    }

    return val;
}

 //коэффициент таунсенда в воздухе с плотностью соответствующей р=1
static double alfa_air(double x)
{
    //x - отношение модуля эл. поля к давлению[атм]
    double ax;
    double val = 0.0;

    ax = fabs(x);

    if (ax < 330.0) {
        if (ax < 97.1) {
            if (ax < 1.0E-08) {
                val = 0.0;
            } else {
                val = 5.04E02*exp(-360.0/ax);
            }
        } else {
            val = 1.28E-02*(ax - 66.0)*(ax - 66.0);
        }
    } else {
        val = 2.44E02*(sqrt(ax) - 14.5);
    }

    return val;
}

 //прилипание к кислороду в воздухе с плотностью соответствующей р=1
static double gamma_air(double x, double p)
{
    //x - отношение модуля эл поля к давлению[атм]
    //p - давление[атм]

    double ax;
    double val = 0.0;

    ax = fabs(x);

    if (ax == 0.0) {
        val = 0.63E08*p/sqrt(1.055);
    } else {
        val = 0.63E08*p/sqrt(1.055 + ax) + 1.3E08*exp(-25.0/(ax + 1.0E-04/p));
    }

    return val;
}

 //отлипание от кислорода в воздухе с плотностью соответствующей р=1

double gammd_air(double x, double p)
{
    //x - отношение модуля эл поля к давлению[атм]
    //p - давление[атм]

    return 17.03;
}

 //коэффициент электронно-ионной рекомбинации
double beta_ei(double x)
{
    //x - отношение модуля эл поля к давлению[атм]

      return 2.0E-07;
}

 //коэффициент ион-ионной рекомбинации
double beta_ii(double p)
{
    //p - давление[атм]
    double val = 0.0;

    val = 2.0E-07 + 5.4E-06*p;

    return val;
}

// вычисление тока проводимости в точке 
// (модель 6 - равновесная для воздуха)
void sub_sigma6(double p, double e, double qc,
                double *ce, double *ci, double *sg)
{
    //p - давление слоя
    //x_air - доля воздуха в смеси воздух + азот
    //e - модуль электрического поля
    //qc - источник вторичных электронов
    //ce, ci - концентрации электронов и ионов
    //sg - проводимость

    double qc1, ce1, ce2, ci1, ci2, ry, py, ay, by;
    double ue1, al1, gm1, gm2, bt1, bt2;
    double ht = spars_ht;

    ce1 = *ce; //концентрация вторичных электронов
    ci1 = *ci; //концентрация отрицательных ионов 
    qc1 = qc; //источник вторичных

    ue1 = u_air(e/p); //модуль скорости вторичных электронов

    al1 = alfa_air(e/p); //коэффициент таунсенда / p
    gm1 = gamma_air(e/p, p); //прилипание к кислороду / p
    gm2 = gammd_air(e/p, p); //отлипание от кислорода / p
    bt1 = beta_ei(e/p); //коэффициент электронно-ионной рекомбинации e + n2+ и e + o2+
    bt2 = beta_ii(p); //коэффициент ион-ионной рекомбинации o2- + n2+ ; o2- + o2+

    //вычисление концентрации вторичных электронов
    ry = p*al1*ue1 - p*gm1 - bt1*(ci1 + ce1);
    py = ry*ht;
    ay = exp(py);

    if (py < 1.0E-10 ) {
        by = ht;
    } else {
        by = -(1.0 - ay)/ry;
    }

    ce2 = ce1*ay + (qc1 + p*gm2*ci1)*by;
    ce1 = 0.5*(ce1 + ce2);

    //вычисление концентрации отрицательных ионов
    ry = -p*gm2 - bt2*(ci1 + ce1);
    py = ry*ht;
    ay = exp(py);

    if (py < 1.0E-10) {
        by = ht;
    } else {
        by = -(1.0 - ay)/ry;
    }

    ci2 = ci1*ay + (p*gm1*ce1)*by;

    *ce = ce2;
    *ci = ci2;

    //вычисление проводимости
    if (e == 0.0) {
        *sg = 0.0;
    } else {
        *sg = _qe*(ce1*ue1/e + (ce1 + 2.0*ci1)*6.E02/p);
    }
}

// вычисление тока проводимости в точке 
// (модель jump - для прыжка)
void sub_sigma_jump(double p, double e, double qc,
                double *ce, double *ci, double *sg)
{
    //p - давление слоя
    //x_air - доля воздуха в смеси воздух + азот
    //e - модуль электрического поля
    //qc - источник вторичных электронов
    //ce, ci - концентрации электронов и ионов
    //sg - проводимость

    double qc1, ce1, ci1;
    double ue1, gm1, bt2;

    ce1 = *ce; //концентрация вторичных электронов
    ci1 = *ci; //концентрация отрицательных ионов 
    qc1 = qc; //источник вторичных

    ue1 = u_air(e/p); //модуль скорости вторичных электронов
    gm1 = gamma_air(e/p, p); //прилипание к кислороду / p
    bt2 = beta_ii(p); //коэффициент ион-ионной рекомбинации o2- + n2+ ; o2- + o2+

    ce1 = qc1/gm1;
    ci1 = sqrt(qc1/bt2);

    *ce = ce1;
    *ci = ci1;

    //вычисление проводимости
    if (e == 0.0) {
        *sg = 0.0;
    } else {
        *sg = _qe*(ce1*ue1/e + (ce1 + 2.0*ci1)*6.E02/p);
    }
}