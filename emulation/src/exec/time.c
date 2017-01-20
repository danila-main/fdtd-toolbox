
/*! \file 
* \brief Системное время.
*
* Функция возвращает системное время в секундах.
* \author Крюков А.А. anton.krv@gmail.com
*/

#if defined(WIN32)
#include <time.h>
#include <windows.h>
#elif defined(LINUX)
#include <sys/time.h>
#else
#include <sys/time.h>
#endif

#include <stdio.h>

/*!
* \brief Системное время в секундах.
*
* Возвращает системное время в секундах, используется для оценки скорости работы вычислительных алгоритмов.
* \return Системное время в секундах.
*/
double dwalltime()
{
    double t = 0;

#if defined(WIN32)

    #if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
    #define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
    #else
    #define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
    #endif

    FILETIME ft;
    unsigned __int64 tmpres = 0;

    GetSystemTimeAsFileTime(&ft);

    tmpres |= ft.dwHighDateTime;
    tmpres <<= 32;
    tmpres |= ft.dwLowDateTime;

    tmpres /= 10;
    tmpres -= DELTA_EPOCH_IN_MICROSECS;

    t = (long)(tmpres / 1000000UL) + (long)(tmpres % 1000000UL) / 1.0e+6;
#elif defined(LINUX)
    struct timeval tval;
    gettimeofday(&tval,NULL);
    t = (double)tval.tv_sec + (double)tval.tv_usec / 1.0e+6;
#else
    struct timeval tval;
    gettimeofday(&tval,NULL);
    t = (double)tval.tv_sec + (double)tval.tv_usec / 1.0e+6;
#endif

  return t;
}
