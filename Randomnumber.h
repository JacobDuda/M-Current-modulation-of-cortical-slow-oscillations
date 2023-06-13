#ifndef randomnumber_H
#define randomnumber_H

// #include "Variables.h"
#include <iostream>
#include <cmath>


using namespace std;

const int numbermonitoringneurons=5;


float ran(long *idum);	/*random number*/
/*************************************************/
/*Return Gaussian distributed numbers*/
double gaussian_generator(double mu, double sigma, long seed)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
  if (call == 1){
    call = !call;
    return (mu + sigma * (double) X2);
  }
  do{
    U1 = -1.0 + ran(&seed) * 2;
    U2 = -1.0 + ran(&seed) * 2;
    W = pow (U1, 2) + pow (U2, 2);
  }
  while (W >= 1 || W == 0);
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
  call = !call;
  return (mu + sigma * (double) X1);
}
/*************************************************/



//%%%%%%%%%%%%%%%%%%%%
//Random Number:
#define IM1	2147483647
#define IM2	2147483647
#define AM (1.0/IM1)
#define AM1	(1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 16807
#define IA2 40692
#define IQ1	127773
#define IQ2 52774
#define IR1	2836
#define IR2 3791
#define NTAB 32
#define NTAB1	32
#define NDIV (1+IMM1/NTAB)
#define NDIV1	(1+(IM1-1)/NTAB1)
#define EPS 1.2e-7
#define EPS1	1.2e-7
#define RNMX (1.0-EPS)
#define RNMX1	(1.0-EPS1)


float ran(long *idum)
{

int j;
long k;
static long idum2 = 123456789, iy = 0, iv[NTAB];
float temp;

if(*idum <= 0)
{
 if(-(*idum) < 1)   *idum = 1;
 else   *idum = -(*idum);
 idum2 = (*idum);
  for(j = NTAB+7;j >= 0;j--)
   {
    k = (*idum)/IQ1;
    *idum = IA1*(*idum-k*IQ1)-k*IR1;
    if(*idum < 0)   *idum += IM1;
    if(j < NTAB)   iv[j] = *idum;
   }
 iy = iv[0];
}

k = (*idum)/IQ1;
*idum = IA1*(*idum-k*IQ1)-k*IR1;
if(*idum < 0)   *idum += IM1;

k = idum2/IQ2;
idum2 = IA2*(idum2-k*IQ2)-k*IR2;
if(idum2 < 0)   idum2 += IM2;

j = iy/NDIV;
iy = iv[j]-idum2;
iv[j] = *idum;
if(iy < 1)   iy += IMM1;

if((temp = AM*iy) > RNMX)   return RNMX;
else   return temp;
}

#endif
