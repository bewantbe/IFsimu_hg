#include "stdafx.h"
#include "random.h"

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran0(long *idum)
{
  long k;
  double ans;

  *idum ^= MASK;
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  ans=AM*(*idum);
  *idum ^= MASK;
  return ans;
}

double ran1(long *idum, long*iy, long**iv, int index)
{
  int j;
  long k;
//  static long iy=0;
//  static long iv[NTAB];
  double temp;

  if (*idum <= 0 || !iy[index]) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    for (j=NTAB+7; j>=0; j--) {
      k=(*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[index][j] = *idum;
    }
    iy[index]=iv[index][0];
  }
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  j=iy[index]/NDIV;
  iy[index]=iv[index][j];
  iv[index][j] = *idum;
  if ((temp=AM*iy[index]) > RNMX) return RNMX;
  else return temp;
}

double gasdev(long *idum, int &iset, double &gset)
{
//  static int iset=0;
//  static double gset;
  double fac,rsq,v1,v2;

  if  (iset == 0) {
    do {
      v1=2.0*ran0(idum)-1.0;
      v2=2.0*ran0(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK
#undef NDIV
#undef EPS
#undef RNMX
