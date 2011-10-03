#ifndef _RAN0_H_
#define _RAN0_H_

// one method to generate a random number which obeys the homogeneous
// distribution in [0,1]
double ran0(long *idum);

// another method to generate a random number which also obeys the
// homogeneous distribution in [0,1]
double ran1(long *idum, long*iy, long**iv, int index);

// this generates a random number with normal distribution, although it
// is not used in this code, but may be useful in the future
double gasdev(long *idum, int *iset, double *gset);

#endif
