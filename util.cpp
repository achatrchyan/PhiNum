#ifndef _UTIL_CPP_
#define _UTIL_CPP_

#include <ctime>
#include <cmath>

#include "parameters.cpp"
#include "lattice.cpp"



int mod(int i,int N)
{
  if (i>=0)
    return i%N;
  else
    return (N+i)%N;
}

pseudo_double sqr(pseudo_double x)
{
  return x*x;
}

pseudo_double EucledianLength(pseudo_double *vector, int length)
{
  pseudo_double res=0.;
  for(int i=0;i<length;i++)
    res+=sqr(vector[i]);
  return sqrt(res);
}


int pos(int x, int y, int z)
{
	return Nc*x + Nc*N_t*y + Nc*N_t*N_t*z;
}

pseudo_double taufunc()
{
	return LocalLattice::tau;
}


#endif
