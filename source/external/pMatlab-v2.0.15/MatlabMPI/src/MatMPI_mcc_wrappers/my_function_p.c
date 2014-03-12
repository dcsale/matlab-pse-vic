#include <stdio.h>
#include <string.h>
#include <math.h>
#include "matlab.h"
#include "multpkg.h"

/*
 * Function prototype; the MATLAB Compiler creates mlfMy_function
 *  from my_function.m
 */

int main (int argc, char *argv[]) /* Programmer written coded to call
mlfMy_function */
{
  mxArray *my_cpu;
  int rank;

  rank=(int)atoi(argv[1]);

  /* fprintf(stdout, "  rank: %d\r\n", rank); */

  multpkgInitialize();

  my_cpu=mlfScalar(rank);

  /* Call the mlfMy_function function. */
  mlfMy_function(my_cpu);

  multpkgTerminate();

  return(0);
}
