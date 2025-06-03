#ifndef MATUTIL
#define MATUTIL

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifndef EXIT_FAILURE
#define EXIT_FAILURE (1)
#endif

#ifndef SCALAR
#define SCALAR double
#endif

typedef SCALAR *vector, **matrix;

/**********************************************************************/
void error(char *message);
vector newvec(int n);
matrix newmat(int nrow, int ncol);
vector new_vector(int n);
matrix new_matrix(int nrow, int ncol);
void free_vector(vector v);
void free_matrix(matrix a);
/**********************************************************************/
#endif
