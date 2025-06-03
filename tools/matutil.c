#include "matutil.h"

/*******************************************************************/
void error(char *message){
  fprintf(stderr, "\n%s\n\n", message); 
  exit(EXIT_FAILURE);
}

vector newvec(int n){
  return malloc(sizeof(SCALAR) * n);
}

matrix newmat(int nrow, int ncol){
  int    i;
  matrix a;
  
  a = malloc((nrow + 1) * sizeof(void *));
  if (a == NULL) return NULL;
  for(i = 0; i < nrow; i++){
    a[i] = malloc(sizeof(SCALAR) * ncol);
    if ( a[i] == NULL) {
      while (--i >= 0) free(a[i]);
      free(a);
      return NULL;
    }
  }
  a[nrow] = NULL;
  return a;
}

vector new_vector(int n){
  vector v;
  
  v = newvec(n);
  if( v == NULL ) error("LACK of AVAILABLE MEMORY! in new_vector()");
  return v;
}

matrix new_matrix(int nrow, int ncol){
  matrix a;
  
  a = newmat(nrow, ncol);
  if(a == NULL) error("LACK of AVAILABLE MEMORY! in new_matrix()");
  return a;
}

void free_vector(vector v){
  free(v);
}

void free_matrix(matrix a){
  matrix b;

  b = a;
  while ( *b != NULL) free(*b++);
  free(a);
}
/*******************************************************************/
