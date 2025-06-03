// Time-stamp: <2025-06-02 18:13:39 shige>
//
// make data directory
// 

#include "head_paste.h"

/*********************************************/
int main( void ){
  char   dir[100], filename[100];
  FILE   *fp;

  // create: making data directory
  if( select == 0 ) 
    sprintf( filename, "EL-phi0=%d-N=%d-IPc=%.2f", select_phi, N, IPc );
  else if( select == 1 ) 
    sprintf( filename, "RK-phi0=%d-N=%d-IPc=%.2f", select_phi, N, IPc );
  else{
    printf( "Input select = 0 (Euler) or 1 (RK)\n");
    exit(1);
  }
  sprintf( dir, "data/%s", filename );

  fp = fopen( "../mkdatadir.sh", "w" );
  fprintf( fp, "#!/bin/bash\n" );
  fprintf( fp, "mkdir %s\n", dir );
  fclose( fp );

  return 0; 
}
