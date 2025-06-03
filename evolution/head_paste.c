// Time-stamp: <2025-06-02 18:07:11 shige>

#include "head_paste.h"

// random double num. in [a, b]
double RND( double a, double b ){
  double r; 
  static int flag;
  
  /* if( flag == 0 ){ */
  /*   srand((unsigned int)time(NULL)); */
  /*   flag = 1; */
  /* } */

  // random double num. in [a, b]
  r = a + ( (double) rand() / (double) RAND_MAX ) * ( b - a );  
  
  return r; 
}

// get volume
double volume( matrix u ){
  double vol; 
  int i, j;

  vol = 0;
  for( i=1; i<=N; i++ ){
    for( j=1; j<=N; j++ ){
      vol += u[i][j] * h * h;
    }
  }

  return vol;
}

// NBC for phi
void NBC( matrix phi ){
  int i, j;

  for( i=0; i<=N+1; i++ ){
    phi[i][0] = phi[i][1];
    phi[i][N+1] = phi[i][N];
  }
  for( j=0; j<=N+1; j++ ){
    phi[0][j] = phi[1][j];
    phi[N+1][j] = phi[N][j];
  }
}

// v --> Lap v (i=1, 2, ..., N; j=1, 2, ..., N)
void Lap( matrix v, matrix Lapv ){
  double Lapx_v, Lapy_v; 
  int i, j;

  for( i=1; i<=N; i++ ){
    for( j=1; j<=N; j++ ){
      Lapx_v = ( v[i+1][j] - 2 * v[i][j] + v[i-1][j] ) / h / h;
      Lapy_v = ( v[i][j+1] - 2 * v[i][j] + v[i][j-1] ) / h / h;
      Lapv[i][j] = Lapx_v + Lapy_v;  
    }
  }
}

// (phi) --> F=F(phi) for solving dphi/dt=F(phi)
void get_F( matrix phi, double t, matrix Fphi ){
  matrix Lapphi;
  int i, j; 

  Lapphi = new_matrix( N + 2, N + 2 );

  // v --> Lap v (i=1, 2, ..., N; j=1, 2, ..., N)
  Lap( phi, Lapphi );

  for( i=1; i<=N; i++ ){
    for( j=1; j<=N; j++ ){
      Fphi[i][j] = D(t) * Lapphi[i][j] + feps( phi[i][j] );
    }
  }

  free_matrix( Lapphi );
}

// evolution by Euler for solving dphi/dt=F(phi)
// old (phi) --> new (phi)
void Euler( matrix old_phi, double t, matrix new_phi ){
  matrix Fphi; 
  int i, j; 

  Fphi = new_matrix( N + 2, N + 2 );

  get_F( old_phi, t, Fphi );
  for( i=1; i<=N; i++ ){
    for( j=1; j<=N; j++ ){
      new_phi[i][j] = old_phi[i][j] + Fphi[i][j] * dt;
    }
  }

  free_matrix( Fphi );
}

// evolution by Runge-Kutta for solving dphi/dt=F(phi)
// old (phi) --> new (phi)
void RK( matrix old_phi, double t, matrix new_phi ){
  matrix Fphi1, Fphi2, Fphi3, Fphi4;
  matrix phi1, phi2, phi3;
  int i, j; 

  Fphi1 = new_matrix( N + 2, N + 2 );
  Fphi2 = new_matrix( N + 2, N + 2 );
  Fphi3 = new_matrix( N + 2, N + 2 );
  Fphi4 = new_matrix( N + 2, N + 2 );
  phi1 = new_matrix( N + 2, N + 2 );
  phi2 = new_matrix( N + 2, N + 2 );
  phi3 = new_matrix( N + 2, N + 2 );
  
  // step 1
  get_F( old_phi, t, Fphi1 );

  // step 2
  for( i=1; i<=N; i++ ){
    for( j=1; j<=N; j++ ){
      phi1[i][j] = old_phi[i][j] + Fphi1[i][j] * dt / 2;
    }
  }
  NBC( phi1 );
  get_F( phi1, t + dt / 2, Fphi2 );

  // step 3
  for( i=1; i<=N; i++ ){
    for( j=1; j<=N; j++ ){
      phi2[i][j] = old_phi[i][j] + Fphi2[i][j] * dt / 2;
    }
  }
  NBC( phi2 );
  get_F( phi2, t + dt / 2, Fphi3 );

  // step 4
  for( i=1; i<=N; i++ ){
    for( j=1; j<=N; j++ ){
      phi3[i][j] = old_phi[i][j] + Fphi3[i][j] * dt;
    }
  }
  NBC( phi3 );
  get_F( phi3, t + dt, Fphi4 );

  // final step
  for( i=1; i<=N; i++ ){
    for( j=1; j<=N; j++ ){
      new_phi[i][j] = old_phi[i][j] 
	+ ( Fphi1[i][j] + 2 * Fphi2[i][j] + 2 * Fphi3[i][j] + Fphi4[i][j] ) * dt / 6;
    }
  }

  free_matrix( Fphi1 );
  free_matrix( Fphi2 );
  free_matrix( Fphi3 );
  free_matrix( Fphi4 );
  free_matrix( phi1 );
  free_matrix( phi2 );
  free_matrix( phi3 );
}

// evolution: (old_phi, t) --> (new_phi)
void evolution( matrix old_phi, double t, matrix new_phi ){

  // Euler or RK
  if( select == 0 ) Euler( old_phi, t, new_phi );
  else if( select == 1 ) RK( old_phi, t, new_phi );

  NBC( new_phi );
}
/*******************************************************************/
