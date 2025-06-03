// Time-stamp: <2025-06-02 18:12:59 shige>
//
// Authors: Yu ICHIDA, Rika YAMADA, Shiori KATO, Yuki KAMAYA, Minami KOSUGE, 
// Mamoru AIZAWA, Takashi Okuda SAKAMOTO, Shigetoshi YAZAKI
// Title: A simple mathematical model for evaluation of non-fragmentation property of 
// injectable calcium-phosphate cement
// Journal: Scientific Reports
// 
// Allen-Cahn type PDE of the order parameter phi = phi(x, y, t)
// 
// d/dt(phi) = D(t).Lap(phi) + f(phi)
// Lap = (d/dx)^2 + (d/dy)^2
// f(phi) = phi.(1 - phi).(phi - phi_bar(c)) / eps^2
// phi_bar(c): monotone dicreasing function w.r.t. c (concentration of IP6)
// (phi_bar(0)=0.5, c=0: water)
// 
// Domain: Omega=[-L, L]x[-L, L], L=1
// Boundary condition: NBC
// Space increment: h = dx = dy = 2 * L / N
// 
// Numerical method: the method of lines
// 

#include "head_paste.h"

/*********************************************/
int main( void ){
  matrix phi;
  vector x, y;
  double t, phimin, phimax, frame, noise_range;
  int    i, j, count, dummyi, cbflag;
  char   file[100], dir[100], filename[100];
  FILE   *fp, *fp1;

  // open matrices
  phi = new_matrix( N + 2, N + 2 );
  x = new_vector( N + 2 );
  y = new_vector( N + 2 );

  // staggered grid
  for( i=0; i<=N+1; i++ ) x[i] = -L + ( i - 0.5 ) * h;
  for( j=0; j<=N+1; j++ ) y[j] = -L + ( j - 0.5 ) * h;

  // initial data
  t = 0; 
  count = 0; 

  phimin = 1000;
  phimax = -1000;
  noise_range = h * 0.05; 
  for( i=1; i<=N; i++ ){
    for( j=1; j<=N; j++ ){
      if( select_phi == 0 ){
	// whole region: phi_bar+RND
	phi[i][j] = phi_bar + phi_init + RND( -noise_range, noise_range );
      }
      else if( select_phi == 1 ){
	// only on the outer boundary: phi_int+phi_init+RND
	if( i==1 || i==N || j==1 || j==N ) 
	  phi[i][j] = phi_bar + phi_init + RND( -noise_range, noise_range );
	else phi[i][j] = phi_bar + RND( -noise_range, noise_range );
      }
      else if( select_phi == 2 ){
	// only on the boundary of a small square: phi_bar+phi_init+RND
	frame = 0.2; // 0 <= frame < 0.5
	if( ( i==(int)(frame*N) && (int)(frame*N) <= j && j<=(int)((1-frame)*N) ) 
	    || ( i==(int)((1-frame)*N) && (int)(frame*N) <= j && j<=(int)((1-frame)*N) )  
	    || ( j==(int)(frame*N) && (int)(frame*N) <= i && i<=(int)((1-frame)*N) ) 
	    || ( j==(int)((1-frame)*N) && (int)(frame*N) <= i && i<=(int)((1-frame)*N) ) )
	  phi[i][j] = phi_bar + phi_init + RND( -noise_range, noise_range );
	else phi[i][j] = phi_bar + RND( -noise_range, noise_range );
      }
      else if( select_phi == 3 ){
	// increments on a circle
	if( sqrt( x[i] * x[i] + y[j] * y[j] ) < L / 2 
	    && ( fabs( y[j] ) > L / 50 || x[i] < - L / 5 ) )
	  phi[i][j] = phi_bar + phi_init + RND( -noise_range, noise_range );
	else
	  phi[i][j] = 0.0;
      }

      if( phi[i][j] < phimin ) phimin = phi[i][j];
      else if( phi[i][j] > phimax ) phimax = phi[i][j];
    }
  }

  // NBC
  NBC( phi );

  // data directory
  if( select == 0 ) 
    sprintf( filename, "EL-phi0=%d-N=%d-IPc=%.2f", select_phi, N, IPc );
  else if( select == 1 ) 
    sprintf( filename, "RK-phi0=%d-N=%d-IPc=%.2f", select_phi, N, IPc );
  else{
    printf( "Input select = 0 (Euler) or 1 (RK)\n");
    exit(1);
  }
  sprintf( dir, "data/%s", filename );

  // save parameter
  sprintf( file, "../%s/%04d-param.dat", dir, count );
  fp = fopen( file, "w" );
  fprintf( fp, "select=%d\n", select );
  fprintf( fp, "IPc=%g\n", IPc );
  fprintf( fp, "D0=%g\n", D0 );
  fprintf( fp, "tmax=%g\n", tmax );
  fprintf( fp, "eps=%g\n", eps );
  fclose( fp );

  // save the initial data
  sprintf( file, "../%s/%04d.dat", dir, count );
  fp = fopen( file, "w" );
  for( i=1; i<=N; i++ ){
    for( j=1; j<=N; j++ ){
      fprintf( fp, "%18.15f %18.15f %18.15f\n", x[i], y[j], phi[i][j] );
    }
    fprintf( fp, "\n" );
  }
  fclose( fp );

  // parameters data 
  sprintf( file, "../%s/0000-data.dat", dir );
  fp = fopen( file, "w" );
  fprintf( fp, "%04d %18.10f %d\n", count, t, N );
  fclose( fp );

  printf( "------\n" );
  if( IPc > 0 ) printf( "tmax=%g, datanum=%d, IPc=%g (IP6), eps=%g\n", tmax, datanum, IPc, eps);
  else printf( "tmax=%g, datanum=%d, IPc=0 (water), eps=%g\n", tmax, datanum, eps );
  // printf( "dT01=%5.4f; Tw0=%5.4f, Tw1=%5.4f, Tp0=%5.4f, Tp1=%5.4f\n", dT01, Tw0, Tw1, Tp0, Tp1 );
  printf( "------\n" );
  printf( "%04d, t=%5.3f, N=%d, phimin=%5.4f, phimax=%5.4f\n", 
	  count, t, N, phimin, phimax );
  count++;

  // evolution
  while( t <= tmax + 0.1 * dt ){

    // evolution: (old_phi, t) --> (new_phi)
    evolution( phi, t, phi );
      
    t += dt; 

    // save data: (x, y, phi)
    if( tmax / (double) datanum * count <= t && t <= tmax + 0.1 * dt ){
      sprintf( file, "../%s/%04d.dat", dir, count );
      fp = fopen( file, "w" );
      phimin = 1000;
      phimax = -1000;
      for( i=1; i<=N; i++ ){
	for( j=1; j<=N; j++ ){
	  fprintf( fp, "%18.15f %18.15f %18.15f\n", x[i], y[j], phi[i][j] );
	  if( phi[i][j] < phimin ) phimin = phi[i][j];
	  else if( phi[i][j] > phimax ) phimax = phi[i][j];
	}
	fprintf( fp, "\n" );
      }
      fclose( fp );

      // parameters data 
      sprintf( file, "../%s/0000-data.dat", dir );
      fp = fopen( file, "a" );
      fprintf( fp, "%04d %18.10f %d\n", count, t, N );
      fclose( fp );

      printf( "%04d, t=%5.3f, N=%d, phimin=%5.4f, phimax=%5.4f\n", 
	      count, t, N, phimin, phimax );

      count++; 
    }
  }
  
  // gnuplot
  sprintf( file, "../%s/0000-data.gnu", dir );
  fp = fopen( file, "w" );
  fprintf( fp, "unset key\n" );
  fprintf( fp, "set border\n" );
  fprintf( fp, "unset xtics\n" );  
  fprintf( fp, "unset ytics\n" );  
  fprintf( fp, "set xrange [%g:%g]\n", -L, L );
  fprintf( fp, "set yrange [%g:%g]\n", -L, L );
  fprintf( fp, "set zrange [0:2]\n" );
  fprintf( fp, "set ticslevel 0\n" );
  fprintf( fp, "set isosamples 100\n" );
  fprintf( fp, "set hidden3d\n" );
  fprintf( fp, "set view equal xy\n" );
  fprintf( fp, "set view map\n" );
  
  //
  // make jpeg files
  // 
  fprintf( fp, "set term jpeg\n" );
  sprintf( file, "../%s/0000-data.dat", dir );

  // Initial data: i = 0
  // When i = 0, zoom in the figure for noise-visible
  // (If cbrange = [0,1], the figure is drawn in one-color only (cannot see noise). 
  if( select_phi > 0 ){
    fprintf( fp, "set cbrange [%g:%g]\n", 
	     phi_bar - noise_range, phi_bar + phi_init + noise_range );
    fprintf( fp, "set palette defined (%g 'black', %g 'white', %g 'yellow')\n", 
	     phi_bar - noise_range, phi_bar, phi_bar + phi_init + noise_range );
  }
  else{
    fprintf( fp, "set cbrange [%g:%g]\n", 
	     phi_bar + phi_init - noise_range, phi_bar + phi_init + noise_range );
    fprintf( fp, "set palette defined (%g 'black', %g 'white', %g 'yellow')\n", 
	     phi_bar + phi_init - noise_range, phi_bar + phi_init, phi_bar + phi_init + noise_range );
  }
  // fprintf( fp, "set palette defined (%g 'blue', %g 'white', %g 'red')\n", 
  // fprintf( fp, "set palette defined (%g 'gray', %g 'black', %g 'white')\n", 

  i = 0; 
  t = 0.0;
  if( IPc > 0 ) fprintf( fp, "set title 'time=%5.3f, IP6_c=%g'\n", t, IPc );
  else fprintf( fp, "set title 'time=%5.3f, IP6_c=0 (water)'\n", t );
  fprintf( fp, "set output '%s/0000_init.jpeg'\n", dir );
  fprintf( fp, "splot '%s/%04d.dat' u 1:2:3 w pm3d\n", dir, i );  

  // all data: i = 0, ..., count-1
  fprintf( fp, "set cbrange [0:1]\n" );
  fprintf( fp, "set palette defined (0 'black', 0.5 'white', 1 'yellow')\n" );
  // fprintf( fp, "set palette defined (0 'gray', 0.5 'black', 1 'white')\n" );
  //fprintf( fp, "set palette defined (0 'blue', 0.5 'white', 1 'red')\n" ); 
  fp1 = fopen( file, "r" );
  for( i=0; i<count; i++ ){
    fscanf( fp1, "%d %lf %d\n", &dummyi, &t, &dummyi );
    if( IPc > 0 ) fprintf( fp, "set title 'time=%5.3f, IP6_c=%g'\n", t, IPc );
    else fprintf( fp, "set title 'time=%5.3f, IP6_c=0 (water)'\n", t );
    fprintf( fp, "set output '%s/%04d.jpeg'\n", dir, i ); 
    fprintf( fp, "splot '%s/%04d.dat' u 1:2:3 w pm3d\n", dir, i );  
  }
  fclose( fp1 );

  // last data: i = count - 1
  i = count - 1; 
  fprintf( fp, "set output '%s/0000_final.jpeg'\n", dir );  
  fprintf( fp, "splot '%s/%04d.dat' u 1:2:3 w pm3d\n", dir, i );  

  fprintf( fp, "unset output\n" );

  //
  // make eps files
  // 
  fprintf( fp, "set terminal epscairo enhanced\n" );
  sprintf( file, "../%s/0000-data.dat", dir );

  // Initial data: i = 0
  // When i = 0, zoom in the figure for noise-visible
  // (If cbrange = [0,1], the figure is drawn in one-color only (cannot see noise). 
  if( select_phi > 0 ){
    fprintf( fp, "set cbrange [%g:%g]\n", 
	     phi_bar - noise_range, phi_bar + phi_init + noise_range );
    fprintf( fp, "set palette defined (%g 'black', %g 'white', %g 'yellow')\n", 
	     phi_bar - noise_range, phi_bar, phi_bar + phi_init + noise_range );
  }
  else{
    fprintf( fp, "set cbrange [%g:%g]\n", 
	     phi_bar + phi_init - noise_range, phi_bar + phi_init + noise_range );
    fprintf( fp, "set palette defined (%g 'black', %g 'white', %g 'yellow')\n", 
	     phi_bar + phi_init - noise_range, phi_bar + phi_init, phi_bar + phi_init + noise_range );
  }

  i = 0; 
  t = 0.0;
  if( IPc > 0 ) fprintf( fp, "set title 'time=%5.3f, IP6_c=%g'\n", t, IPc );
  else fprintf( fp, "set title 'time=%5.3f, IP6_c=0 (water)'\n", t );
  fprintf( fp, "set output '%s/0000_init.eps'\n", dir );
  fprintf( fp, "splot '%s/%04d.dat' u 1:2:3 w pm3d\n", dir, i );  

  // all data: i = 0, ..., count-1
  fprintf( fp, "set cbrange [0:1]\n" );
  fprintf( fp, "set palette defined (0 'black', 0.5 'white', 1 'yellow')\n" );
  fp1 = fopen( file, "r" );
  for( i=0; i<count; i++ ){
    fscanf( fp1, "%d %lf %d\n", &dummyi, &t, &dummyi );
    if( IPc > 0 ) fprintf( fp, "set title 'time=%5.3f, IP6_c=%g'\n", t, IPc );
    else fprintf( fp, "set title 'time=%5.3f, IP6_c=0 (water)'\n", t );
    fprintf( fp, "set output '%s/%04d.eps'\n", dir, i ); 
    fprintf( fp, "splot '%s/%04d.dat' u 1:2:3 w pm3d\n", dir, i );  
  }
  fclose( fp1 );

  // last data: i = count - 1
  i = count - 1; 
  fprintf( fp, "set output '%s/0000_final.eps'\n", dir );  
  fprintf( fp, "splot '%s/%04d.dat' u 1:2:3 w pm3d\n", dir, i );  

  fprintf( fp, "unset output\n" );

  //
  // make gif animation
  // 
  fprintf( fp, "set terminal gif animate\n" );
  sprintf( file, "../%s/0000-data.dat", dir );
  fprintf( fp, "set output '%s/0-%s-data-phi.gif'\n", dir, filename );

  fprintf( fp, "set cbrange [0:1]\n" );
  fprintf( fp, "set palette defined (0 'black', 0.5 'white', 1 'yellow')\n" );
  fp1 = fopen( file, "r" );
  for( i=0; i<count; i++ ){
    fscanf( fp1, "%d %lf %d\n", &dummyi, &t, &dummyi );
    if( IPc > 0 ) fprintf( fp, "set title 'time=%5.3f, IP6_c=%g'\n", t, IPc );
    else fprintf( fp, "set title 'time=%5.3f, IP6_c=0 (water)'\n", t );
    fprintf( fp, "splot '%s/%04d.dat' u 1:2:3 w pm3d\n", dir, i );  
  }
  fclose( fp1 );
  fprintf( fp, "unset output\n" );

  //
  // create movie files
  //
  fp = fopen( "../gnuplot.sh", "w" );
  fprintf( fp, "#!/bin/bash\n" );
  fprintf( fp, "gnuplot %s/0000-data.gnu\n", dir );
  // fprintf( fp, "ffmpeg -y -i %s/%s %s/0-%s-data-phi.mp4\n", dir, "%04d.jpeg", dir, filename );
  fprintf( fp, "ffmpeg -y -i %s/%s %s/0-%s-data-phi.mpeg\n", dir, "%04d.jpeg", dir, filename );
  fprintf( fp, "#open -n /Applications/VLC.app %s/0-%s-data.mp4\n", dir, filename );
  fprintf( fp, "#open -n /Applications/VLC.app %s/0-%s-data.mpeg\n", dir, filename );
  fprintf( fp, "\n" );

  free_matrix( phi );
  free_vector( x );
  free_vector( y );

  return 0; 
}
