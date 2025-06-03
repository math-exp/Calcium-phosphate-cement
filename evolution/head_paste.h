// Time-stamp: <2025-06-02 18:22:08 shige>
#ifndef HEAD_PASTE
#define HEAD_PASTE

#include "../tools/matutil.h"

#include "../head_Sphi_IPc.h"
// 0: whole region: phi_bar+RND
// 1: only on the outer boundary: phi_int+phi_init+RND
// 2: only on the boundary of a small square: phi_bar+phi_init+RND 
// #define select_phi (0) 
// #define IPc (0.0) // concentration of IP6 (0: water, infty: max concentration)

// parameter settings in the paper
//
// select_phi=1
// 1: only on the outer boundary: phi_int+phi_init+RND
// 
// RK, L=1, N=200, D0=0.01, tmax=1, T0=0.8, eps=0.1
// dt=0.1/N^2, h=2L/N
// c: IPc=0, 0.02, 0.03
// phi_bar=0.5

#define select (1) // 0: Euler, 1: RK
#define N (128)
#define D0 (0.01)
// #define D0 (0.005)
#define tmax (1.0)
#define T0 (0.8 * tmax)
#define D(t) ( t < T0 ? D0 : 0.0 ) // step function: common in IP6 and water cases

// eps << 1 --> phase-separation is fast
// transition layer: order(eps)
// Note. MCF: V = -curvature <-- time scale of the order(eps^2 * t)
#define eps (0.1)
#define L (1.0) // [-L, L] x [-L, L]
#define h (2 * L / N)
#define phi_init ( select_phi == 0 ? 0 : 0.01 )
#define phi_bar (0.5)
// #define feps0(IPc) ( IPc > 0 ? log(1+IPc)/IPc : 1.0 )
#define feps0(IPc) ( 1.0 / ( 1.0 + IPc / 1000.0 ) )
#define feps(phi) ( phi * (1.0-phi) * (phi-(phi_bar * feps0(IPc))) / eps / eps )
// #define feps(phi) ( phi * (1.0-phi) * (phi-phi_bar) / eps / eps )

#define dt (0.1 / N / N)
#define datanum (400)

/**********************************************************************/
// random double num. in [a, b]
double RND( double a, double b );
// NBC for phi
void NBC( matrix phi );
// v --> Lap v (i=1, 2, ..., N; j=1, 2, ..., N)
void Lap( matrix v, matrix Lapv );
// (phi) --> F=F(phi) for solving dphi/dt=F(phi)
void get_F( matrix phi, double t, matrix Fphi );
// evolution by Euler for solving dphi/dt=F(phi)
// old (phi) --> new (phi)
void Euler( matrix old_phi, double t, matrix new_phi );
// evolution by Runge-Kutta for solving dphi/dt=F(phi)
// old (phi) --> new (phi)
void RK( matrix old_phi, double t, matrix new_phi );
// evolution: (old_phi, t) --> (new_phi)
void evolution( matrix old_phi, double t, matrix new_phi );
/**********************************************************************/
#endif

