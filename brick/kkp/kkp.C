#include "kkp.h"
#include <iostream>

using std::cout;
using std::endl;

void kkp(int ih, int iset, double x, double q, double dh[10]) {
  kkp_(ih,iset,x,q,dh); // call fortran function
}


/*
  Wrapper for kkp, to allow use as a TF1
*/
double kkp_func(double *x, double *par) {
  // Note: par[0] is Q
  //       par[1] is parton flavour (integer)
  //     Parton label:
  //     0    1    2    3    4    5    6    7    8     9    10
  //     g    u   ubar  d   dbar  s   sbar  c   cbar   b   bbar
  //       par[2] is hadron type ih
/*
c     ih, iset, x, qs are input; dh is output.
c     ih   = 1 : (pi^+ + pi^-)  /2
c     ih   = 2 : (K^+ + K^-)    /2
c     ih   = 3 : (K^0 + K^0_bar)/2
c     ih   = 4 : (p + p_bar)    /2
c     ih   = 5 : (pi^0)
c     ih   = 6 : (n + n_bar)    /2
c     ih   = 7 : (h^+ + h^-)         [as sum of pions, kaons and protons]
c
c     iset = 0 : LO
c     iset = 1 : NLO
*/
  //
  double dh[10];
  if (par[0] == 0)
    cout << "ERROR: Q==0 in kkp_func" << endl;
  int ih = int(par[2]+0.001);
  if (ih < 1 || ih > 7) {cout << "kkp_func: invalid hadron type: " << ih << endl; ih = 7; }
  // 0 is LO; 1 is NLO
  kkp(ih,0,x[0],par[0],dh);
  int part_type = int(par[1]+0.0001);
  return dh[part_type]; 
}
