extern "C" void kkp_(int &ih, int &iset, double &x, double &q, double dh[10]);
void kkp(int ih, int iset, double x, double q, double dh[10]);

//  Wrapper for kkp, to allow use as a TF1
double kkp_func(double *x, double *par);

/*
from the FORTRAN source:
c=====================================================================
c
c     ------------------------------------------------------------
c     Fragmentation functions for: Pions, Kaons, Protons, Neutrons
c            (includes mass-threshholds for c and b quarks)
c
c     Reference: B.A.Kniehl, G.Kramer, B.Potter, NPB582 (2000) 514
c     ------------------------------------------------------------
c
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
c
c     x    = longitudinal-momentum fraction
c     qs   = fragmentation scale (in GeV)
c
c     Parton label:
c     0    1    2    3    4    5    6    7    8     9    10
c     g    u   ubar  d   dbar  s   sbar  c   cbar   b   bbar
c
c     Lambda_QCD (in GeV):
c     0.088 in LO
c     0.213 in NLO
c
c=====================================================================
*/
