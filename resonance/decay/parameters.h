#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <complex>
#include <cmath>

#define VERBOSE						2
#define USE_CHEBYSHEV				1
#define USE_AZIMUTHAL_SYMMETRY		true

using namespace std;

const int n_s_pts = 21;
const int n_v_pts = 21;
const int n_zeta_pts = 21;

// Particle information
const int Maxparticle=400;
const int Maxdecaychannel=13;
const int Maxdecaypart=5;

// Momentum information
const int n_pT_pts = 15;
const int n_pphi_pts_preset = 36;
const int n_pphi_pts = ( USE_AZIMUTHAL_SYMMETRY ) ? 1 : n_pphi_pts_preset;
const int tmp_n_pY_pts = 101;
const int n_pY_pts = 1001;

// Misc.
const double resonanceThreshold = 0.0;

const int n_xi_pts = 101;
const int n_k_pts = 101;

const double xi_max = 5.0;
const double k_max = 5.0;
const double xi_min = -xi_max;
const double k_min = -k_max;

const double pT_min = 0.0;
const double pT_max = 4.0;
const double pphi_min = 0.0;
const double pphi_max = 2.0*M_PI;
const double Del_pY_min = -50.0;
const double Del_pY_max = 50.0;

//other misc.
extern int n_resonances;

#endif
