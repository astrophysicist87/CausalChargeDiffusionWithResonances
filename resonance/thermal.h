#ifndef THERMAL_H
#define THERMAL_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <complex>
#include <cmath>

using namespace std;

#include "decay/parameters.h"
#include "decay/readindata.h"

extern vector<double> pT_pts, pT_wts;
extern vector<double> pphi_pts, pphi_wts;
extern vector<double> Del_pY_pts, Del_pY_wts;
extern double * x_pts, * x_wts;

extern const int n_r_pts;
extern const double TFO;

//define some parameters for the exact emission function
const double Rad = 5.0, Del_tau = 1.0, tau0 = 5.0, etaf = 0.6;

inline double Hfactor(double r, double tau)
{
	return (
			exp( -r*r/(2.0*Rad*Rad) - (tau-tau0)*(tau-tau0)/(2.0*Del_tau*Del_tau) ) / (M_PI*Del_tau)
			);
}

inline double Hfactor(double r)
{
	return (
			exp( -r*r/(2.0*Rad*Rad) )
			);
}

inline double eta_t(double r)
{
	return ( etaf*r/Rad );
}

double Cal_dN_dypTdpTdphi_toy_func(
		int local_pid, vector<readindata::particle_info> * all_particles,
		double pT_loc, double pphi_loc, double pY_loc )
{
	const double hbarC = 0.197327053;

	// set particle information
	double degen = (*all_particles)[local_pid].gspin;
	double localmass = (*all_particles)[local_pid].mass;

	// set some freeze-out surface information that's constant the whole time
	double prefactor = degen/sqrt(8.0*M_PI*M_PI*M_PI)/M_PI/hbarC;
	double Tdec = TFO;	//back to MeV
	double one_by_Tdec = 1./Tdec;

	double ch_pY = cosh(pY_loc);
	double sh_pY = sinh(pY_loc);

	double result = 0.0;

	for (int ir = 0; ir < n_r_pts; ++ir)
	{
		double mT = sqrt(pT_loc*pT_loc+localmass*localmass);

		//set r-point inside pT loop, since optimal distribution of integration points
		// depends on value of MT
		double rmin = 0.0, rmax = 25.0 * Rad / sqrt(1.0 + mT*one_by_Tdec*etaf*etaf);
		double hw = 0.5 * (rmax - rmin), cen = 0.5 * (rmax + rmin);
		double rpt = cen + hw * x_pts[ir];

		double local_H = Hfactor(rpt);
		double ch_eta_t = cosh(eta_t(rpt));
		double sh_eta_t = sinh(eta_t(rpt));

		double alpha = mT * one_by_Tdec * ch_eta_t;
		double beta = one_by_Tdec * pT_loc * sh_eta_t;
		if (alpha > 700.0)
			continue;

		double K1 = 2.0 * gsl_sf_bessel_K1(alpha);
		double I0 = gsl_sf_bessel_I0(beta);

		double S_p_with_weight = hw*x_wts[ir]*mT*tau0*rpt*prefactor*local_H*I0;

		result += S_p_with_weight * K1;
	}

	return (result);
}

#endif
