#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <complex>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>

using namespace std;

#include "decay/parameters.h"
#include "lib.h"
#include "decay/readindata.h"
#include "gauss/gauss_quadrature.h"
#include "decay/decay.h"
#include "defs.h"


//space-time integration grid
const int n_tau_pts = 51;
const int n_r_pts = 101;
const int n_phi_pts = 51;
double * tau_pts, * tau_wts;
double * x_pts, * x_wts;
double * phi_pts, * phi_wts;

vector<double> xi_pts, xi_wts;
vector<double> k_pts, k_wts;

int n_resonances = -1;
const double AT = 1.0, tauFINAL = 10.0, TFO = 150.0 / 197.33;
const double chiQ = 2.0/3.0;

vector<complex<double> > scriptFn_vector;
vector<double> thermal_resonance_spectra_re, thermal_resonance_spectra_im;
vector<double> full_resonance_spectra_re, full_resonance_spectra_im;
vector<double> pT_pts, pT_wts;
vector<double> pphi_pts, pphi_wts;
vector<double> Del_pY_pts, Del_pY_wts;

//function prototypes
void set_up_misc();
void set_scriptFn_vector();
double f(double x, double alpha);
double get_integral(double kpt, double alpha);
complex<double> thermal_FT_spectra(readindata::particle_info * particle, double kpt, double PT, double PPHI, double PY);
vector<double> * copy_spectra(vector<double> * spectra, int ik);
double Cal_dN_dypTdpTdphi_toy_func(
		int local_pid, vector<readindata::particle_info> * all_particles,
		double pT_loc, double pphi_loc, double pY_loc );

///////////////////////////////////////////////////////////////
//////////////////////
//Start of main
//////////////////////

double g(double x, double width)
{
	return (sin(x)*exp(-x*width));
}

int main(int argc, char *argv[])
{
	/*pT_pts = vector<double>(n_pT_pts);
	pT_wts = vector<double>(n_pT_pts);

	double width = 10.0;
	//gauss_quadrature(n_pT_pts, 5, 0.0, 0.0, 0.0, width, pT_pts, pT_wts);
	gauss_quadrature(101, 1, 0.0, 0.0, 0.0, 5.0/width, pT_pts, pT_wts);

	double integral = 0.0;
	for (int ipt = 0; ipt < n_pT_pts; ++ipt)
		integral += pT_wts[ipt]*g(pT_pts[ipt], width);

	cout << "result = " << integral << endl;

	if (1) exit(1);*/

	//initializes a few things
	set_up_misc();
	gsl_set_error_handler_off();

	//load resonance information
	const int particle_idx = 1;		// pi^+
	vector<readindata::particle_info> all_particles;
	vector<int> chosen_resonance_indices;

	readindata::load_resonance_info(all_particles, chosen_resonance_indices, TFO, particle_idx);
	n_resonances = chosen_resonance_indices.size();

	/*
	thermal_resonance_spectra_re.resize(n_k_pts*n_resonances*n_pT_pts*n_pphi_pts*n_pY_pts);
	thermal_resonance_spectra_im.resize(n_k_pts*n_resonances*n_pT_pts*n_pphi_pts*n_pY_pts);

	//start some actual calculations here
	//load thermal spectra for all needed resonances
	for (int ik = 0; ik < n_k_pts; ik++)
	for (int iRes = 0; iRes < n_resonances; iRes++)
	for (int ipT = 0; ipT < n_pT_pts; ipT++)
	for (int ipphi = 0; ipphi < n_pphi_pts; ipphi++)
	for (int ipY = 0; ipY < n_pY_pts; ipY++)
	{
		complex<double> z = thermal_FT_spectra( &(all_particles[iRes]xxxxxxxxxxxx),
									k_pts[ik], pT_pts[ipT],
									pphi_pts[ipphi], Del_pY_pts[ipY] );
		thermal_resonance_spectra_re[res_vector_indexer(ik, iRes, ipT, ipphi, ipY)]
			= z.real();
		thermal_resonance_spectra_re[res_vector_indexer(ik, iRes, ipT, ipphi, ipY)]
			= z.imag();
	}

	//do not alter thermal_resonance_spectra with resonance feeddown
	full_resonance_spectra_re = thermal_resonance_spectra_re;
	full_resonance_spectra_im = thermal_resonance_spectra_im;

	//get resonance feeddown (spectra updated and returned)
	for (int ik = 0; ik < n_k_pts; ++ik)
	{
		vector<double> * spec_REAL = copy_spectra(&full_resonance_spectra_re, ik);
		vector<double> * spec_IMAG = copy_spectra(&full_resonance_spectra_im, ik);
		
		decay::Compute_phase_space_integrals(all_particles, chosen_resonance_indices, spec_REAL, spec_IMAG, particle_idx);
	}

	//integrate pion "spectra" over d^2pT to get script F_n(k) in notes
	//set_scriptFn_vector();

	//for (int i = 0; i < (int)scriptFn_vector.size(); ++i)
	//	cout << scriptFn_vector[i] << endl;
	*/

	vector<double> full_resonance_spectra_re (n_resonances*n_pT_pts*n_pphi_pts*n_pY_pts);
	vector<double> full_resonance_spectra_im (n_resonances*n_pT_pts*n_pphi_pts*n_pY_pts);

	//for (int iRes = 0; iRes < n_resonances; iRes++)
	//	cout << iRes << "   " << chosen_resonance_indices[iRes] << endl;
	//if (1) exit (1);

	tau_pts = new double [n_tau_pts];
	tau_wts = new double [n_tau_pts];
	x_pts = new double [n_r_pts];
	x_wts = new double [n_r_pts];
	phi_pts = new double [n_phi_pts];
	phi_wts = new double [n_phi_pts];
	gauss_quadrature(n_tau_pts, 1, 0.0, 0.0, 0.0, 10.0, tau_pts, tau_wts);
	gauss_quadrature(n_r_pts, 1, 0.0, 0.0, -1.0, 1.0, x_pts, x_wts);
	gauss_quadrature(n_phi_pts, 1, 0.0, 0.0, 0.0, 2.0*M_PI, phi_pts, phi_wts);


	//set thermal spectra for all needed resonances
	cout << "Starting thermal calculations..." << endl;
	#pragma omp parallel for
	for (int iRes = 0; iRes < n_resonances; iRes++)
	{
		printf("Hello World (call #%d) from thread = %d, nthreads = %d\n",
				iRes, omp_get_thread_num(), omp_get_num_threads());
		/*cout << "Thread = " << omp_get_thread_num()
				<< " of " << omp_get_num_threads() << " currently devoted to "
				<< all_particles[chosen_resonance_indices[iRes]].name
				<< " spectra:\n";
		cout << " * working on "
				<< all_particles[chosen_resonance_indices[iRes]].name
				<< " spectra..." << endl;*/
		//loop over momenta
		for (int ipT = 0; ipT < n_pT_pts; ipT++)
		for (int ipphi = 0; ipphi < n_pphi_pts; ipphi++)
		for (int ipY = 0; ipY < n_pY_pts; ipY++)
		{
			/*cout << "  --> computing "
					<< all_particles[chosen_resonance_indices[iRes]].name
					<< " spectra at pT==" << pT_pts[ipT]
					<< ", pphi==" << pphi_pts[ipphi]
					<< ", pY==" << Del_pY_pts[ipY] << endl;*/
			full_resonance_spectra_re[res_FIX_K_vector_indexer(iRes, ipT, ipphi, ipY)]
				= Cal_dN_dypTdpTdphi_toy_func(
					chosen_resonance_indices[iRes],
					&all_particles,
					pT_pts[ipT], pphi_pts[ipphi], Del_pY_pts[ipY] );
		}
		/*cout << "Finished working on "
				<< all_particles[chosen_resonance_indices[iRes]].name
				<< " spectra!" << endl;*/
	}
	cout << "Finished thermal calculations!" << endl;

	//do the resonance feeddown and update spectra appropriately
	cout << "Starting resonance feeddown..." << endl;
	decay::Compute_phase_space_integrals(
			all_particles,
			chosen_resonance_indices,
			&full_resonance_spectra_re,
			&full_resonance_spectra_im,
			particle_idx );
	cout << "Finished resonance feeddown!" << endl <<
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl << endl;

	//print out results
	for (int iRes = 0; iRes < n_resonances; iRes++)
	for (int ipT = 0; ipT < n_pT_pts; ipT++)
	for (int ipphi = 0; ipphi < n_pphi_pts; ipphi++)
	for (int ipY = 0; ipY < n_pY_pts; ipY++)
		cout << all_particles.at(chosen_resonance_indices.at(iRes)).name << ":   "
				<< pT_pts.at(ipT) << "   " << pphi_pts.at(ipphi) << "   " << Del_pY_pts.at(ipY) << "   "
				<< full_resonance_spectra_re.at(res_FIX_K_vector_indexer(iRes, ipT, ipphi, ipY)) << endl;

	chosen_resonance_indices.clear();

	//cout << "Still having problems?" << endl;

	//clean up
	delete [] tau_pts;
	delete [] tau_wts;
	delete [] x_pts;
	delete [] x_wts;
	delete [] phi_pts;
	delete [] phi_wts;

	return 0;
}

//////////////////////
//End of main
//////////////////////








void set_up_misc()
{
	xi_pts = vector<double>(n_xi_pts);
	xi_wts = vector<double>(n_xi_pts);
	gauss_quadrature(n_xi_pts, 1, 0.0, 0.0, xi_min, xi_max, xi_pts, xi_wts);

	k_pts = vector<double>(n_k_pts);
	k_wts = vector<double>(n_k_pts);
	gauss_quadrature(n_k_pts, 1, 0.0, 0.0, k_min, k_max, k_pts, k_wts);

	pT_pts = vector<double>(n_pT_pts);
	pT_wts = vector<double>(n_pT_pts);
	//gauss_quadrature(n_pT_pts, 1, 0.0, 0.0, pT_min, pT_max, pT_pts, pT_wts);
	gauss_quadrature(n_pT_pts, 5, 0.0, 0.0, 0.0, 13.0*n_pT_pts/15.0, pT_pts, pT_wts);

	pphi_pts = vector<double>(n_pphi_pts);
	pphi_wts = vector<double>(n_pphi_pts);
	gauss_quadrature(n_pphi_pts, 1, 0.0, 0.0, pphi_min, pphi_max, pphi_pts, pphi_wts);

	Del_pY_pts = vector<double>(n_pY_pts);
	//Del_pY_wts = vector<double>(n_pY_pts);
	//gauss_quadrature(n_pY_pts, 1, 0.0, 0.0, Del_pY_min, Del_pY_max, Del_pY_pts, Del_pY_wts);
	linspace(Del_pY_pts, Del_pY_min, Del_pY_max);
}




/*void set_scriptFn_vector()
{
	for (int ik = 0; ik < n_k_pts; ik++)
	for (int ipY = 0; ipY < n_pY_pts; ipY++)
	{
		scriptFn_vector[scriptFn_vector_indexer(ik, ipY)] = 0.0;

		for (int ipT = 0; ipT < n_pT_pts; ipT++)
		for (int ipphi = 0; ipphi < n_pphi_pts; ipphi++)
			scriptFn_vector[scriptFn_vector_indexer(ik, ipY)]
				+= pT_wts[ipT] * pphi_wts[ipphi] * pT_pts[ipT]	//note extra factor of pT for polar integration
					* full_resonance_spectra[res_vector_indexer(chosen_particle_idx, ik, ipT, ipphi, ipY)];
	}

	return;
}*/

double f(double x, double alpha)
{
	return (
			cosh(x) * exp(-alpha * cosh(x) )
		);
}

double get_integral(double kpt, double alpha)
{
	double result = 0.0;

	for (int ixi = 0; ixi < n_xi_pts; ixi++)
		result += xi_wts[ixi] * cos(kpt * xi_pts[ixi]) * f(xi_pts[ixi], alpha);

	return (result);	// == -2 K'_{i kpt}(alpha)
}

complex<double> thermal_FT_spectra(readindata::particle_info * particle, double kpt, double PT, double PPHI, double PY)
{
	double M = particle->mass;
	double d = particle->gspin;
	double Q = particle->charge;
	double MT = sqrt(PT*PT + M*M);
	double prefactor = d*Q*AT*tauFINAL/(8.0*M_PI*M_PI*M_PI*M_PI*TFO*chiQ);
	
	return (
		prefactor*MT*exp(i*kpt*PY)*get_integral(kpt, MT/TFO)
	);
}

vector<double> * copy_spectra(vector<double> * spectra, int ik)
{
	vector<double>::const_iterator first
		= spectra->begin() + res_vector_indexer( ik, 0, 0, 0, 0 );
	vector<double>::const_iterator last
		= spectra->begin() + res_vector_indexer( ik, n_resonances,
													n_pT_pts, n_pphi_pts, n_pY_pts );
	static vector<double> reso_spectra(first, last);

	return ( & reso_spectra );
}




//////////////////////////////////////////////////////////////
//some functions to help check the phase space integrations
//////////////////////////////////////////////////////////////



//define some parameters for the exact emission function
const double Rad = 5.0, Del_tau = 1.0, tau0 = 5.0, etaf = 0.6;

inline double Hfactor(double r, double tau)
{
	return (
			exp( -r*r/(2.0*Rad*Rad) - (tau-tau0)*(tau-tau0)/(2.0*Del_tau*Del_tau) ) / (M_PI*Del_tau)
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
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/hbarC;
	double Tdec = TFO;
	double one_by_Tdec = 1./Tdec;

	double ch_pY = cosh(pY_loc);
	double sh_pY = sinh(pY_loc);

	double result = 0.0;

	for (int ir = 0; ir < n_r_pts; ++ir)
	for (int itau = 0; itau < n_tau_pts; ++itau)
	{
		double tau = tau_pts[itau];

		double mT = sqrt(pT_loc*pT_loc+localmass*localmass);

		//set r-point inside pT loop, since optimal distribution of integration points
		// depends on value of MT
		double rmin = 0.0, rmax = 15.0 * Rad / sqrt(1.0 + mT*one_by_Tdec*etaf*etaf);
		double hw = 0.5 * (rmax - rmin), cen = 0.5 * (rmax + rmin);
		double rpt = cen + hw * x_pts[ir];

		double local_H = Hfactor(rpt, tau);
		double ch_eta_t = cosh(eta_t(rpt));
		double sh_eta_t = sinh(eta_t(rpt));

		double alpha = mT * one_by_Tdec * ch_eta_t;
		double I1 = 2.0 * gsl_sf_bessel_K1(alpha);

		for (int iphi = 0; iphi < n_phi_pts; ++iphi)
		{
			double phipt = phi_pts[iphi];
			double xpt = rpt*cos(phipt), ypt = rpt*sin(phipt);

			if (one_by_Tdec*pT_loc*sh_eta_t*cos(phipt - pphi_loc) > 700.0)
				continue;

			double S_p_with_weight = tau_wts[itau]*hw*x_wts[ir]*phi_wts[iphi]
										*mT*tau*rpt*prefactor*local_H
										*exp( one_by_Tdec*pT_loc*sh_eta_t*cos(phipt - pphi_loc) );

			result += S_p_with_weight * I1;
		}
	}

	return (result);
}




// End of file
