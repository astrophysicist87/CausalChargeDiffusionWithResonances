#ifndef MAIN_H
#define MAIN_H

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
#include "thermal.h"
#include "chebyshev.h"
#include "Stopwatch.h"

extern const int operation_mode;

//space-time integration grid
const int n_r_pts = 51;
double * x_pts, * x_wts;

vector<double> xi_pts, xi_wts;
vector<double> k_pts, k_wts;

int n_resonances = -1;
const double AT = 1.0, tauFINAL = 10.0, TFO = 0.15444512;	//GeV
const double chiQ = 2.0/3.0;

vector<complex<double> > scriptFn_vector;
vector<double> thermal_resonance_spectra_re, thermal_resonance_spectra_im;
vector<double> full_resonance_spectra_re, full_resonance_spectra_im;
vector<double> tmp_full_resonance_spectra_re;
vector<double> tmp_full_resonance_spectra_im;

vector<double> pT_pts, pT_wts;
vector<double> pphi_pts, pphi_wts;
vector<double> Del_pY_pts, Del_pY_wts;
vector<double> tmp_Del_pY_pts, tmp_Del_pY_wts;	//use these and Chebyshev to construct actual Del_pY points...

//function prototypes
void operation_mode_0();
void operation_mode_1();

//
double g(double x, double width);
void set_up_misc();
void set_scriptFn_vector();
double f(double x, double alpha);
double get_integral(double kpt, double alpha);
complex<double> thermal_FT_spectra(readindata::particle_info * particle, double kpt, double PT, double PPHI, double PY);
vector<double> * copy_spectra(vector<double> * spectra, int ik);
double Cal_dN_dypTdpTdphi_toy_func(
		int local_pid, vector<readindata::particle_info> * all_particles,
		double pT_loc, double pphi_loc, double pY_loc );
void set_thermal_spectra( vector<readindata::particle_info> all_particles,
							vector<int> chosen_resonance_indices );
void do_chebyshev_interpolation();
void print_results(vector<readindata::particle_info> * all_particles_ptr, vector<int> * chosen_resonance_indices_ptr);


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

void operation_mode_0()
{
	Stopwatch sw;
	sw.Start();
	//initializes a few things
	set_up_misc();
	//gsl_set_error_handler_off();

	//load resonance information
	const int particle_idx = 1;		// pi^+
	vector<readindata::particle_info> all_particles;
	vector<int> chosen_resonance_indices;

	readindata::load_resonance_info(all_particles, chosen_resonance_indices, TFO, particle_idx);
	n_resonances = chosen_resonance_indices.size();

	tmp_full_resonance_spectra_re.resize(n_resonances*n_pT_pts*n_pphi_pts*tmp_n_pY_pts);
	tmp_full_resonance_spectra_im.resize(n_resonances*n_pT_pts*n_pphi_pts*tmp_n_pY_pts);

	sw.Stop();
	cout << "Set up calculation in " << sw.printTime() << " sec." << endl;
	sw.Reset();

	//one of the main functions
	sw.Start();
	set_thermal_spectra(all_particles, chosen_resonance_indices);
	sw.Stop();
	cout << "Got thermal spectra in " << sw.printTime() << " sec." << endl;
	sw.Reset();

	//print out results
	for (int iRes = 0; iRes < n_resonances; iRes++)
	for (int ipT = 0; ipT < n_pT_pts; ipT++)
	for (int ipphi = 0; ipphi < n_pphi_pts; ipphi++)
	for (int ipY = 0; ipY < tmp_n_pY_pts; ipY++)
		cout << all_particles.at(chosen_resonance_indices.at(iRes)).name << ":   "
				<< pT_pts.at(ipT) << "   " << pphi_pts.at(ipphi) << "   " << tmp_Del_pY_pts.at(ipY) << "   "
				<< tmp_full_resonance_spectra_re.at(res_FIX_K_vector_indexer(iRes, ipT, ipphi, ipY)) << endl;

//if (1) exit(8);

	full_resonance_spectra_re.resize(n_resonances*n_pT_pts*n_pphi_pts*n_pY_pts);
	full_resonance_spectra_im.resize(n_resonances*n_pT_pts*n_pphi_pts*n_pY_pts);

	///////////////
	//one of the main functions
	sw.Start();
	do_chebyshev_interpolation();
	sw.Stop();
	cout << "Did Chebyshev interpolation in " << sw.printTime() << " sec." << endl;
	sw.Reset();
	///////////////

	///////////////
	//one of the main functions
	sw.Start();
	//do the resonance feeddown and update spectra appropriately
	cout << "Starting resonance feeddown..." << endl;
	decay::Compute_phase_space_integrals(
			all_particles,
			chosen_resonance_indices,
			&full_resonance_spectra_re,
			&full_resonance_spectra_im,
			particle_idx );
	sw.Stop();
	cout << "Did resonance feeddown in " << sw.printTime() << " sec." << endl;
	sw.Reset();
	cout << "Finished resonance feeddown!" << endl <<
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl << endl;
	///////////////

	///////////////
	print_results(&all_particles, &chosen_resonance_indices);

	chosen_resonance_indices.clear();

	//clean up
	delete [] x_pts;
	delete [] x_wts;

	return;
}


void operation_mode_1()
{

	//initializes a few things
	set_up_misc();
	//gsl_set_error_handler_off();

	//load resonance information
	const int particle_idx = 1;		// pi^+
	vector<readindata::particle_info> all_particles;
	vector<int> chosen_resonance_indices;

	readindata::load_resonance_info(all_particles, chosen_resonance_indices, TFO, particle_idx);
	n_resonances = chosen_resonance_indices.size();

	
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
		complex<double> z
			= thermal_FT_spectra(
				&(all_particles[chosen_resonance_indices[iRes]]),
				k_pts[ik], pT_pts[ipT],
				pphi_pts[ipphi], Del_pY_pts[ipY] );
		thermal_resonance_spectra_re
			[res_vector_indexer(ik, iRes, ipT, ipphi, ipY)]
			= z.real();
		thermal_resonance_spectra_im
			[res_vector_indexer(ik, iRes, ipT, ipphi, ipY)]
			= z.imag();

		cout << all_particles.at
				(chosen_resonance_indices.at(iRes)).name << ":   "
				<< k_pts.at(ik) << "   "
				<< pT_pts.at(ipT) << "   "
				<< pphi_pts.at(ipphi) << "   "
				<< Del_pY_pts.at(ipY) << "   "
				<< thermal_resonance_spectra_re
					[res_vector_indexer(ik, iRes, ipT, ipphi, ipY)] << "   "
				<< thermal_resonance_spectra_im
					[res_vector_indexer(ik, iRes, ipT, ipphi, ipY)]
				<< endl;
	}

if (1) exit (8);

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
	

	tmp_full_resonance_spectra_re.resize(n_resonances*n_pT_pts*n_pphi_pts*tmp_n_pY_pts);
	tmp_full_resonance_spectra_im.resize(n_resonances*n_pT_pts*n_pphi_pts*tmp_n_pY_pts);

	full_resonance_spectra_re.resize(n_resonances*n_pT_pts*n_pphi_pts*n_pY_pts);
	full_resonance_spectra_im.resize(n_resonances*n_pT_pts*n_pphi_pts*n_pY_pts);

	//one of the main functions
	do_chebyshev_interpolation();

	//one of the main functions
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

	//clean up
	delete [] x_pts;
	delete [] x_wts;


	return;
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
double g(double x, double width)
{
	return (sin(x)*exp(-x*width));
}

void set_up_misc()
{
	xi_pts = vector<double>(n_xi_pts);
	xi_wts = vector<double>(n_xi_pts);
	gauss_quadrature(n_xi_pts, 1, 0.0, 0.0, xi_min, xi_max, xi_pts, xi_wts);

	k_pts = vector<double>(n_k_pts);
	k_wts = vector<double>(n_k_pts);
	gauss_quadrature(n_k_pts, 1, 0.0, 0.0, k_min, k_max, k_pts, k_wts);

	x_pts = new double [n_r_pts];
	x_wts = new double [n_r_pts];
	gauss_quadrature(n_r_pts, 1, 0.0, 0.0, -1.0, 1.0, x_pts, x_wts);

	pT_pts = vector<double>(n_pT_pts);
	pT_wts = vector<double>(n_pT_pts);
	//gauss_quadrature(n_pT_pts, 1, 0.0, 0.0, pT_min, pT_max, pT_pts, pT_wts);
	gauss_quadrature(n_pT_pts, 5, 0.0, 0.0, 0.0, 13.0*n_pT_pts/15.0, pT_pts, pT_wts);

	pphi_pts = vector<double>(n_pphi_pts);
	pphi_wts = vector<double>(n_pphi_pts);
	gauss_quadrature(n_pphi_pts, 1, 0.0, 0.0, pphi_min, pphi_max, pphi_pts, pphi_wts);

	tmp_Del_pY_pts = vector<double>(tmp_n_pY_pts);
	tmp_Del_pY_wts = vector<double>(tmp_n_pY_pts);
	Del_pY_pts = vector<double>(n_pY_pts);
	gauss_quadrature(tmp_n_pY_pts, 1, 0.0, 0.0, Del_pY_min, Del_pY_max, tmp_Del_pY_pts, tmp_Del_pY_wts);
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
													n_pT_pts, n_pphi_pts, tmp_n_pY_pts );
	static vector<double> reso_spectra(first, last);

	return ( & reso_spectra );
}


/////////////////////////////////////////////////////
//main routines called in driver function
/////////////////////////////////////////////////////

void set_thermal_spectra( vector<readindata::particle_info> all_particles,
							vector<int> chosen_resonance_indices )
{
	cout << "Starting thermal calculations..." << endl;
	for (int iRes = 0; iRes < n_resonances; iRes++)
	{
		cout << " * working on "
				<< all_particles[chosen_resonance_indices[iRes]].name
				<< " spectra..." << endl;

		for (int ipT = 0; ipT < n_pT_pts; ipT++)
		for (int ipphi = 0; ipphi < n_pphi_pts; ipphi++)
		for (int ipY = 0; ipY < tmp_n_pY_pts; ipY++)
		{
			tmp_full_resonance_spectra_re[res_FIX_K_vector_indexer(iRes, ipT, ipphi, ipY)]
				= Cal_dN_dypTdpTdphi_toy_func(
					chosen_resonance_indices[iRes],
					&all_particles,
					pT_pts[ipT], pphi_pts[ipphi], tmp_Del_pY_pts[ipY] );
		}
		cout << "Finished working on "
				<< all_particles[chosen_resonance_indices[iRes]].name
				<< " spectra!" << endl;
	}
	cout << "Finished thermal calculations!" << endl;

	return;
}


void do_chebyshev_interpolation()
{
	vector<double> old_vals( tmp_n_pY_pts );
	vector<double> new_vals( n_pY_pts );

	cheb_int::set_up( tmp_n_pY_pts );

	Stopwatch sw_DCI;

	for (int iRes = 0; iRes < n_resonances; iRes++)
	for (int ipT = 0; ipT < n_pT_pts; ipT++)
	for (int ipphi = 0; ipphi < n_pphi_pts; ipphi++)
	{
		sw_DCI.Start();
		for (int ipY = 0; ipY < tmp_n_pY_pts; ++ipY)
		{
			old_vals.at(ipY) = tmp_full_resonance_spectra_re[res_FIX_K_vector_indexer(iRes, ipT, ipphi, ipY)];
			//imag[ipY] = tmp_full_resonance_spectra_im[res_FIX_K_vector_indexer(iRes, ipT, ipphi, ipY)];
			//cout << "old: " << iRes << "   " << ipT << "   " << ipphi << "   " << ipY << "   " << old_vals.at(ipY) << endl;
		}

		cheb_int::chebyshev_interpolate(&tmp_Del_pY_pts, &old_vals, &Del_pY_pts, &new_vals);
		//cheb_int::chebyshev_interpolate(&tmp_Del_pY_pts, &imag, &Del_pY_pts, &interp_imag);

		for (int ipY = 0; ipY < n_pY_pts; ++ipY)
		{
			full_resonance_spectra_re[( ( iRes * n_pT_pts + ipT ) * n_pphi_pts + ipphi ) * n_pY_pts + ipY]
				 = new_vals.at(ipY);
			//full_resonance_spectra_im[res_FIX_K_vector_indexer(iRes, ipT, ipphi, ipY)]
			//	 = interp_imag[ipY];
			//cout << "new: " << iRes << "   " << ipT << "   " << ipphi << "   " << ipY << "   " << new_vals.at(ipY) << endl;			
		}

		sw_DCI.Stop();
		cout << "cheb_int::Finished (" << iRes << ", " << ipT << ", " << ipphi << ") loop in " << sw_DCI.printTime() << " sec." << endl;
		sw_DCI.Reset();

//if (1) exit (8);
	}

	cheb_int::clean_up();

	return;
}



void print_results(vector<readindata::particle_info> * all_particles_ptr, vector<int> * chosen_resonance_indices_ptr)
{
	//print out results
	for (int iRes = 0; iRes < n_resonances; iRes++)
	for (int ipT = 0; ipT < n_pT_pts; ipT++)
	for (int ipphi = 0; ipphi < n_pphi_pts; ipphi++)
	for (int ipY = 0; ipY < n_pY_pts; ipY++)
		cout << all_particles_ptr->at(chosen_resonance_indices_ptr->at(iRes)).name << ":   "
				<< pT_pts.at(ipT) << "   " << pphi_pts.at(ipphi) << "   " << Del_pY_pts.at(ipY) << "   "
				<< full_resonance_spectra_re.at(res_FIX_K_vector_indexer(iRes, ipT, ipphi, ipY)) << endl;
	return;
}


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



// End of file

#endif
