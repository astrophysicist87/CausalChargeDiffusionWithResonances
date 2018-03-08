#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <complex>
#include <algorithm>

#include <gsl/gsl_sf_gamma.h>

using namespace std;

#include "defs.h"

//const int particle_to_study
	// 1 - pion
	// 2 - proton
	// 3 - kaon
int particle_to_study;

bool print_dndn_k = false;
bool print_dndn_Dxi = false;
//white noise is default
bool white_noise = true;
bool white_Green = true;

const double hbarC = 197.33;
const double xi_infinity = 5.0;
const double k_infinity = 2000.0;

const int n_Dy = 51;
double Delta_y_step = 0.1;

const double DQ = 0.162035;	//fm (rough estimate!)
double vQ2, tauQ;

long n_interp;

double vs, Neff, tauf, taui, Tf, Ti, nu, nuVB, ds, sf;
double alpha0, psi0;
double chi_tilde_mu_mu, chi_tilde_T_mu, chi_tilde_T_T, chi_mu_mu, chi_T_mu, chi_T_T, Delta;

double exp_delta, exp_gamma, exp_nu;
double T0, mu0, Tc, Pc, nc, sc, wc, muc;
double A0, A2, A4, C0, B, mui, muf, xi0, xibar0, etaBYs, RD, sPERn, Nf, qD, si, ni;
double a_at_tauf, vs2_at_tauf, vn2_at_tauf, vsigma2_at_tauf;

const int n_xi_pts = 5000;
const int n_k_pts = 500;	//# of k points should be even to avoid poles in 1F1, etc.!!!
const int n_tau_pts = 1001;
double * xi_pts_minf_inf, * xi_wts_minf_inf;
double * k_pts, * k_wts;
double * tau_pts, * tau_wts;
double * T_pts;
double * running_integral_array;

///////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
	//set parameters from command-line input
	// set particles to study...
	particle1.index = atoi(argv[1]);
	particle2.index = atoi(argv[2]);
	Ti = 350.0 / hbarC;
	//Ti = 250.0 / hbarC;

	//set speed of sound and correlation timescale
	vQ2 = atof(argv[3]);
	tauQ = DQ/vQ2;

	//update white noise and white Green's functions boolean variables
	//assume both are the same for now...
	white_noise = (vQ2 > 20.0);
	//white_noise = true;
	white_Green = (vQ2 > 20.0);

	resultsPath = argv[4];

	fraction_of_evolution = (argc > 5) ? atof(argv[5]) : 1.0;

	set_phase_diagram_and_EOS_parameters();

	switch (particle1.index)
	{
		case 1:	//pion
			particle1.mass = 139.57 / hbarC;
			particle1.name = "pi";
			break;
		case 2:	//proton
			particle1.mass = 939.0 / hbarC;
			particle1.name = "p";
			break;
		case 3:	//kaon
			particle1.mass = 493.68 / hbarC;
			particle1.name = "K";
			break;
		default:
			cerr << "Not a supported particle!" << endl;
			exit(1);
			break;
	}
	switch (particle2.index)
	{
		case 1:	//pion
			particle2.mass = 139.57 / hbarC;
			particle2.name = "pi";
			break;
		case 2:	//proton
			particle2.mass = 939.0 / hbarC;
			particle2.name = "p";
			break;
		case 3:	//kaon
			particle2.mass = 493.68 / hbarC;
			particle2.name = "K";
			break;
		default:
			cerr << "Not a supported particle!" << endl;
			exit(1);
			break;
	}

	taui = 0.5;		//fm/c

	si = s_vs_T(Ti);

	//fixing Tf explicitly instead of
	//calculating it from P=0 curve
	Tf = 150.0 / hbarC;
	//Tf = 106.4226 / hbarC;

	sf = s_vs_T(Tf);
	tauf = si * taui / sf;

	// added this to allow "snapshots" throughout evolution
	tauf = taui + fraction_of_evolution * (tauf - taui);
	sf = s_vs_tau(tauf);
	//Tf = compute_newTf(tauf);
	Tf = guess_T(tauf);
	// rest of code runs the same as before

    // initialize other parameters
    ds = 2.0;

	string filename1, filename2, filename3;
	set_outfilenames(filename1, filename2, filename3);

	//set the susceptibilities
	chi_mu_mu = chi_mumu(Tf);
	chi_T_mu = chi_Tmu(Tf);
	chi_T_T = chi_TT(Tf);
	Delta = chi_mu_mu * chi_T_T - chi_T_mu * chi_T_mu;
	chi_tilde_mu_mu = chi_mu_mu / Delta;
	chi_tilde_T_mu = chi_T_mu / Delta;
	chi_tilde_T_T = chi_T_T / Delta;

    // set up grid points for integrations
    xi_pts_minf_inf = new double [n_xi_pts];
    tau_pts = new double [n_tau_pts];
    k_pts = new double [n_k_pts];
    xi_wts_minf_inf = new double [n_xi_pts];
    tau_wts = new double [n_tau_pts];
    k_wts = new double [n_k_pts];

    int tmp = gauss_quadrature(n_xi_pts, 1, 0.0, 0.0, -xi_infinity, xi_infinity, xi_pts_minf_inf, xi_wts_minf_inf);
    //tmp = gauss_quadrature(n_k_pts, 6, 0.0, 0.0, 0.0, 0.0005/vQ2, k_pts, k_wts);
    //tmp = gauss_quadrature(n_k_pts, 5, 0.0, 0.0, 0.0, 0.5/vQ2, k_pts, k_wts);
    tmp = gauss_quadrature(n_k_pts, 1, 0.0, 0.0, 0.0, k_infinity, k_pts, k_wts);
	double inverse_smearing_width = 1.0;
    tmp = gauss_quadrature(n_tau_pts, 1, 0.0, 0.0, taui, tauf, tau_pts, tau_wts);

	T_pts = new double [n_tau_pts];

	//computes tau-dependence of T for remainder of calculation
	for (int it = 0; it < n_tau_pts; ++it)
		T_pts[it] = guess_T(tau_pts[it]);

	/*
	for (int ik = 0; ik < n_k_pts; ++ik)
	for (int it = 0; it < n_tau_pts; ++it)
	{
		complex<double> k = k_pts[ik];
		complex<double> x = tau_pts[it] / tauQ;
		complex<double> lambda = sqrt(0.25 - vQ2*k*k);
		complex<double> mlambda = -sqrt(0.25 - vQ2*k*k);

		complex<double> result1 = Hypergeometric1F1(lambda+0.5, 2.0*lambda+1.0, x);
		complex<double> result2 = Hypergeometric1F1(mlambda+0.5, 2.0*mlambda+1.0, x);
		complex<double> result3 = Hypergeometric1F1(lambda+1.5, 2.0*lambda+1.0, x);
		complex<double> result4 = Hypergeometric1F1(mlambda+1.5, 2.0*mlambda+1.0, x);

		cout << setprecision(15) << (lambda+0.5).real() << "   " << (lambda+0.5).imag() << "   "
				<< (2.0*lambda+1.0).real() << "   " << (2.0*lambda+1.0).imag() << "   "
				<< x.real() << "   " << x.imag() << "   "
				<< result1.real() << "   " << result1.imag() << endl;
		cout << setprecision(15) << (mlambda+0.5).real() << "   " << (mlambda+0.5).imag() << "   "
				<< (2.0*mlambda+1.0).real() << "   " << (2.0*mlambda+1.0).imag() << "   "
				<< x.real() << "   " << x.imag() << "   "
				<< result2.real() << "   " << result2.imag() << endl;
		cout << setprecision(15) << (lambda+1.5).real() << "   " << (lambda+1.5).imag() << "   "
				<< (2.0*lambda+1.0).real() << "   " << (2.0*lambda+1.0).imag() << "   "
				<< x.real() << "   " << x.imag() << "   "
				<< result3.real() << "   " << result3.imag() << endl;
		cout << setprecision(15) << (mlambda+1.5).real() << "   " << (mlambda+1.5).imag() << "   "
				<< (2.0*mlambda+1.0).real() << "   " << (2.0*mlambda+1.0).imag() << "   "
				<< x.real() << "   " << x.imag() << "   "
				<< result4.real() << "   " << result4.imag() << endl;
	}
	if (1) return (0);
	*/

	//get the ensemble averaged spectra
	double norm = integrate_1D(norm_int, xi_pts_minf_inf, xi_wts_minf_inf, n_xi_pts, &particle1);	//by definition of charge balance function (CBF)

	//vectorize calculations to make them run a little faster
	vector<complex<double> > Ftn_particle1_vec, Ftn_particle2_vec;
	vector<complex<double> > Ctnn_vec, Ctnn_no_SC_vec, SC_vec;
	double * SC_arr = new double [n_k_pts];

	for (int ik = 0; ik < n_k_pts; ++ik)
	{
		double k = k_pts[ik];
		current_ik = ik;
		Ftn_particle1_vec.push_back(Ftilde_n(k, &particle1));
		Ftn_particle2_vec.push_back(Ftilde_n(k, &particle2));
	}

	//compute running integral of transport/medium effects
	running_integral_array = new double [n_tau_pts];
	set_running_transport_integral(running_integral_array);

	ofstream output_dndn_k, output_dndn_Dxi;	//files for Fourier output, coordinate-space output
	string dndn_k_filename = filename2;
	string dndn_Dxi_filename = filename3;
	if (print_dndn_k)
		output_dndn_k.open(dndn_k_filename.c_str());
	if (print_dndn_Dxi)
		output_dndn_Dxi.open(dndn_Dxi_filename.c_str());

	double Ctnn_shift = 0.0;
	for (int ik = 0; ik < n_k_pts; ++ik)
	{
		//cout << "ik = " << ik << "..." << endl;
		double k = k_pts[ik];
		current_ik = ik;
		vector<complex<double> > results(3);

		//get density-density correlator
		Ctilde_n_n(k, &results);

		Ctnn_vec.push_back(results[0]);
		Ctnn_no_SC_vec.push_back(results[1]);
		SC_vec.push_back(results[2]);
		SC_arr[ik] = SC_vec[ik].real();

		if (print_dndn_k)
			output_dndn_k << k << "   " << setprecision(15)
							<< results[0].real() << "   "
							<< results[1].real() << "   "
							<< results[2].real() << endl;

		double k_cutoff = 250.0;
		if ( abs(k_pts[ik]) >= k_cutoff )
		{
			//try to get average after curve has asymptoted...
			Ctnn_shift += k_wts[ik] * Ctnn_no_SC_vec[ik].real() / (k_infinity-k_cutoff);
		}
	}

	bool dndn_Dxi_mode = true;
	bool SC_Dxi_mode = !dndn_Dxi_mode;	//do one or the other

	if (dndn_Dxi_mode)
	{
		for (int ixi = 0; ixi < 5001; ++ixi)
		{
			double Delta_xi = -1.0 + (double)ixi * 2.0/5000.0;
			double Ctnn_no_SC_sum = 0.0;
			for (int ik = 0; ik < n_k_pts; ++ik)
			{
				if (ixi == 0)
					cerr << k_pts[ik] << "   " << (Ctnn_no_SC_vec[ik]-Ctnn_shift).real() << "   " << Ctnn_no_SC_vec[ik].real() << "   " << Ctnn_shift << "   " << Ctnn_vec[ik].real() << "   " << SC_vec[ik].real() << endl;
				Ctnn_no_SC_sum += k_wts[ik] * 2.0 * cos(k_pts[ik] * Delta_xi) * (Ctnn_no_SC_vec[ik]-Ctnn_shift).real();
			}
			cout << Delta_xi << "   " << Ctnn_no_SC_sum << endl;
		}
		if (1) exit (0);
	}
	

	if (SC_Dxi_mode)
	{
		//convert self-correlations to Delta_xi space
		const int n_SC_k_pts = 5000;
		double * SC_k_pts = new double [n_SC_k_pts];
		double * SC_k_wts = new double [n_SC_k_pts];
		//tmp = gauss_quadrature(n_SC_k_pts, 6, 0.0, 0.0, 0.0, 0.005/vQ2, SC_k_pts, SC_k_wts);
		//tmp = gauss_quadrature(n_SC_k_pts, 5, 0.0, 0.0, 0.0, 0.5/vQ2, SC_k_pts, SC_k_wts);
		tmp = gauss_quadrature(n_k_pts, 1, 0.0, 0.0, 0.0, k_infinity, SC_k_pts, SC_k_wts);
		for (int ixi = 0; ixi < 5001; ++ixi)
		{
			complex<double> SC_xi_sum = 0.0;
			double Delta_xi = -0.25 + (double)ixi * 0.0001;
			for (int ik = 0; ik < n_SC_k_pts; ++ik)
			{
				double SC_k = SC_k_pts[ik];
				double SC_loc = 0.0;
				if (SC_k < k_pts[0])
				{
					double lSC0 = log(SC_vec[0].real());
					double lSC1 = log(SC_vec[1].real());
					double slope = (lSC1 - lSC0) / ( log(abs(k_pts[1])) - log(abs(k_pts[0])) );
					SC_loc = exp( lSC0 + slope * ( log(abs(SC_k)) - log(abs(k_pts[0])) ) );
					SC_loc = 0.0;	//for now
				}
				else if (SC_k > k_pts[n_k_pts-1])
				{
					double lSC0 = log(SC_vec[n_k_pts-2].real());
					double lSC1 = log(SC_vec[n_k_pts-1].real());
					double slope = (lSC1 - lSC0) / ( log(abs(k_pts[n_k_pts-1])) - log(abs(k_pts[n_k_pts-2])) );
					SC_loc = exp(lSC0 + slope * ( log(abs(SC_k)) - log(abs(k_pts[n_k_pts-2])) ) );
					SC_loc = 0.0;	//for now
				}
				else
				{
					SC_loc = interpolate1D(k_pts, SC_arr, SC_k, n_k_pts, 1, false);
				}
				if (ixi == 0)
					cerr << SC_k << "   " << SC_loc << endl;
				//SC_xi_sum += SC_k_wts[ik] * exp(i * SC_k * Delta_xi) * SC_loc;
				SC_xi_sum += SC_k_wts[ik] * 2.0 * cos(SC_k * Delta_xi) * SC_loc;
			}
			cout << Delta_xi << "   " << SC_xi_sum.real() << endl;
		}
	}

	/*
	//convert self-correlations to Delta_xi space
	double lSC0 = log(SC_vec[0].real());
	double lSC1 = log(SC_vec[1].real());
	double lower_slope = abs((lSC1 - lSC0) / ( log(abs(k_pts[1])) - log(abs(k_pts[0])) ));
	lSC0 = log(SC_vec[n_k_pts-2].real());
	lSC1 = log(SC_vec[n_k_pts-1].real());
	double upper_slope = abs((lSC1 - lSC0) / ( log(abs(k_pts[n_k_pts-1])) - log(abs(k_pts[n_k_pts-2])) ));

	for (int ixi = 0; ixi < 501; ++ixi)
	{
		double Delta_xi = -0.1 + (double)ixi * 0.0004;

		complex<double> lower_region
				= 2.0 * (
									pow(k_infinity, 1.0-lower_slope)
									* Hypergeometric1F2(0.5*(1.0-lower_slope), 0.5, 0.5*(3.0-lower_slope), -0.25*k_infinity*k_infinity*Delta_xi*Delta_xi)
									/ (-1.0+lower_slope)
									+ pow(abs(Delta_xi), -1.0+lower_slope) * gsl_sf_gamma(1.0-lower_slope) * sin(0.5*M_PI*lower_slope)
									);
		complex<double> upper_region
				= 2.0 * (
									pow(k_infinity, 1.0-upper_slope)
									* Hypergeometric1F2(0.5*(1.0-upper_slope), 0.5, 0.5*(3.0-upper_slope), -0.25*k_infinity*k_infinity*Delta_xi*Delta_xi)
									/ (-1.0+upper_slope)
									+ pow(abs(Delta_xi), -1.0+upper_slope) * gsl_sf_gamma(1.0-upper_slope) * sin(0.5*M_PI*upper_slope)
									);

		complex<double> SC_xi_sum = lower_region + upper_region;

		for (int ik = 0; ik < n_k_pts; ++ik)
		{
			double k = k_pts[ik];
			double SC_loc = SC_vec[ik].real();
			SC_xi_sum += k_wts[ik] * exp(i * k * Delta_xi) * SC_loc;
		}
		cout << Delta_xi << "   " << SC_xi_sum.real() << endl;
	}
	*/

	if (print_dndn_k)
		output_dndn_k.close();
	if (print_dndn_Dxi)
		output_dndn_Dxi.close();

	return 0;
}
