#ifndef DEFS_H
#define DEFS_H

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
#include "gauss/gauss_quadrature.h"
#include "lib.h"
#include "asymptotics.h"

/*USAGE:
debugger(__LINE__, __FILE__);
*/
void inline debugger(int cln, const char* cfn)
{
	cerr << "You made it to " << cfn << ":" << cln << "!" << endl;
	return;
}

inline int res_vector_indexer(const int ik, const int iRes, const int ipT, const int ipphi, const int ipY)
{
	return (
		( ( ( ik * n_resonances + iRes ) * n_pT_pts + ipT ) * n_pphi_pts + ipphi ) * tmp_n_pY_pts + ipY
	);
}

inline int res_FIX_K_vector_indexer(const int iRes, const int ipT, const int ipphi, const int ipY)
{
	return (
		( ( iRes * n_pT_pts + ipT ) * n_pphi_pts + ipphi ) * tmp_n_pY_pts + ipY
	);
}

inline int scriptFn_vector_indexer(const int ik, const int ipY)
{
	return (
		ik * tmp_n_pY_pts + ipY
	);
}

inline int mom_indexer(const int ipT, const int ipphi, const int ipY)
{
	return (
		( ipT * n_pphi_pts + ipphi ) * tmp_n_pY_pts + ipY
	);
}

/*
string truestring = "true";
string falsestring = "false";

extern bool white_noise;
extern bool white_Green;
string resultsPath;

inline string return_boolean_string(bool test){return (test ? truestring : falsestring);}

struct chosen_particle
{
	int index;
	double mass;
	string name;
};

extern const double hbarC;

extern const double DQ;
extern double tauQ, vQ2;

double fraction_of_evolution;

extern long n_interp;

extern double vs, Neff, taui, tauf, Ti, Tf;
extern double chi_tilde_mu_mu, chi_tilde_T_mu, chi_tilde_T_T, Delta;

extern double chi_T_T, chi_T_mu, chi_mu_mu;

extern const int n_xi_pts;
extern const int n_k_pts;
extern const int n_tau_pts;
extern double * xi_pts_minf_inf, * xi_wts_minf_inf;
extern double * k_pts, * k_wts;
extern double * tau_pts, * tau_wts;
extern double * T_pts;
extern double * running_integral_array;

chosen_particle particle1;
chosen_particle particle2;

int current_ik = -1;

extern double T0, mu0, Tc, Pc, nc, sc, wc, muc;
extern double A0, A2, A4, C0, B, xi0, xibar0, etaBYs, RD, sPERn, Nf, qD, si;

//general functions
inline double incompleteGamma3(double x)
{
	return (exp(-x)*(2.+x*(2.+x)));
}

inline double incompleteGamma4(double x)
{
	return (exp(-x)*(6. + x*(6. + x*( 3. + x))));
}

inline double incompleteGamma5(double x)
{
	return (exp(-x)*(24.+x*(24.+x*(12.+x*(4.+x)))));
}

//equation of state and other thermodynamic relations
inline double P(double T)
{
	return(
		A4*pow(T, 4.0) - C0*T*T - B
	);
}

inline void set_phase_diagram_and_EOS_parameters()
{
	Nf = 3.0;											//number of massless flavors
	T0 = 170.0 / hbarC;											//phase transition curve, T scale
	mu0 = 1218.48 / hbarC;										//phase transition curve, mu scale
	A4 = M_PI*M_PI*(16.0 + 10.5*Nf) / 90.0;				//coeff in P(T) (i.e., EOS)
	A2 = Nf / 18.0;										//same
	A0 = Nf / (324.0 * M_PI * M_PI);					//same
	//C0 = mu0*mu0*( A2 - 2.0*A0*mu0*mu0 / (T0*T0) );		//same
	C0 = 0.0;											//same
	B = 0.8 * pow(T0, 4.0);								//same

	return;
}

// functions to guess seed value for T,
// followed by functions to iteratively solve equations
// for T

inline double guess_T(double tau)
{
	return (Ti * pow(taui / tau, 1.0/3.0));
}

///////////////////////////////////////////////////////////
//two separate definitions of s, for convenience
///////////////////////////////////////////////////////////
inline double s_vs_tau(double tau)
{
	return (si * taui / tau);
}

inline double s_vs_T(double T)
{
	return (
		-2.0 * C0 * T + 4.0 * A4 * T*T*T
	);
}

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

inline double w(double T)
{
	return (T * s_vs_T(T));
}

inline double chi_TT(double T)
{
	return ( 2.0 * ( -C0 + 6.0 * A4 * T * T ) );
}

inline double chi_Tmu(double T)
{
	T;
	return ( 0.0 );
}

inline double chi_mumu(double T)
{
	return ( 4.0 * A2 * T * T );
}

inline double norm_int(double x, void * p)
{
	struct chosen_particle * params = (struct chosen_particle *)p;
	double cx = cosh(x);
	double mByT = (params->mass) / Tf;
	return (incompleteGamma3(mByT * cx) / (cx*cx));
}

inline double Fn(double x, void * p)
{
	struct chosen_particle * params = (struct chosen_particle *)p;
	double cx = cosh(x);
	
	double c2 = chi_tilde_T_T;

	double mByT = (params->mass) / Tf;

	return ( c2 * incompleteGamma3(mByT * cx) / (cx*cx) );
}

inline complex<double> Ftilde_n(double k, void * p)
{
	return (integrate_1D_FT(Fn, xi_pts_minf_inf, xi_wts_minf_inf, n_xi_pts, k, p));
}











inline complex<double> psi_plus(complex<double> k, complex<double> x)
{
	complex<double> one_fourth = 0.25;
	complex<double> lambda = sqrt(one_fourth - vQ2*k*k);
	return (
			pow(x,lambda-0.5)
			* exp(-x) * Hypergeometric1F1(lambda+1.5, 2.0*lambda+1.0, x)
			);
}

inline complex<double> psi_minus(complex<double> k, complex<double> x)
{
	complex<double> one_fourth = 0.25;
	complex<double> mlambda = -sqrt(one_fourth - vQ2*k*k);
	return (
			pow(x,mlambda-0.5)
			* exp(-x) * Hypergeometric1F1(mlambda+1.5, 2.0*mlambda+1.0, x)
			);
}

inline complex<double> psi_dot_plus(complex<double> k, complex<double> x)
{
	complex<double> one_fourth = 0.25;
	complex<double> lambda = sqrt(one_fourth - vQ2*k*k);
	return (
			0.5 * (2.0*lambda-1.0) * pow(x,lambda-1.5)
			* exp(-x) * Hypergeometric1F1(lambda+0.5, 2.0*lambda+1.0, x)
			);
}

inline complex<double> psi_dot_minus(complex<double> k, complex<double> x)
{
	complex<double> one_fourth = 0.25;
	complex<double> mlambda = -sqrt(one_fourth - vQ2*k*k);
	complex<double> tmp = Hypergeometric1F1(mlambda+0.5, 2.0*mlambda+1.0, x);
	return (
			0.5 * (2.0*mlambda-1.0) * pow(x,mlambda-1.5)
			* exp(-x) * Hypergeometric1F1(mlambda+0.5, 2.0*mlambda+1.0, x)
			);
}

inline double mediumN(double tau_loc)
{
	//double T_loc = interpolate1D(tau_pts, T_pts, tau_loc, n_tau_pts, 0, false);
	double T_loc = guess_T(tau_loc);
	double s_loc = s_vs_T(T_loc);
	double chi_Q = chi_mumu(T_loc);
	double numerator = 2.0 * DQ * chi_Q * T_loc;
	double denominator = tau_loc * s_loc * s_loc;	//assume transverse area A already accounted for...
	return ( numerator / denominator );
}

inline complex<double> Gtilde_n_white(double k, double tau, double taup)
{
	double arg = DQ * k * k * ((1.0/tau) - (1.0/taup));	//since tau==tauf always in this calculation,
														//exp(arg) will always be Gaussian in k
	return ( i * k * exp(arg) ); //dimensionless
}



inline complex<double> asymptotic_Gtilde_n_color(complex<double> k, double tau, double taup)
{
	double x = tau / tauQ;
	double xp = taup / tauQ;
	complex<double> one_fourth = 0.25;
	complex<double> lambda = sqrt(one_fourth - vQ2*k*k);

	complex<double> prefactor = 0.25*M_PI*sqrt( taup / tau )
						* exp( 0.5 * (xp - x) )
						* csc(M_PI * lambda);
	double vQ = sqrt(vQ2);

	complex<double> I_lambda_x_by_2 = asymptotics::I(lambda, 0.5*x);
	complex<double> I_prime_lambda_x_by_2 = asymptotics::Iprime(lambda, 0.5*x);
	complex<double> I_mlambda_x_by_2 = asymptotics::I(-lambda, 0.5*x);
	complex<double> I_prime_mlambda_x_by_2 = asymptotics::Iprime(-lambda, 0.5*x);

	complex<double> I_lambda_xp_by_2 = asymptotics::I(lambda, 0.5*xp);
	complex<double> I_mlambda_xp_by_2 = asymptotics::I(-lambda, 0.5*xp);

	complex<double> term1 = ( (x + 1.0) * I_lambda_x_by_2 + x * I_prime_lambda_x_by_2 )
							* I_mlambda_xp_by_2;
	complex<double> term2 = ( (x + 1.0) * I_mlambda_x_by_2 + x * I_prime_mlambda_x_by_2 )
							* I_lambda_xp_by_2;

	complex<double> mainfactor = term1 - term2;

	return ( i*k*prefactor*mainfactor );
}




inline complex<double> Gtilde_n_color(double k, double tau, double taup)
{
	double xp = taup / tauQ;
	complex<double> one_fourth = 0.25;
	complex<double> lambda = sqrt(one_fourth - vQ2*k*k);

	complex<double> psi_plus_at_tau = psi_plus(k, tau / tauQ);
	complex<double> psi_minus_at_tau = psi_minus(k, tau / tauQ);
	complex<double> psi_plus_at_taup = psi_plus(k, taup / tauQ);
	complex<double> psi_minus_at_taup = psi_minus(k, taup / tauQ);
	complex<double> psi_dot_plus_at_taup = psi_dot_plus(k, taup / tauQ);
	complex<double> psi_dot_minus_at_taup = psi_dot_minus(k, taup / tauQ);

	complex<double> n1 = psi_plus_at_tau * psi_dot_minus_at_taup;
	complex<double> n2 = psi_minus_at_tau * psi_dot_plus_at_taup;
	complex<double> numerator = n1 - n2;
	complex<double> denominator = -2.0 * lambda * exp(-xp) / (xp*xp);

	complex<double> result = numerator / denominator;

	//if the numerator and denominator are too small, assume loss of significance
	// --> just use the white noise Green's function instead
	if ( abs(denominator) < 1.e-15 && abs(numerator) < 1.e-15
			&& 2.0*abs(n1-n2)/( abs(n1)+abs(n2) ) < 1.e-10 )
	{
		if (vQ2*k*k < 1.0)	//play with exactly which point to turn on pure white noise approximation
		{
			result = exp(DQ * k * k * ((1.0/tau) - (1.0/taup)) );
		}
		else
		{
			return ( asymptotic_Gtilde_n_color(k, tau, taup) );
		}
		//cerr << "Warning: " << k << "   " << abs(denominator) << "   " << abs(numerator)
		//		<< "   " << 2.0*abs(n1-n2)/( abs(n1)+abs(n2) ) << endl;
	}

	return ( i * k * result ); //dimensionless
}



inline void tau_integration_WhiteGreen(
					complex<double> (*Gtilde_X)(double, double, double),
					complex<double> (*Gtilde_Y)(double, double, double),
					double k, vector<complex<double> > * results)
{
	complex<double> result(0,0);
	double * local_x_pts = new double [n_tau_pts];
	double * local_x_wts = new double [n_tau_pts];
	gauss_quadrature(n_tau_pts, 1, 0.0, 0.0, -1.0, 1.0, local_x_pts, local_x_wts);
	double hw_loc = 0.5 * (tauf - taui);
	double cen_loc = 0.5 * (tauf + taui);
	
	for (int it = 0; it < n_tau_pts; ++it)
	{
		//double tau_loc = tau_pts[it];
		double tau_loc = cen_loc + hw_loc * local_x_pts[it];
		//double T_loc = T_pts[it];
		double T_loc = guess_T(tau_loc);
		double s_loc = s_vs_T(T_loc);
		double chi_Q = chi_mumu(T_loc);
		complex<double> tmp_result = hw_loc * local_x_wts[it] * ( 2.0 * DQ * chi_Q * T_loc / tau_loc )
										* (*Gtilde_X)(k, tauf, tau_loc) * (*Gtilde_Y)(-k, tauf, tau_loc);
		result += tmp_result;
	}
	
	double local_Tf = Tf;

	//need two different tau_integration routines because of this!!!
	double self_correlation = chi_mumu(local_Tf) * local_Tf * tauf;	//allows evaluating at other times

	delete [] local_x_pts;
	delete [] local_x_wts;
//cerr << "self_correlation = " << k << "   " << self_correlation << endl;
//cout << "Sanity check: " << k << "   " << self_correlation << "   " << result.real() << endl;
//if (1) return (0);

	(*results)[0] = result;
	(*results)[1] = self_correlation - result;
	(*results)[2] = self_correlation;

	//cout << k << "   " << result.real() << "   " << self_correlation - result.real() << "   " << self_correlation << endl;

	return;
}

inline void tau_integration_coloredGreen(
					complex<double> (*Gtilde_X)(double, double, double),
					complex<double> (*Gtilde_Y)(double, double, double),
					double k, vector<complex<double> > * results)
{
	complex<double> result(0,0), SC_result(0,0);
	double * local_x_pts = new double [n_tau_pts];
	double * local_x_wts = new double [n_tau_pts];
	gauss_quadrature(n_tau_pts, 1, 0.0, 0.0, -1.0, 1.0, local_x_pts, local_x_wts);
	double hw_loc = 0.5 * (tauf - taui);
	double cen_loc = 0.5 * (tauf + taui);

	double self_correlation = chi_mumu(Tf) * Tf * tauf;

	for (int it = 0; it < n_tau_pts; ++it)
	{
		double tau_loc = cen_loc + hw_loc * local_x_pts[it];
		double T_loc = guess_T(tau_loc);
		double chi_Q = chi_mumu(T_loc);

		complex<double> tmp_result = hw_loc * local_x_wts[it] * ( 2.0 * DQ * chi_Q * T_loc / tau_loc ) *
										(*Gtilde_X)(k, tauf, tau_loc) * (*Gtilde_Y)(-k, tauf, tau_loc);
		result += tmp_result;
		//SC_result += hw_loc * local_x_wts[it] * ( 2.0 * DQ * chi_Q * T_loc / tau_loc ) * GG_self_correlations(k, tau_loc, tau_loc);
	}

	delete [] local_x_pts;
	delete [] local_x_wts;

	(*results)[0] = result;
	(*results)[1] = self_correlation - result;
	(*results)[2] = self_correlation;
	//(*results)[1] = SC_result - result;
	//(*results)[2] = SC_result;

	return;
}

inline void set_running_transport_integral(double * run_int_array)
{
	const int n_x_pts = 11;
	double * x_pts = new double [n_x_pts];
	double * x_wts = new double [n_x_pts];
	gauss_quadrature(n_x_pts, 1, 0.0, 0.0, -1.0, 1.0, x_pts, x_wts);

	run_int_array[0] = 0.0;
	for (int it = 0; it < n_tau_pts-1; ++it)
	{
		double sum = 0.0;
		double t0 = tau_pts[it], t1 = tau_pts[it+1];
		double cen = 0.5*(t0+t1);
		double hw = 0.5*(t1-t0);
		for (int ix = 0; ix < n_x_pts; ++ix)
		{
			double t_loc = cen + hw * x_pts[ix];
			sum += x_wts[ix] * hw * exp(t_loc / tauQ) * mediumN(t_loc);
			//cout << "Medium check: " << t_loc << "   " << mediumN(t_loc) << endl;
		}
		run_int_array[it+1] = exp(-t1 / tauQ) * ( tauQ * exp(t0 / tauQ) * run_int_array[it] + sum ) / tauQ;	//this array contains eta(x)
		//cout << t0 << "   " << run_int_array[it] << endl;
	}
	delete [] x_pts;
	delete [] x_wts;

	return;
}

inline void colored_tau_integration(
					complex<double> (*Gtilde_X)(double, double, double),
					complex<double> (*Gtilde_Y)(double, double, double),
					double k, vector<complex<double> > * results)
{
	complex<double> locsum = 0.0, SC_sum = 0.0;

	double * local_x_pts = new double [n_tau_pts];
	double * local_x_wts = new double [n_tau_pts];
	gauss_quadrature(n_tau_pts, 1, 0.0, 0.0, -1.0, 1.0, local_x_pts, local_x_wts);
	double new_hw_loc = 0.5 * (tauf - taui);
	double new_cen_loc = 0.5 * (tauf + taui);

	complex<double> delta_n_delta_n_self_corr = 0.0;
	for (int itp = 0; itp < n_tau_pts; ++itp)
	{
		double tau_p_local = new_cen_loc + new_hw_loc * local_x_pts[itp];
		delta_n_delta_n_self_corr += new_hw_loc * local_x_wts[itp] / tau_p_local
										* ( exp(-(tauf - tau_p_local)/tauQ)
											- exp(-(tauf - taui + tau_p_local)/tauQ) )
										* (*Gtilde_X)(-k, tauf, tau_p_local);
	}
	delta_n_delta_n_self_corr *= -tauf*tauf*chi_mumu(Tf)*Tf/(tauQ*i*k);

	const int n_x_pts = 201;	//try this
	double * x_pts = new double [n_x_pts];
	double * x_wts = new double [n_x_pts];
	gauss_quadrature(n_x_pts, 1, 0.0, 0.0, -1.0, 1.0, x_pts, x_wts);

	//this bounds the interval where the integrand is large-ish
	double delta_tau_lower = -20.0 * tauQ, delta_tau_upper = 20.0 * tauQ;

	for (int itp = 0; itp < n_tau_pts; ++itp)
	{
		double tX_loc = new_cen_loc + new_hw_loc * local_x_pts[itp];
		double TX_loc = guess_T(tX_loc);

		double sX_loc = s_vs_T(TX_loc);

		complex<double> factor_X = sX_loc * (*Gtilde_X)(k, tauf, tX_loc);		//extra factor of entropy!!!

		//if lower limit goes before beginning of lifetime, just start at tau0
		double tau_lower = max(taui, tX_loc + delta_tau_lower);
		//if upper limit goes past end of lifetime, just end at tauf
		double tau_upper = min(tauf, tX_loc + delta_tau_upper);

		double hw_loc = 0.5 * (tau_upper - tau_lower);
		double cen_loc = 0.5 * (tau_upper + tau_lower);

		complex<double> sum_X = 0.0, sum_SC_X = 0.0;
		for (int ix = 0; ix < n_x_pts; ++ix)
		{
			double tY_loc = cen_loc + hw_loc * x_pts[ix];
			double TY_loc = guess_T(tY_loc);

			double sY_loc = s_vs_T(TY_loc);

			complex<double> factor_Y = sY_loc * (*Gtilde_Y)(-k, tauf, tY_loc);	//extra factor of entropy!!!

			double min_tp_tpp = min(tX_loc, tY_loc);
			double eta_at_min_tp_tpp = interpolate1D(tau_pts, running_integral_array, min_tp_tpp, n_tau_pts, 0, false, true);

			double sum_XY = exp(-abs(tX_loc - tY_loc) / tauQ) * eta_at_min_tp_tpp / (2.0*tauQ);
			sum_X += hw_loc * x_wts[ix] * factor_X * factor_Y * sum_XY;
			//sum_SC_X += hw_loc * x_wts[ix] * sX_loc * sY_loc * GG_self_correlations(k, tX_loc, tY_loc) * sum_XY;
		}
		locsum += new_hw_loc * local_x_wts[itp] * sum_X;
		//SC_sum += new_hw_loc * local_x_wts[itp] * sum_SC_X;
	}

	delete [] x_pts;
	delete [] x_wts;
	delete [] local_x_pts;
	delete [] local_x_wts;

	(*results)[0] = locsum;
	(*results)[1] = delta_n_delta_n_self_corr - locsum;
	(*results)[2] = delta_n_delta_n_self_corr;
	//(*results)[1] = SC_sum - locsum;
	//(*results)[2] = SC_sum;

	//cout << k << "   " << locsum.real() << "   " << delta_n_delta_n_self_corr.real() - locsum.real() << "   " << delta_n_delta_n_self_corr.real() << "   "
	//		<< SC_sum.real() - locsum.real() << "   " << SC_sum.real() << endl;

	return;
}

inline void Ctilde_n_n(double k, vector<complex<double> > * results)	//allows to calculate at other points
{
	if ( white_Green )
	{
		tau_integration_WhiteGreen(Gtilde_n_white, Gtilde_n_white, k, results);
	}
	else	//otherwise, use colored Green's function
	{
		if ( white_noise )
			tau_integration_coloredGreen(Gtilde_n_color, Gtilde_n_color, k, results);
		else
			colored_tau_integration(Gtilde_n_color, Gtilde_n_color, k, results);
	}

	(*results)[0] /= tauf*tauf;	//fm^2; note that we must include extra factor of tauf^-2 for consistency with manuscript
	(*results)[1] /= tauf*tauf;	//fm^2; note that we must include extra factor of tauf^-2 for consistency with manuscript
	(*results)[2] /= tauf*tauf;	//do the same for self-correlations...

	return;
}

//Miscellaneous stuff
void replace_character(std::string & tempstring, char swap_out, char swap_in)
{
	int len = tempstring.length();
	for(unsigned int i = 0; i < len; i++)
	{
		char c = tempstring[i];
		if (c == swap_out)
			tempstring[i] = swap_in;
	}
	
	if (tempstring[len - 1] == swap_in)
		tempstring.erase( len - 1 );
	
	return;
}

string get_vQ2_string()
{
	ostringstream vQ2stream;
	vQ2stream << setprecision(3) << std::fixed << vQ2;
	string vQ2string = vQ2stream.str();
	replace_character(vQ2string, '.', '_');

	return (vQ2string);
}

string get_snapshot_string()
{
	ostringstream snapshotstream;
	snapshotstream << setprecision(2) << std::fixed << fraction_of_evolution;
	string snapshotstring = snapshotstream.str();
	replace_character(snapshotstring, '.', '_');

	return (snapshotstring);
}

void set_outfilenames(
						string & mainResultsFilename,
						string & dndn_k_Filename,
						string & dndn_Dxi_Filename
						)
{
	string noiseStem = white_noise ? "white" : "color";
	string GreenStem = white_Green ? "white" : "color";

	string vQ2string = ( white_noise and white_Green ) ? "" : "_vQ2_" + get_vQ2_string();
	string snapshotstring = "_snapshot_" + get_snapshot_string();

	mainResultsFilename = resultsPath + "/twoPC_" + particle1.name + particle2.name + "_"
							+ noiseStem + "_noise_" + GreenStem + "_Green" + vQ2string + snapshotstring + ".dat";
	dndn_k_Filename = resultsPath + "/dndn_k_" + particle1.name + particle2.name + "_"
							+ noiseStem + "_noise_" + GreenStem + "_Green" + vQ2string + snapshotstring + ".dat";
	dndn_Dxi_Filename = resultsPath + "/dndn_Dxi_" + particle1.name + particle2.name + "_"
							+ noiseStem + "_noise_" + GreenStem + "_Green" + vQ2string + snapshotstring + ".dat";

	return;
}
*/

// End of file

#endif
