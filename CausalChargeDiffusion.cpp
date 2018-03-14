#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <complex>
#include <algorithm>

using namespace std;

#include "resonance/sfn.h"
#include "CCDparams.h"
#include "defs.h"

using namespace CCD;

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
	//Tf = 106.4226 / hbarC;	//old FO temperature, for comparison

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

	//output some parameters from calculation
	ofstream output_results;
	output_results.open(filename1.c_str());

	output_results 
			<< "#########################################################" << endl
			<< "# Using following parameters:" << endl
			<< "# resultsPath = " << resultsPath << endl
			<< "# fraction_of_evolution = " << fraction_of_evolution << endl
			<< "# taui = " << taui << " fm/c, tauf = " << tauf << " fm/c" << endl
			<< "# si = " << si << ", sf = " << sf << endl
			<< "# Ti = " << Ti*hbarC << " MeV, Tf = " << Tf*hbarC << " MeV" << endl
			<< "# v_Q^2 = " << vQ2 << ", D_Q = " << DQ << " fm/c, tau_Q = " << tauQ << " fm/c" << endl
			<< "# chi_{T,T} = " << chi_T_T << ", chi_{T,mu} = chi_{mu,T} = " << chi_T_mu << ", chi_{mu,mu} = " << chi_mu_mu << endl
			<< "#########################################################" << endl;

    // set up grid points for integrations
    xi_pts_minf_inf = new double [n_xi_pts];
    tau_pts = new double [n_tau_pts];
    k_pts = new double [n_k_pts];
    xi_wts_minf_inf = new double [n_xi_pts];
    tau_wts = new double [n_tau_pts];
    k_wts = new double [n_k_pts];

    int tmp = gauss_quadrature(n_xi_pts, 1, 0.0, 0.0, -xi_infinity, xi_infinity, xi_pts_minf_inf, xi_wts_minf_inf);
    //tmp = gauss_quadrature(n_k_pts, 1, 0.0, 0.0, -k_infinity, k_infinity, k_pts, k_wts);
    tmp = gauss_quadrature(n_k_pts, 1, 0.0, 0.0, 0.0, k_infinity, k_pts, k_wts);
    tmp = gauss_quadrature(n_tau_pts, 1, 0.0, 0.0, taui, tauf, tau_pts, tau_wts);

	T_pts = new double [n_tau_pts];

	//computes tau-dependence of T for remainder of calculation
	for (int it = 0; it < n_tau_pts; ++it)
		T_pts[it] = guess_T(tau_pts[it]);

	//////////////////////////////////////////////////////////
	//Check hypergeometric function implementations
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
	//Check colored Green's function implementation
	/*
	for (int ik = 0; ik < n_k_pts; ++ik)
	for (int it = 0; it < n_tau_pts; ++it)
	{
		double k = k_pts[ik];
		double tau = tau_pts[it];
		double x = tau / tauQ;
		complex<double> Gt = Gtilde_n_color(k, tauf, tau);
		complex<double> aGt = asymptotic_Gtilde_n_color(k, tauf, tau);
		cout << setprecision(20) << k << "   " << x << "   " << Gt.real() << "   " << Gt.imag() << "   " << aGt.real() << "   " << aGt.imag() << endl;
	}
	if (1) return (0);
	*/
	//////////////////////////////////////////////////////////

	//get the ensemble averaged spectra
	double norm = integrate_1D(norm_int, xi_pts_minf_inf, xi_wts_minf_inf, n_xi_pts, &particle1);	//by definition of charge balance function (CBF)

//cout << "norm = " << norm << endl;
//if (1) exit (8);

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

		if (print_dndn_k)
			output_dndn_k << k << "   " << setprecision(15)
							<< results[0].real() << "   "
							<< results[1].real() << "   "
							<< results[2].real() << endl;
	}

	//start computing actual charge balance functions here
	for (int iDy = 0; iDy < n_Dy; iDy++)
	{
		double Delta_y = (double)iDy * Delta_y_step;
		double Delta_xi = Delta_y;

		complex<double> sum(0,0), sum_no_SC(0,0);
		complex<double> sum_Dxi(0,0), sum_Dxi_no_SC(0,0), SC_Dxi(0,0);
		for (int ik = 0; ik < n_k_pts; ++ik)
		{
			double k = k_pts[ik];

			complex<double> Ftn1 = Ftn_particle1_vec[ik];
			complex<double> Ftn2 = Ftn_particle2_vec[ik];
			complex<double> Ctnn = Ctnn_vec[ik];
			complex<double> Ctnn_no_SC = Ctnn_no_SC_vec[ik];

			//sum += k_wts[ik] * exp(i * k * Delta_y)
			//		* ( Ftn1 * conj(Ftn2) * Ctnn );
			//sum_no_SC += k_wts[ik] * exp(i * k * Delta_y)
			//		* ( Ftn1 * conj(Ftn2) * Ctnn_no_SC );
			sum += 2.0 * k_wts[ik] * cos(k * Delta_y)
					* ( Ftn1 * conj(Ftn2) * Ctnn );
			sum_no_SC += 2.0 * k_wts[ik] * cos(k * Delta_y)
					* ( Ftn1 * conj(Ftn2) * Ctnn_no_SC );

			//sum_Dxi += k_wts[ik] * exp(i * k * Delta_xi) * Ctnn;
			//sum_Dxi_no_SC += k_wts[ik] * exp(i * k * Delta_xi) * Ctnn_no_SC;
			//SC_Dxi += k_wts[ik] * exp(i * k * Delta_xi) * SC_vec[ik];
			sum_Dxi += 2.0 * k_wts[ik] * cos(k * Delta_xi) * Ctnn;
			sum_Dxi_no_SC += 2.0 * k_wts[ik] * cos(k * Delta_xi) * Ctnn_no_SC;
			SC_Dxi += 2.0 * k_wts[ik] * cos(k * Delta_xi) * SC_vec[ik];
		}

		//omit smearing function factors to just F.T. back to xi-space
		if (print_dndn_Dxi)
			output_dndn_Dxi << Delta_xi << "   " << setprecision(15)
							<< sum_Dxi.real() << "   "
							<< sum_Dxi_no_SC.real() << "   "
							<< SC_Dxi.real() << endl;

		complex<double> result = (ds*tauf*Tf / (4.0*M_PI*M_PI * norm)) * sum;
		complex<double> result_no_SC = (ds*tauf*Tf / (4.0*M_PI*M_PI * norm)) * sum_no_SC;

		output_results
				<< setprecision(15) << Delta_y << "   "
				<< result.real() << "   " << result.imag() << "   "
				<< result_no_SC.real() << "   " << result_no_SC.imag() << endl;
	}

	//one last check...
	bool run_final_check = true;
	if (run_final_check)
	{
		vector<double> Delta_y_pts(n_Dy);
		for (int iDy = 0; iDy < n_Dy; iDy++)
			Delta_y_pts[iDy] = Delta_y_min + (double)iDy * Delta_y_step;

		vector<double> k_pts_vec(n_k_pts);
		for (int ik = 0; ik < n_k_pts; ++ik)
			k_pts_vec[ik] = k_pts[ik];

		vector<complex<double> > sFn = sfn::operation_mode_1(Delta_y_pts, k_pts_vec, false);
		for (int iDy = 0; iDy < n_Dy; iDy++)
		{
			//double Delta_y = Delta_y_min + (double)iDy * Delta_y_step;
			double Delta_y = Delta_y_pts[iDy];
			for (int ik = 0; ik < n_k_pts; ++ik)
			{
//		if (iDy > 0)
//			continue;
//		if (ik > 0)
//			continue;
				double k = k_pts[ik];

				complex<double> Ftn1 = Ftn_particle1_vec[ik];
				complex<double> Ftn2 = Ftn_particle2_vec[ik];

				complex<double> z = exp(i * k * Delta_y) * Ftn1 * conj(Ftn2);

				int idx_ik_iDpY0 = ik * n_Dy + (n_Dy-1)/2;
				int idx_ik_iDpY = ik * n_Dy + iDy;
				complex<double> zcomp
					= sFn[idx_ik_iDpY0] * sFn[idx_ik_iDpY];

				cout << "sFn: "
						<< Delta_y << "   " << k << "   "
						<< z.real() << "   " << z.imag() << "   "
						<< zcomp.real() << "   " << zcomp.imag() << endl;
			}
		}
	}


	if (print_dndn_k)
		output_dndn_k.close();
	if (print_dndn_Dxi)
		output_dndn_Dxi.close();
	output_results.close();

	return 0;
}
