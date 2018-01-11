#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <complex>
#include <algorithm>
#include <cmath>

using namespace std;

#include "decay.h"

extern vector<double> pT_pts, pT_wts;
extern vector<double> pphi_pts, pphi_wts;
extern vector<double> Del_pY_pts, Del_pY_wts;

namespace decay
{
	const double PTCHANGE = 1.0;	//GeV
	int target_pid = -1, particle_monval = -1;
	int n_resonance = -1;
	vector<double> current_resonance_decay_masses;
	double current_resonance_mass, current_resonance_mu, current_resonance_Gamma, current_resonance_total_br, current_resonance_direct_br;
	double current_daughter_mass, current_daughter_Gamma, current_m2_Gamma, current_m3_Gamma;
	vector<double> spectra_re, spectra_im;
	int current_reso_nbody = -1, current_parent_resonance = -1, current_resonance_pid = -1;
	int current_decay_channel_idx = -1, current_resonance_idx = -1;
	string current_decay_channel_string = "";

	const double vmin = -1.0, vmax = 1.0;
	const double zetamin = 0.0, zetamax = M_PI;
	vector<double> x_pts, x_wts, s_pts, s_wts, v_pts, v_wts, zeta_pts, zeta_wts;

	vector<readindata::particle_info> all_particles;
	vector<decay::decay_info> decay_channels;
	vector<int> chosen_resonance_indices;

	bool thermal_pions_only = false;
	double VEC_n2_s_factor;
	vector<double> VEC_n2_P_Y, VEC_n2_v_factor, VEC_n2_zeta_factor, VEC_n2_PPhi_tilde, VEC_n2_PPhi_tildeFLIP, VEC_n2_PT;
	vector<double> VEC_n3_P_Y, VEC_n3_s_factor, VEC_n3_v_factor, VEC_n3_zeta_factor, VEC_n3_PPhi_tilde, VEC_n3_PPhi_tildeFLIP, VEC_n3_PT;
	vector<vector<double> > VEC_n2_Ppm, VEC_n3_Ppm;

	double m, m2, m3, M, br, Qfunc, mT, pT, Gamma;

	/////////////////////////////////////////////
	//end of declarations
	/////////////////////////////////////////////
	
	void Compute_phase_space_integrals(vector<readindata::particle_info> all_particles_in,
										vector<int> chosen_resonance_indices_in,
										vector<double> * spectra_re_in,
										vector<double> * spectra_im_in,
										int target_particle_id_in)
	{
		all_particles = all_particles_in;

		//load decay channel information
		int n_decay_channels = load_decay_channel_info(all_particles, chosen_resonance_indices_in);
		chosen_resonance_indices = chosen_resonance_indices_in;

		target_pid = target_particle_id_in;

		if (thermal_pions_only)
		{
			cout << "Thermal pions only: no phase-space integrals need to be computed." << endl;
			return;
		}

		Allocate_decay_channel_info();

		for (int idc = 1; idc <= n_decay_channels; ++idc)
		{
			if (decay_channels[idc-1].resonance_particle_id == target_pid || thermal_pions_only)
				break;
			else if (!Do_this_decay_channel(idc))
				continue;
	
			Set_current_particle_info(idc);

			for (int idc_DI = 0; idc_DI < current_reso_nbody; ++idc_DI)
			{
				int daughter_resonance_particle_id = -1;
				if (!Do_this_daughter_particle(idc, idc_DI, &daughter_resonance_particle_id))
					continue;

				Set_current_daughter_info(idc, idc_DI);

				Do_resonance_integrals( idc, current_resonance_pid,
										daughter_resonance_particle_id,
										spectra_re_in, spectra_im_in );
			}
		}			
									// END of decay channel loop
		return;
	}

	void Do_resonance_integrals(int decay_channel,
								int parent_pid, int daughter_pid,
								vector<double> * spec_re_ptr,
								vector<double> * spec_im_ptr)
	{
		int n_body = current_reso_nbody;

		if (n_body == 2)
			two_body_decay( decay_channel,
							parent_pid, daughter_pid,
							spec_re_ptr, spec_im_ptr );
		else
			three_body_decay( decay_channel,
								parent_pid, daughter_pid,
								spec_re_ptr, spec_im_ptr );

		return;
	}

	double get_Q()
	{
		double smin = (m2+m3)*(m2+m3);
		double smax = (M-m)*(M-m);
		double sum = 0.;
	
		for (int is = 0; is < n_s_pts; ++is)
		{
			double sp = s_pts[is];
			double f1 = (M+m)*(M+m) - sp;
			double f2 = smax - sp;
			double f3 = smin - sp;
			double f4 = (m2-m3)*(m2-m3) - sp;
			sum += s_wts[is]*sqrt(f1*f2*f3*f4)/(sp+1.e-15);
		}

		return sum;
	}

	double g(double s)
	{
		int n_body = current_reso_nbody;
		double gs_pstar = sqrt( ((M+m)*(M+m) - s)*((M-m)*(M-m) - s) )/(2.0*M);
		double g_res = br/(4.*M_PI*gs_pstar);
		if (n_body == 3 || n_body == 4)		//both set up to work the same way
		{
			double pre_f = (M*br)/(2.*M_PI*s);
			double num = sqrt( (s - (m2+m3)*(m2+m3)) * (s - (m2-m3)*(m2-m3)) );
			double den = Qfunc;
			g_res = pre_f * num / den;
		}

		return g_res;
	}

	///////////////////////////////////////////////////////////////
	void two_body_decay( int dc_idx, int parent_pid, int daughter_pid,
							vector<double> * spec_re, vector<double> * spec_im )
	{
		//useful for interpolation
		vector<double> log_spec_re(spec_re->size());
		vector<double> log_spec_im(spec_im->size());
		vector<double> sign_spec_re(spec_re->size());
		vector<double> sign_spec_im(spec_im->size());
		for (int i = 0; i < (int)spec_re->size(); ++i)
		{
			log_spec_re[i] = log(abs(spec_re->at(i))+1.e-100);
			sign_spec_re[i] = sgn(spec_re->at(i));
		}
		for (int i = 0; i < (int)spec_im->size(); ++i)
		{
			log_spec_im[i] = log(abs(spec_im->at(i))+1.e-100);
			sign_spec_im[i] = sgn(spec_im->at(i));
		}

		//set indices in arrays to know which elements to access/update
		int parent_idx = lookup_resonance_idx_from_particle_id( parent_pid );
		int daughter_idx = lookup_resonance_idx_from_particle_id( daughter_pid );

		//loop over daughter momentum grid
		for (int ipT = 0; ipT < n_pT_pts; ++ipT)
		for (int ipY = 0; ipY < n_pY_pts; ++ipY)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		{
			double local_pT = pT_pts[ipT];
			double local_pphi = pphi_pts[ipphi];
			double local_pY = Del_pY_pts[ipY];
			Load_decay_channel_info_nb2(dc_idx, local_pT, local_pphi, local_pY);	// set decay channel information

			//then g(s) is delta-function, skip s-integration entirely
			//double s = m2*m2;
			double vsum_re = 0.0, vsum_im = 0.0;
			for (int iv = 0; iv < n_v_pts; ++iv)
			{
				double zetasum_re = 0.0, zetasum_im = 0.0;
				for (int izeta = 0; izeta < n_zeta_pts; ++izeta)
				{
					double Csum_re = 0.0, Csum_im = 0.0;
					double PTR = VEC_n2_PT[NB2_indexer(iv,izeta)];
					double PYR = VEC_n2_P_Y[iv];
					double PphiR = VEC_n2_PPhi_tilde[NB2_indexer(iv,izeta)];

					for (int tempidx = 0; tempidx <= 1; ++tempidx)
					{
						if (tempidx != 0)
							PphiR = VEC_n2_PPhi_tildeFLIP[NB2_indexer(iv,izeta)];		//also takes Pp --> Pm

						//spectra
						Edndp3(PTR, PphiR, PYR, parent_idx, &Csum_re, spec_re, &log_spec_re, &sign_spec_re);
						Edndp3(PTR, PphiR, PYR, parent_idx, &Csum_im, spec_im, &log_spec_im, &sign_spec_im);
					}												// end of tempidx sum
					zetasum_re += VEC_n2_zeta_factor[NB2_indexer(iv,izeta)]*Csum_re;
					zetasum_im += VEC_n2_zeta_factor[NB2_indexer(iv,izeta)]*Csum_im;
				}													// end of zeta sum
				vsum_re += VEC_n2_v_factor[iv]*zetasum_re;
				vsum_im += VEC_n2_v_factor[iv]*zetasum_im;
			}														// end of v sum
			double ssum_re = M*VEC_n2_s_factor*vsum_re;
			double ssum_im = M*VEC_n2_s_factor*vsum_im;

			(*spec_re)[res_vector_indexer(daughter_idx, ipT, ipphi, ipY)] += ssum_re;
			(*spec_im)[res_vector_indexer(daughter_idx, ipT, ipphi, ipY)] += ssum_im;
		}											// end of pT, pphi, pY loops
		return;
	}




	///////////////////////////////////////////////////////////////
	void three_body_decay( int dc_idx, int parent_pid, int daughter_pid,
							vector<double> * spec_re, vector<double> * spec_im )
	{
		//useful for interpolation
		vector<double> log_spec_re(spec_re->size());
		vector<double> log_spec_im(spec_im->size());
		vector<double> sign_spec_re(spec_re->size());
		vector<double> sign_spec_im(spec_im->size());

		for (int i = 0; i < (int)spec_re->size(); ++i)
		{
			log_spec_re[i] = log(abs(spec_re->at(i))+1.e-100);
			sign_spec_re[i] = sgn(spec_re->at(i));
		}
		for (int i = 0; i < (int)spec_im->size(); ++i)
		{
			log_spec_im[i] = log(abs(spec_im->at(i))+1.e-100);
			sign_spec_im[i] = sgn(spec_im->at(i));
		}

		//set indices in arrays to know which elements to access/update
		int parent_idx = lookup_resonance_idx_from_particle_id( parent_pid );
		int daughter_idx = lookup_resonance_idx_from_particle_id( daughter_pid );

		//mT = sqrt(m*m + pT*pT);
		//double s_min = (m2 + m3)*(m2 + m3);
		//double s_max = (M - m)*(M - m);
		//Qfunc = get_Q();

		//loop over daughter momentum grid
		for (int ipT = 0; ipT < n_pT_pts; ++ipT)
		for (int ipY = 0; ipY < n_pY_pts; ++ipY)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		{
			double local_pT = pT_pts[ipT];
			double local_pphi = pphi_pts[ipphi];
			double local_pY = Del_pY_pts[ipY];
			Load_decay_channel_info_nb3(dc_idx, local_pT, local_pphi, local_pY);

			double ssum_re = 0.0, ssum_im = 0.0;
			for (int is = 0; is < n_s_pts; ++is)
			{
				double vsum_re = 0.0, vsum_im = 0.0;
				for (int iv = 0; iv < n_v_pts; ++iv)
				{
					double zetasum_re = 0.0, zetasum_im = 0.0;
					for (int izeta = 0; izeta < n_zeta_pts; ++izeta)
					{
						double Csum_re = 0.0, Csum_im = 0.0;
						double PTR = VEC_n3_PT[NB3_indexer(is,iv,izeta)];
						double PYR = VEC_n3_P_Y[is*n_v_pts+iv];
						double PphiR = VEC_n3_PPhi_tilde[NB3_indexer(is,iv,izeta)];

						for (int tempidx = 0; tempidx <= 1; ++tempidx)
						{
							if (tempidx != 0)
								PphiR = VEC_n3_PPhi_tildeFLIP[NB3_indexer(is,iv,izeta)];		//also takes Pp --> Pm

							//spectra
							Edndp3(PTR, PphiR, PYR, parent_idx, &Csum_re, spec_re, &log_spec_re, &sign_spec_re);
							Edndp3(PTR, PphiR, PYR, parent_idx, &Csum_im, spec_im, &log_spec_im, &sign_spec_im);
						}										// end of tempidx sum
						zetasum_re += VEC_n3_zeta_factor[NB3_indexer(is,iv,izeta)]*Csum_re;
						zetasum_im += VEC_n3_zeta_factor[NB3_indexer(is,iv,izeta)]*Csum_im;
					}											// end of zeta sum
					vsum_re += VEC_n3_v_factor[is*n_v_pts+iv]*zetasum_re;
					vsum_im += VEC_n3_v_factor[is*n_v_pts+iv]*zetasum_im;
				}												// end of v sum
				ssum_re += M*VEC_n3_s_factor[is]*vsum_re;
				ssum_im += M*VEC_n3_s_factor[is]*vsum_im;
			}													// end of s sum

			(*spec_re)[res_vector_indexer(daughter_idx, ipT, ipphi, ipY)] += ssum_re;
			(*spec_im)[res_vector_indexer(daughter_idx, ipT, ipphi, ipY)] += ssum_im;
		}								// end of pT, pphi, pY loops
		return;
	}


	///////////////////////////////////////////////////////
	///////////////////////////////////////////////////////
	void Load_decay_channel_info_nb2(int dc_idx, double pT, double pphi, double p_y)
	{
		//set decay information
		decay_info * current_decay = &decay_channels[dc_idx];

		M = current_resonance_mass;
		Gamma = current_resonance_Gamma;
		//one_by_Gamma_Mres = 1./(Gamma*M + 1.e-25);	//keeps calculation safe when Gamma == 0
		//N.B. - no need for hbarc, since this will only multiply something with GeV^2 units in the end
		m = current_daughter_mass;
		br = current_resonance_direct_br;	//doesn't depend on target daughter particle, just parent resonance and decay channel
		m2 = current_resonance_decay_masses[0];
		m3 = current_resonance_decay_masses[1];
		int n_body = current_decay->nbody;

		// some particles may decay to particles with more total mass than originally
		// --> broaden with resonance widths
		while ((m + m2) > M)
		{
			M += 0.25 * Gamma;
			m -= 0.5 * current_daughter_Gamma;
			m2 -= 0.5 * current_m2_Gamma;
		}

		mT = sqrt(m*m + pT*pT);

		//set up vectors of points to speed-up integrals...
		double s = m2*m2;
		double pstar = sqrt( ((M+m)*(M+m) - s)*((M-m)*(M-m) - s) )/(2.0*M);
		double g_s = g(s);	//for n_body == 2, doesn't actually use s since result is just a factor * delta(...); just returns factor
		double Estar = sqrt(m*m + pstar*pstar);
		double psBmT = pstar / mT;
		double DeltaY = log(psBmT + sqrt(1.+psBmT*psBmT));

		VEC_n2_s_factor = br/(4.*M_PI*pstar);	//==g_s

		reset(&VEC_n2_v_factor);
		reset(&VEC_n2_zeta_factor);
		reset(&VEC_n2_PT);
		reset(&VEC_n2_PPhi_tilde);
		reset(&VEC_n2_PPhi_tildeFLIP);
		reset(&VEC_n2_P_Y);

		for(int iv = 0; iv < n_v_pts; ++iv)
		{
			double v = v_pts[iv];
			double P_Y = p_y + v*DeltaY;
			double mT_ch_P_Y_p_y = mT*cosh(v*DeltaY);
			double x2 = mT_ch_P_Y_p_y*mT_ch_P_Y_p_y - pT*pT;
			double MTbar = Estar*M*mT_ch_P_Y_p_y/x2;
			double DeltaMT = M*pT*sqrt(Estar*Estar - x2)/x2;

			VEC_n2_P_Y[iv] = P_Y;
			VEC_n2_v_factor[iv] = v_wts[iv]*DeltaY/sqrt(x2);

			for(int izeta = 0; izeta < n_zeta_pts; ++izeta)
			{
				double zeta = zeta_pts[izeta];
				double MT = MTbar + cos(zeta)*DeltaMT;
				double PT = sqrt(MT*MT - M*M);
				double temp_cos_PPhi_tilde = (mT*MT*cosh(P_Y-p_y) - Estar*M)/(pT*PT);
				//assume that PPhi_tilde is +ve in next step...
				double temp_sin_PPhi_tilde = sqrt(1. - temp_cos_PPhi_tilde*temp_cos_PPhi_tilde);
				double PPhi_tilde = place_in_range( atan2(temp_sin_PPhi_tilde, temp_cos_PPhi_tilde), pphi_min, pphi_max);

				VEC_n2_zeta_factor[NB2_indexer(iv, izeta)] = zeta_wts[izeta]*MT;
				VEC_n2_PPhi_tilde[NB2_indexer(iv, izeta)] = place_in_range( pphi + PPhi_tilde, pphi_min, pphi_max);
				VEC_n2_PPhi_tildeFLIP[NB2_indexer(iv, izeta)] = place_in_range( pphi - PPhi_tilde, pphi_min, pphi_max);
				VEC_n2_PT[NB2_indexer(iv, izeta)] = PT;
				//set P^+ components
				VEC_n2_Ppm[NB2_indexer(iv, izeta)*2+0][0] = MT * cosh(P_Y);
				VEC_n2_Ppm[NB2_indexer(iv, izeta)*2+0][1] = PT * cos(pphi + PPhi_tilde);
				VEC_n2_Ppm[NB2_indexer(iv, izeta)*2+0][2] = PT * sin(pphi + PPhi_tilde);
				VEC_n2_Ppm[NB2_indexer(iv, izeta)*2+0][3] = MT * sinh(P_Y);
				//set P^- components
				VEC_n2_Ppm[NB2_indexer(iv, izeta)*2+1][0] = MT * cosh(P_Y);
				VEC_n2_Ppm[NB2_indexer(iv, izeta)*2+1][1] = PT * cos(pphi - PPhi_tilde);
				VEC_n2_Ppm[NB2_indexer(iv, izeta)*2+1][2] = PT * sin(pphi - PPhi_tilde);
				VEC_n2_Ppm[NB2_indexer(iv, izeta)*2+1][3] = MT * sinh(P_Y);
			}
		}

		return;
	}

	///////////////////////////////////////////////////////
	///////////////////////////////////////////////////////
	void Load_decay_channel_info_nb3(int dc_idx, double pT, double pphi, double p_y)
	{
		//set decay information
		decay_info * current_decay = &decay_channels[dc_idx];

		M = current_resonance_mass;
		Gamma = current_resonance_Gamma;
		//one_by_Gamma_Mres = 1./(Gamma*M + 1.e-25);	//keeps calculation safe when Gamma == 0
		//N.B. - no need for hbarc, since this will only multiply something with GeV^2 units in the end
		m = current_daughter_mass;
		br = current_resonance_direct_br;	//doesn't depend on target daughter particle, just parent resonance and decay channel
		m2 = current_resonance_decay_masses[0];
		m3 = current_resonance_decay_masses[1];
		int n_body = current_decay->nbody;

		mT = sqrt(m*m + pT*pT);
		double s_min_temp = (m2 + m3)*(m2 + m3);
		double s_max_temp = (M - m)*(M - m);
		double cen = 0.5 * (s_max_temp + s_min_temp);
		double hw = 0.5 * (s_max_temp - s_min_temp);
		//gauss_quadrature(n_s_pts, 1, 0.0, 0.0, s_min_temp, s_max_temp, s_pts, s_wts);
		for (int is = 0; is < n_s_pts; ++is)
		{
			s_pts[is] = cen + hw * x_pts[is];
			s_wts[is] = hw * x_wts[is];
		}
		Qfunc = get_Q();

		reset(&VEC_n3_s_factor);
		reset(&VEC_n3_v_factor);
		reset(&VEC_n3_zeta_factor);
		reset(&VEC_n3_PT);
		reset(&VEC_n3_PPhi_tilde);
		reset(&VEC_n3_PPhi_tildeFLIP);
		reset(&VEC_n3_P_Y);

		// s-loop
		for (int is = 0; is < n_s_pts; ++is)
		{
			double s = s_pts[is];
			double g_s = g(s);
			double pstar = sqrt(((M+m)*(M+m) - s)*((M-m)*(M-m) - s))/(2.0*M);
			double Estar = sqrt(m*m + pstar*pstar);
			double psBmT = pstar / mT;
			double DeltaY = log(psBmT + sqrt(1.+psBmT*psBmT));

			VEC_n3_s_factor[is] = s_wts[is]*g_s;

			// v-loop
			for(int iv = 0; iv < n_v_pts; ++iv)
			{
				double v = v_pts[iv];
				double P_Y = p_y + v*DeltaY;
				double mT_ch_P_Y_p_y = mT*cosh(v*DeltaY);
				double x2 = mT_ch_P_Y_p_y*mT_ch_P_Y_p_y - pT*pT;
				double MTbar = Estar*M*mT_ch_P_Y_p_y/x2;
				double DeltaMT = M*pT*sqrt(Estar*Estar - x2)/x2;

				VEC_n3_P_Y[is * n_v_pts + iv] = P_Y;
				VEC_n3_v_factor[is * n_v_pts + iv] = v_wts[iv]*DeltaY/sqrt(x2);

				// zeta-loop
				for(int izeta = 0; izeta < n_zeta_pts; ++izeta)
				{
					double zeta = zeta_pts[izeta];
					double MT = MTbar + cos(zeta)*DeltaMT;
					double PT = sqrt(MT*MT - M*M);
					double temp_cos_PPhi_tilde = (mT*MT*cosh(P_Y-p_y) - Estar*M)/(pT*PT);
					//assume that PPhi_tilde is +ve in next step...
					double temp_sin_PPhi_tilde = sqrt(1. - temp_cos_PPhi_tilde*temp_cos_PPhi_tilde);
					double PPhi_tilde = place_in_range( atan2(temp_sin_PPhi_tilde, temp_cos_PPhi_tilde), pphi_min, pphi_max);

					VEC_n3_zeta_factor[NB3_indexer(is, iv, izeta)] = zeta_wts[izeta]*MT;
					VEC_n3_PPhi_tilde[NB3_indexer(is, iv, izeta)] = place_in_range( pphi + PPhi_tilde, pphi_min, pphi_max);
					VEC_n3_PPhi_tildeFLIP[NB3_indexer(is, iv, izeta)] = place_in_range( pphi - PPhi_tilde, pphi_min, pphi_max);
					VEC_n3_PT[NB3_indexer(is, iv, izeta)] = PT;
					//set P^+ components
					VEC_n3_Ppm[NB3_indexer(is, iv, izeta)*2+0][0] = MT * cosh(P_Y);
					VEC_n3_Ppm[NB3_indexer(is, iv, izeta)*2+0][1] = PT * cos(pphi + PPhi_tilde);
					VEC_n3_Ppm[NB3_indexer(is, iv, izeta)*2+0][2] = PT * sin(pphi + PPhi_tilde);
					VEC_n3_Ppm[NB3_indexer(is, iv, izeta)*2+0][3] = MT * sinh(P_Y);
					//set P^- components
					VEC_n3_Ppm[NB3_indexer(is, iv, izeta)*2+1][0] = MT * cosh(P_Y);
					VEC_n3_Ppm[NB3_indexer(is, iv, izeta)*2+1][1] = PT * cos(pphi - PPhi_tilde);
					VEC_n3_Ppm[NB3_indexer(is, iv, izeta)*2+1][2] = PT * sin(pphi - PPhi_tilde);
					VEC_n3_Ppm[NB3_indexer(is, iv, izeta)*2+1][3] = MT * sinh(P_Y);
				}
			}
		}

		return;
	}

	///////////////////////////////////////////////////////
	///////////////////////////////////////////////////////
	void Edndp3(double ptr, double pphir, double pyr, int parent_idx, double * result,
				vector<double> * loc_spectra_ptr, vector<double> * loc_log_spectra_ptr,
				vector<double> * loc_sign_spectra_ptr)
	{
		double phi0, phi1, py0, py1;
		double f11, f21, f12, f22;

		int npT_max = n_pT_pts - 1;
		int npphi_max = n_pphi_pts - 1;
		int npy_max = n_pY_pts - 1;

		// locate pT interval
		int npt = 1;
		while ((ptr > pT_pts[npt]) &&
				(npt < npT_max)) ++npt;
		double pT0 = pT_pts[npt-1];
		double pT1 = pT_pts[npt];

		// locate pphi interval
		int nphi = 1, nphim1 = 0;
		if(pphir < pphi_pts[0])			//if angle is less than minimum angle grid point
		{
			phi0 = pphi_pts[npphi_max] - 2. * M_PI;
			phi1 = pphi_pts[0];
			nphi = 0;
			nphim1 = npphi_max;
		}
		else if(pphir > pphi_pts[npphi_max])	//if angle is greater than maximum angle grid point
		{
			phi0 = pphi_pts[npphi_max];
			phi1 = pphi_pts[0] + 2. * M_PI;
			nphi = 0;
			nphim1 = npphi_max;
		}
		else						//if angle is within grid range
		{
			while ((pphir > pphi_pts[nphi]) &&
					(nphi < npphi_max)) ++nphi;
			nphim1 = nphi - 1;
			phi0 = pphi_pts[nphim1];
			phi1 = pphi_pts[nphi];
		}

		// locate py interval
		long npy = 1, npym1 = 0;
		if(pyr > Del_pY_max)	//if rapidity is greater than maximum rapidity grid point
		{
			py0 = Del_pY_max;
			py1 = Del_pY_max;
			if (VERBOSE > 2) cerr << "WARNING in Edndp3(): " << pyr << " > " << Del_pY_max << endl;
			return;
		}
		else if(pyr < Del_pY_min)	//else if rapidity is less than minimum rapidity grid point
		{								//this can't happen when USE_RAPIDITY_SYMMETRY is true, since pyr = abs(spyr) >= 0
			py0 = Del_pY_min;
			py1 = Del_pY_min;
			if (VERBOSE > 2) cerr << "WARNING in Edndp3(): " << pyr << " < " << Del_pY_min << endl;
			return;
		}
		else						//if rapidity is within grid range
		{
			while ((pyr > Del_pY_pts[npy]) &&
					(npy < npy_max)) ++npy;
			npym1 = npy - 1;
			py0 = Del_pY_pts[npym1];
			py1 = Del_pY_pts[npy];
		}

		if (pT0==pT1 || phi0==phi1 || py0==py1)
		{
			cerr << "ERROR in Edndp3(): pT, pphi and/or py values equal!" << endl;
			exit(1);
		}

		double one_by_pTdiff = 1./(pT1 - pT0), one_by_pphidiff = 1./(phi1 - phi0), one_by_pYdiff = 1./(py1 - py0 + 1.e-10);
		double del_ptr_pt0 = ptr - pT0, del_phir_phi0 = pphir - phi0, del_pyr_py0 = pyr - py0;

		// choose pt-pphi slice of resonance info arrays
		double log_f111 = (*loc_log_spectra_ptr)[res_vector_indexer(parent_idx,npt-1,nphim1,npym1)];
		double log_f121 = (*loc_log_spectra_ptr)[res_vector_indexer(parent_idx,npt-1,nphi,npym1)];
		double log_f211 = (*loc_log_spectra_ptr)[res_vector_indexer(parent_idx,npt,nphim1,npym1)];
		double log_f221 = (*loc_log_spectra_ptr)[res_vector_indexer(parent_idx,npt,nphi,npym1)];
		double log_f112 = (*loc_log_spectra_ptr)[res_vector_indexer(parent_idx,npt-1,nphim1,npy)];
		double log_f122 = (*loc_log_spectra_ptr)[res_vector_indexer(parent_idx,npt-1,nphi,npy)];
		double log_f212 = (*loc_log_spectra_ptr)[res_vector_indexer(parent_idx,npt,nphim1,npy)];
		double log_f222 = (*loc_log_spectra_ptr)[res_vector_indexer(parent_idx,npt,nphi,npy)];

		double f111 = (*loc_spectra_ptr)[res_vector_indexer(parent_idx,npt-1,nphim1,npym1)];
		double f121 = (*loc_spectra_ptr)[res_vector_indexer(parent_idx,npt-1,nphi,npym1)];
		double f211 = (*loc_spectra_ptr)[res_vector_indexer(parent_idx,npt,nphim1,npym1)];
		double f221 = (*loc_spectra_ptr)[res_vector_indexer(parent_idx,npt,nphi,npym1)];
		double f112 = (*loc_spectra_ptr)[res_vector_indexer(parent_idx,npt-1,nphim1,npy)];
		double f122 = (*loc_spectra_ptr)[res_vector_indexer(parent_idx,npt-1,nphi,npy)];
		double f212 = (*loc_spectra_ptr)[res_vector_indexer(parent_idx,npt,nphim1,npy)];
		double f222 = (*loc_spectra_ptr)[res_vector_indexer(parent_idx,npt,nphi,npy)];

		/////////////////////////////////////////////////////////////////
		// interpolate over pT values first
		/////////////////////////////////////////////////////////////////
		if(ptr > PTCHANGE)				// if pT interpolation point is larger than PTCHANGE (currently 1.0 GeV)
		{
			double sign_of_f111 = (*loc_sign_spectra_ptr)[res_vector_indexer(parent_idx,npt-1,nphim1,npym1)];
			double sign_of_f121 = (*loc_sign_spectra_ptr)[res_vector_indexer(parent_idx,npt-1,nphi,npym1)];
			double sign_of_f211 = (*loc_sign_spectra_ptr)[res_vector_indexer(parent_idx,npt,nphim1,npym1)];
			double sign_of_f221 = (*loc_sign_spectra_ptr)[res_vector_indexer(parent_idx,npt,nphi,npym1)];
			double sign_of_f112 = (*loc_sign_spectra_ptr)[res_vector_indexer(parent_idx,npt-1,nphim1,npy)];
			double sign_of_f122 = (*loc_sign_spectra_ptr)[res_vector_indexer(parent_idx,npt-1,nphi,npy)];
			double sign_of_f212 = (*loc_sign_spectra_ptr)[res_vector_indexer(parent_idx,npt,nphim1,npy)];
			double sign_of_f222 = (*loc_sign_spectra_ptr)[res_vector_indexer(parent_idx,npt,nphi,npy)];
		
			//*******************************************************************************************************************
			// set f11
			//*******************************************************************************************************************
			// if using extrapolation and spectra at pT1 has larger magnitude than at pT0 (or the signs are different), just return zero
			if ( ptr > pT1 && ( log_f211 > log_f111 || sign_of_f111 * sign_of_f211 < 0 ) )
				f11 = 0.0;
			else if (sign_of_f111 * sign_of_f211 > 0)	// if the two points have the same sign in the pT direction, interpolate logs
				f11 = sign_of_f111 * exp( lin_int(ptr-pT0, one_by_pTdiff, log_f111, log_f211) );
			else					// otherwise, just interpolate original vals
				f11 = lin_int(ptr-pT0, one_by_pTdiff, f111, f211);
			
			//*******************************************************************************************************************
			// set f21
			//*******************************************************************************************************************
			if ( ptr > pT1 && ( log_f221 > log_f121 || sign_of_f121 * sign_of_f221 < 0 ) )
				f21 = 0.0;
			else if (sign_of_f121 * sign_of_f221 > 0)	// if the two points have the same sign in the pT direction, interpolate logs
				f21 = sign_of_f121 * exp( lin_int(ptr-pT0, one_by_pTdiff, log_f121, log_f221) );
			else					// otherwise, just interpolate original vals
				f21 = lin_int(ptr-pT0, one_by_pTdiff, f121, f221);
			//*******************************************************************************************************************
		
			//*******************************************************************************************************************
			// set f12
			//*******************************************************************************************************************
			// if using extrapolation and spectra at pT1 has larger magnitude than at pT0 (or the signs are different), just return zero
			if ( ptr > pT1 && ( log_f212 > log_f112 || sign_of_f112 * sign_of_f212 < 0 ) )
				f12 = 0.0;
			else if (sign_of_f112 * sign_of_f212 > 0)	// if the two points have the same sign in the pT direction, interpolate logs
				f12 = sign_of_f112 * exp( lin_int(ptr-pT0, one_by_pTdiff, log_f112, log_f212) );
			else					// otherwise, just interpolate original vals
				f12 = lin_int(ptr-pT0, one_by_pTdiff, f112, f212);
			
			//*******************************************************************************************************************
			// set f22
			//*******************************************************************************************************************
			if ( ptr > pT1 && ( log_f222 > log_f122 || sign_of_f122 * sign_of_f222 < 0 ) )
				f22 = 0.0;
			else if (sign_of_f122 * sign_of_f222 > 0)	// if the two points have the same sign in the pT direction, interpolate logs
				f22 = sign_of_f122 * exp( lin_int(ptr-pT0, one_by_pTdiff, log_f122, log_f222) );
			else					// otherwise, just interpolate original vals
				f22 = lin_int(ptr-pT0, one_by_pTdiff, f122, f222);
			//*******************************************************************************************************************
		}
		else						// if pT is smaller than PTCHANGE, just use linear interpolation, no matter what
		{
			f11 = lin_int(ptr-pT0, one_by_pTdiff, f111, f211);
			f21 = lin_int(ptr-pT0, one_by_pTdiff, f121, f221);
			f12 = lin_int(ptr-pT0, one_by_pTdiff, f112, f212);
			f22 = lin_int(ptr-pT0, one_by_pTdiff, f122, f222);
		}
				
		// now, interpolate f1 and f2 over the pphi direction
		double f1 = lin_int(pphir-phi0, one_by_pphidiff, f11, f21);
		double f2 = lin_int(pphir-phi0, one_by_pphidiff, f12, f22);

		//interpolate over pY last
		*result += lin_int(pyr-py0, one_by_pYdiff, f1, f2);

		return;
	}

	int load_decay_channel_info(
			vector<readindata::particle_info> all_particles,
			vector<int> & chosen_resonances )
	{
		int Nparticle = (int)all_particles.size();

		//initialize current resonance info
		current_resonance_mass = 0.0;
		current_resonance_mu = 0.0;
		current_resonance_Gamma = 0.0;
		current_resonance_total_br = 0.0;
		current_resonance_decay_masses.resize( 2 );
		current_resonance_decay_masses[0] = 0.0;
		current_resonance_decay_masses[1] = 0.0;
	
		int n_decay_channels = -1;

		//set all decay channel info
		if (chosen_resonances.size() == 0)
		{
			n_decay_channels = 1;
			n_resonance = (int)chosen_resonances.size();
			thermal_pions_only = true;
			if (VERBOSE > 0) cout << "Thermal pion(+) only!" << endl;
			decay_channels.resize( n_decay_channels );
			(decay_channels[0].resonance_decay_masses).resize( Maxdecaypart );	// Maxdecaypart == 5
		}
		else
		{
			//n_decay_channels is actually total number of decay channels which can generate pions
			//from chosen decay_channels
			n_decay_channels = readindata::get_number_of_decay_channels(chosen_resonances);
			n_resonance = (int)chosen_resonances.size();

			if (VERBOSE > 0) cout << "Computed n_decay_channels = " << n_decay_channels << endl
								<< "Computed n_resonance = " << n_resonance << endl;

			//decay_channels = new decay_info [n_decay_channels];
			decay_channels.resize( n_decay_channels );
			int temp_idx = 0;
			for (int icr = 0; icr < n_resonance; icr++)
			{
				readindata::particle_info particle_temp = all_particles[chosen_resonances[icr]];
				if (VERBOSE > 0) cout << "Loading resonance: " << particle_temp.name
						<< ", chosen_resonances[" << icr << "] = " << chosen_resonances[icr] << endl;

				for (int idecay = 0; idecay < particle_temp.decays; idecay++)
				{
					if (VERBOSE > 0) cout << "Current temp_idx = " << temp_idx << endl;
					if (temp_idx == n_decay_channels)	//i.e., all contributing decay channels have been loaded
						break;
					decay_channels[temp_idx].resonance_name = particle_temp.name;		// set name of resonance

					//check if effective branching is too small for inclusion in source variances
					bool effective_br_is_too_small = false;
					if (particle_temp.decays_effective_branchratio[idecay] <= 1.e-12)
						effective_br_is_too_small = true;

					decay_channels[temp_idx].resonance_particle_id = chosen_resonances[icr];	// set index of resonance in all_particles
					decay_channels[temp_idx].resonance_idx = icr;					// set index of resonance in chosen_resonances
					(decay_channels[temp_idx].resonance_decay_masses).resize( Maxdecaypart );	// Maxdecaypart == 5
					(decay_channels[temp_idx].resonance_decay_monvals).resize( Maxdecaypart );	// Maxdecaypart == 5
					(decay_channels[temp_idx].resonance_decay_Gammas).resize( Maxdecaypart );	// Maxdecaypart == 5

					//Set resonance decay masses
					for (int ii = 0; ii < Maxdecaypart; ii++)
					{
						decay_channels[temp_idx].resonance_decay_monvals[ii] = particle_temp.decays_part[idecay][ii];
						if (particle_temp.decays_part[idecay][ii] == 0)
						{
							decay_channels[temp_idx].resonance_decay_masses[ii] = 0.0;
							decay_channels[temp_idx].resonance_decay_Gammas[ii] = 0.0;

						}
						else
						{
							int tempID = readindata::lookup_particle_id_from_monval(
											particle_temp.decays_part[idecay][ii]);

							decay_channels[temp_idx].resonance_decay_masses[ii] = all_particles[tempID].mass;
							decay_channels[temp_idx].resonance_decay_Gammas[ii] = all_particles[tempID].width;

						}
					}
					decay_channels[temp_idx].resonance_mu = particle_temp.mu;
					decay_channels[temp_idx].resonance_gspin = particle_temp.gspin;
					decay_channels[temp_idx].resonance_sign = particle_temp.sign;
					decay_channels[temp_idx].resonance_mass = particle_temp.mass;
					decay_channels[temp_idx].nbody = abs(particle_temp.decays_Npart[idecay]);
					decay_channels[temp_idx].resonance_Gamma = particle_temp.width;
					decay_channels[temp_idx].resonance_total_br = particle_temp.decays_effective_branchratio[idecay];
					decay_channels[temp_idx].resonance_direct_br = particle_temp.decays_branchratio[idecay];
				
					bool lifetime_is_too_long = false;

					// if decay channel parent resonance is not too long-lived
					// and decay channel contains at least one target daughter particle,
					// include channel
					decay_channels[temp_idx].include_channel = !lifetime_is_too_long && !effective_br_is_too_small;

					if (VERBOSE > 0) cout
							<< "Resonance = " << decay_channels[temp_idx].resonance_name
							<< ", decay channel " << idecay + 1
							<< ": mu=" << decay_channels[temp_idx].resonance_mu
							<< ", gs=" << decay_channels[temp_idx].resonance_gspin
							<< ", sign=" << decay_channels[temp_idx].resonance_sign
							<< ", M=" << decay_channels[temp_idx].resonance_mass
							<< ", nbody=" << decay_channels[temp_idx].nbody
							<< ", Gamma=" << decay_channels[temp_idx].resonance_Gamma
							<< ", total br=" << decay_channels[temp_idx].resonance_total_br
							<< ", direct br=" << decay_channels[temp_idx].resonance_direct_br << endl;

					if (VERBOSE > 0) cout
							<< "Resonance = " << decay_channels[temp_idx].resonance_name << ": ";
					for (int decay_part_idx = 0; decay_part_idx < decay_channels[temp_idx].nbody; decay_part_idx++)
						if (VERBOSE > 0) cout
							<< "m[" << decay_part_idx << "] = "
							<< decay_channels[temp_idx].resonance_decay_masses[decay_part_idx] << "   "
							<< decay_channels[temp_idx].resonance_decay_monvals[decay_part_idx] << "   ";
					if (VERBOSE > 0) cout << endl << endl;

					temp_idx++;
				}
			}
		}

		return (n_decay_channels);
	}

	void Get_current_decay_string(int dc_idx, string * decay_string)
	{
		// N.B. - dc_idx == 0 is thermal pion(+)s in calling loop, dc_idx > 0 gives resonance decays
		//      ==> need to use dc_idx - 1 here
		*decay_string = "\t" + decay_channels[dc_idx - 1].resonance_name + "\t\t--->> ";
		int temp_monval, tempID;
		for (int decay_part_idx = 0; decay_part_idx < decay_channels[dc_idx - 1].nbody; decay_part_idx++)
		{
			temp_monval = decay_channels[dc_idx - 1].resonance_decay_monvals[decay_part_idx];
			//if (VERBOSE > 0) cout << "Get_current_decay_string(): temp_monval = " << temp_monval << endl;
			if (temp_monval == 0)
				continue;
			else
			{
				tempID = readindata::lookup_particle_id_from_monval(temp_monval);
				*decay_string += all_particles[tempID].name;
				if (decay_part_idx < decay_channels[dc_idx - 1].nbody - 1) *decay_string += " + ";
			}
		}
		return;
	}

	bool Do_this_decay_channel(int dc_idx)
	{
		if (dc_idx == 0)
		{
			return true;
		}

		string local_name = decay_channels[dc_idx-1].resonance_name;

		Get_current_decay_string(dc_idx, &current_decay_channel_string);

		bool tmp_bool = decay_channels[dc_idx-1].include_channel;
		if (!tmp_bool && VERBOSE > 1)
			cout << endl << current_decay_channel_string << " (SKIPPED)" << endl;

		return (tmp_bool);
	}

	// ************************************************************
	// Checks whether to do daughter particle for any given decay channel
	// ************************************************************
	bool Do_this_daughter_particle(int dc_idx, int daughter_idx, int * daughter_resonance_pid)
	{
		int Nparticle = (int)all_particles.size();
		// assume dc_idx > 0
		string local_name = decay_channels[dc_idx-1].resonance_name;

		// look up daughter particle info
		int temp_monval = decay_channels[dc_idx-1].resonance_decay_monvals[daughter_idx];

		if (temp_monval == 0)
			return false;

		int temp_ID = readindata::lookup_particle_id_from_monval(temp_monval);

		*daughter_resonance_pid = temp_ID;

		// if daughter was found in chosen_resonances or is pion(+),
		// then this is correct
		readindata::particle_info temp_daughter = all_particles[temp_ID];

		//if can't find daughter pid, we have problems...
		if ( *daughter_resonance_pid < 0
				&& temp_daughter.monval != particle_monval
				&& temp_daughter.effective_branchratio >= 1.e-12 )
		{
			cout
				<< "Couldn't find " << temp_daughter.name
				<< " in chosen_resonances!  Results are probably not reliable..."
				<< endl;
			exit (1);
		}

		//check if daughter contributes to pions
		bool daughter_does_not_contribute
				= ( (temp_daughter.decays_Npart[0] == 1
						|| temp_daughter.effective_branchratio < 1.e-12)
					&& temp_daughter.monval != particle_monval );

		return (daughter_does_not_contribute);
	}

	void Set_current_particle_info(int dc_idx)
	{
		if (dc_idx == 0)
		{
			current_resonance_pid = target_pid;
		
			return;
		}
		else
		{
			// assume dc_idx > 0
			string local_name = decay_channels[dc_idx-1].resonance_name;

			if (VERBOSE > 0)
				cout << current_decay_channel_string << endl;

			//cerr << "Setting current decay channel information for dc_idx = " << dc_idx << endl;
			current_decay_channel_idx = dc_idx;
			current_resonance_idx = decay_channels[dc_idx-1].resonance_idx;
			current_resonance_pid = decay_channels[dc_idx-1].resonance_particle_id;
			current_resonance_mass = decay_channels[dc_idx-1].resonance_mass;
			current_resonance_Gamma = decay_channels[dc_idx-1].resonance_Gamma;
			current_resonance_total_br = decay_channels[dc_idx-1].resonance_total_br;
			current_resonance_direct_br = decay_channels[dc_idx-1].resonance_direct_br;
			current_reso_nbody = decay_channels[dc_idx-1].nbody;
		}
	
		return;
	}

	void Set_current_daughter_info(int dc_idx, int daughter_idx)
	{
		current_daughter_mass = decay_channels[dc_idx-1].resonance_decay_masses[daughter_idx];
		current_daughter_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[daughter_idx];

		// set non-daughter decay masses for computing contributions to spectra of daughter
		double m2ex = 0.0, m3ex = 0.0, m4ex = 0.0;
		switch(current_reso_nbody)
		{
			case 1:
				break;
			case 2:
				current_resonance_decay_masses[1] = 0.0;
				current_m3_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[0];
				if (daughter_idx == 0)
				{
					current_resonance_decay_masses[0] = decay_channels[dc_idx-1].resonance_decay_masses[1];
					current_m2_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[1];
				}
				else
				{
					current_resonance_decay_masses[0] = decay_channels[dc_idx-1].resonance_decay_masses[0];
					current_m2_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[0];
				}
				break;
			case 3:
				if (daughter_idx == 0)
				{
					current_resonance_decay_masses[0] = decay_channels[dc_idx-1].resonance_decay_masses[1];
					current_resonance_decay_masses[1] = decay_channels[dc_idx-1].resonance_decay_masses[2];
					current_m2_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[1];
					current_m3_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[2];
				}
				else if (daughter_idx == 1)
				{
					current_resonance_decay_masses[0] = decay_channels[dc_idx-1].resonance_decay_masses[0];
					current_resonance_decay_masses[1] = decay_channels[dc_idx-1].resonance_decay_masses[2];
					current_m2_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[0];
					current_m3_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[2];
				}
				else
				{
					current_resonance_decay_masses[0] = decay_channels[dc_idx-1].resonance_decay_masses[0];
					current_resonance_decay_masses[1] = decay_channels[dc_idx-1].resonance_decay_masses[1];
					current_m2_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[0];
					current_m3_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[1];
				}
				break;
			case 4:
				if (daughter_idx == 0)
				{
					m2ex = decay_channels[dc_idx-1].resonance_decay_masses[1];
					m3ex = decay_channels[dc_idx-1].resonance_decay_masses[2];
					m4ex = decay_channels[dc_idx-1].resonance_decay_masses[3];
				}
				else if (daughter_idx == 1)
				{
					m2ex = decay_channels[dc_idx-1].resonance_decay_masses[0];
					m3ex = decay_channels[dc_idx-1].resonance_decay_masses[2];
					m4ex = decay_channels[dc_idx-1].resonance_decay_masses[3];
				}
				else if (daughter_idx == 2)
				{
					m2ex = decay_channels[dc_idx-1].resonance_decay_masses[0];
					m3ex = decay_channels[dc_idx-1].resonance_decay_masses[1];
					m4ex = decay_channels[dc_idx-1].resonance_decay_masses[3];
				}
				else
				{
					m2ex = decay_channels[dc_idx-1].resonance_decay_masses[0];
					m3ex = decay_channels[dc_idx-1].resonance_decay_masses[1];
					m4ex = decay_channels[dc_idx-1].resonance_decay_masses[2];
				}
				current_resonance_decay_masses[0] = m2ex;
				current_resonance_decay_masses[1] = 0.5 * (m3ex + m4ex + current_resonance_mass - current_daughter_mass - m2ex);
				break;
			default:
				cerr << "Set_current_daughter_info(): shouldn't have ended up here, bad value of current_reso_nbody!" << endl;
				exit(1);
		}
	}

	void Allocate_decay_channel_info()
	{
		//if (VERBOSE > 2) cout << "Reallocating memory for decay channel information..." << endl;

		//set resonance phase-space integral points
		v_pts = vector<double>(n_v_pts);
		v_wts = vector<double>(n_v_pts);
		gauss_quadrature(n_v_pts, 1, 0.0, 0.0, vmin, vmax, v_pts, v_wts);

		zeta_pts = vector<double>(n_zeta_pts);
		zeta_wts = vector<double>(n_zeta_pts);
		gauss_quadrature(n_zeta_pts, 1, 0.0, 0.0, zetamin, zetamax, zeta_pts, zeta_wts);

		s_pts = vector<double>(n_s_pts);
		s_wts = vector<double>(n_s_pts);
		x_pts = vector<double>(n_s_pts);
		x_wts = vector<double>(n_s_pts);
		gauss_quadrature(n_s_pts, 1, 0.0, 0.0, -1.0, 1.0, x_pts, x_wts);

		//set//allocate other stuff
		VEC_n2_v_factor.resize(n_v_pts);
		VEC_n2_P_Y.resize(n_v_pts);
		VEC_n2_zeta_factor.resize(n_v_pts * n_zeta_pts);
		VEC_n2_PPhi_tilde.resize(n_v_pts * n_zeta_pts);
		VEC_n2_PPhi_tildeFLIP.resize(n_v_pts * n_zeta_pts);
		VEC_n2_PT.resize(n_v_pts * n_zeta_pts);
		VEC_n2_Ppm.resize(n_v_pts * n_zeta_pts * 2);
		for(int ii = 0; ii < n_v_pts * n_zeta_pts * 2; ++ii)
			VEC_n2_Ppm[ii].resize(4);	//four corresponds to space-time components

		//VEC_n3_g_s.resize(n_s_pts);
		VEC_n3_s_factor.resize(n_s_pts);
		VEC_n3_v_factor.resize(n_s_pts * n_v_pts);
		VEC_n3_zeta_factor.resize(n_s_pts * n_v_pts * n_zeta_pts);
		VEC_n3_P_Y.resize(n_s_pts * n_v_pts);
		VEC_n3_PPhi_tilde.resize(n_s_pts * n_v_pts * n_zeta_pts);
		VEC_n3_PPhi_tildeFLIP.resize(n_s_pts * n_v_pts * n_zeta_pts);
		VEC_n3_PT.resize(n_s_pts * n_v_pts * n_zeta_pts);
		VEC_n3_Ppm.resize(n_s_pts * n_v_pts * n_zeta_pts * 2);
		for(int ii = 0; ii < n_s_pts * n_v_pts * n_zeta_pts * 2; ++ii)
			VEC_n3_Ppm[ii].resize(4);	//four corresponds to space-time components

		//if (VERBOSE > 2) cout << "Reallocated memory for decay channel information." << endl;

		return;
	}

	int lookup_resonance_idx_from_particle_id(int pid)
	{
		// pid - particle index in all_particles array
		// looks up location in chosen_resonance_indices of given value particle_id
		int result = -1;

		for (int ii = 0; ii < (int)chosen_resonance_indices.size(); ii++)
		{
			if (chosen_resonance_indices[ii] == pid)
			{
				result = ii;
				break;
			}
		}

		// if pid is not one of the chosen_resonance_indices,
		// is not the target daughter (pion(+)),
		// is not stable and has a non-zero effective branching ratio
		if (result < 0 && pid != target_pid && all_particles[pid].stable == 0 && all_particles[pid].effective_branchratio >= 1.e-12)
		{
			cout << " *** lookup_resonance_idx_from_particle_id(): Particle_id = " << pid
					<< " (" << all_particles[pid].name
					<<" ) not found in chosen_resonance_indices!" << endl
					<< " *** br = " << all_particles[pid].effective_branchratio << endl
					<< " *** Can only choose from: " << endl;
			for (int ii = 0; ii < (int)chosen_resonance_indices.size(); ii++)
				cout << all_particles[chosen_resonance_indices[ii]].name
						<< ": pid = " << chosen_resonance_indices[ii] << endl;

		}
		return (result);
	}

}

//End of file
