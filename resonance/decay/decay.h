#ifndef DECAY_H
#define DECAY_H

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

#include "../gauss/gauss_quadrature.h"
#include "parameters.h"
#include "readindata.h"

namespace decay
{
	typedef struct
	{
		int resonance_particle_id;		// keeps track of current resonance's index in all_particles array
		int resonance_idx;				// keeps track of current resonance's index in chosen_resonances vector
		int nbody;
		int resonance_sign;
		double resonance_mass;
		double resonance_mu;
		double resonance_gspin;
		double resonance_Gamma;
		double resonance_total_br;
		double resonance_direct_br;
		vector<double> resonance_decay_masses;
		vector<int> resonance_decay_monvals;
		vector<double> resonance_decay_Gammas;
		string resonance_name;
		bool include_channel;
	} decay_info;

	//load all decay info
	int load_decay_channel_info(
			vector<readindata::particle_info> all_particles,
			vector<int> & chosen_resonances );

	//driver functions
	void Compute_phase_space_integrals( vector<readindata::particle_info> all_particles,
										vector<int> chosen_resonance_indices,
										vector<double> * spectra_re_in,
										vector<double> * spectra_im_in,
										int target_particle_id_in,
										bool assume_azimuthal_symmetry = false );
	void Do_resonance_integrals( int decay_channel,
									int parent_pid,
									int daughter_pid,
									vector<double> * spec_re,
									vector<double> * spec_im,
									bool assume_azimuthal_symmetry );

	//decay routines
	inline int NB2_indexer(const int iv, const int izeta);
	inline int NB3_indexer(const int is, const int iv, const int izeta);
	int lookup_resonance_idx_from_particle_id(int pid);
	double get_Q();
	double g(double s);

	//load decay info
	void Load_decay_channel_info_nb2(int dc_idx, double K_T_local, double K_phi_local, double K_y_local);
	void Load_decay_channel_info_nb3(int dc_idx, double K_T_local, double K_phi_local, double K_y_local);

	//help with decay logic
	bool Do_this_decay_channel(int dc_idx);
	bool Do_this_daughter_particle(int dc_idx, int daughter_idx, int * daughter_resonance_pid);
	void Set_current_particle_info(int dc_idx);
	void Set_current_daughter_info(int dc_idx, int daughter_idx);

	//decay routines
	void two_body_decay( int dc_idx, int parent_idx, int daughter_idx,
							vector<double> * spec_re, vector<double> * spec_im,
							bool assume_azimuthal_symmetry );
	void three_body_decay( int dc_idx, int parent_idx, int daughter_idx,
							vector<double> * spec_re, vector<double> * spec_im,
							bool assume_azimuthal_symmetry );

	//spectra interpolator
	void Edndp3(double ptr, double pphir, double pyr, int parent_idx, double * result,
				vector<double> * loc_spectra_ptr, vector<double> * loc_log_spectra_ptr,
				vector<double> * loc_sign_spectra_ptr);
	void Edndp3(double ptr, double pyr, int parent_idx, double * result,
				vector<double> * loc_spectra_ptr, vector<double> * loc_log_spectra_ptr,
				vector<double> * loc_sign_spectra_ptr);

	//memory management
	void Allocate_decay_channel_info();

	//misc
	inline void reset(vector<double> * v, double val = 0.0)
	{
		fill(v->begin(), v->end(), val);
	}

	//stuff I have to define in a header file
	template <typename T> int sgn(T val)
	{
		return (T(0) < val) - (val < T(0));
	}

	inline int mom_indexer(const int ipT, const int ipphi, const int ipY)
	{
		return (
			( ipT * n_pphi_pts + ipphi ) * n_pY_pts + ipY
		);
	}

	inline int res_vector_indexer(const int iRes, const int ipT, const int ipphi, const int ipY)
	{
		return (
			( ( iRes * n_pT_pts + ipT ) * n_pphi_pts + ipphi ) * n_pY_pts + ipY
		);
	}

	inline int NB2_indexer(const int iv, const int izeta)
	{
		return (
			iv * n_zeta_pts + izeta
		);
	}

	inline int NB3_indexer(const int is, const int iv, const int izeta)
	{
		return (
			( is * n_v_pts + iv ) * n_zeta_pts + izeta
		);
	}

	inline double place_in_range(double phi, double min, double max)
	{
		while (phi < min || phi > max)
		{
			if (phi < min) phi += 2.0*M_PI;
			else phi -= 2.0*M_PI;
		}

		return (phi);
	}

	inline double lin_int(double x_m_x1, double one_by_x2_m_x1, double f1, double f2)
	{
		return ( f1 + (f2 - f1) * x_m_x1 * one_by_x2_m_x1 );
	}

}


#endif
