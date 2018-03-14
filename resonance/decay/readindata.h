#ifndef READINDATA_H
#define READINDATA_H

#include <string>
#include <fstream>
#include <ctime>
#include <vector>
#include <algorithm>
#include <iostream>

#include <limits.h>

#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

#include "parameters.h"

using namespace std;
using namespace parameters;

namespace readindata
{
	typedef struct
	{
		int monval;			// Montecarlo number according PDG
		string name;
		double mass;
		double width;
		int gspin;			// spin degeneracy
		int baryon;
		int strange;
		int charm;
		int bottom;
		int gisospin;			// isospin degeneracy
		int charge;
		int decays;			// amount of decays listed for this resonance
		int stable;			// defines whether this particle is considered as stable
		int decays_Npart[Maxdecaychannel];
		double decays_branchratio[Maxdecaychannel];
		int decays_part[Maxdecaychannel][Maxdecaypart];
		double mu;
		double thermal_yield;
		double percent_contribution;	//used to compute approximate percentage contribution of particle to net yield of some daughter particle
		double effective_branchratio;	//the effective branching ratio of generic resonance to some specified daughter particle
					//N.B. - effective branching ratio may be larger than 1
		double decays_effective_branchratio[Maxdecaychannel];
					//similar to effective_branchratio, but specific to each decay channel
		int sign;       //Bose-Einstein or Dirac-Fermi statistics
	} particle_info;

	////////////////////////////////////////////////////////
	int get_filelength( string filepath );
	int get_filewidth( string filepath );
	void read_resonance( vector<particle_info> & all_particles );

	void compute_total_contribution_percentages( int stable_particle_idx );

	void calculate_thermal_particle_yield( double Temperature );

	double b_j_to_i( int j, int i, int verbose_monval = 0 );

	int count_targets( int * decay_channel_particles,
						particle_info * i );

	int count_stable( int * decay_channel_particles );

	bool is_stable( int monval );

	int set_stable_particle_monval();

	int lookup_particle_id_from_monval( int monval );

	void print_particle_stability();

	int get_number_of_decay_channels(vector<int> & chosen_resonances);

	void get_important_resonances(
			int chosen_target_particle_idx,
			vector<int> & chosen_resonance_indices,
			double threshold,
			double &running_total_percentage );

	void get_all_descendants(
			vector<int> & chosen_resonance_indices );

	void sort_by_mass(
			vector<int> & chosen_resonance_indices );

	void load_resonance_info(
			vector<particle_info> & all_particles_in,
			vector<int> & chosen_resonance_indices_ptr,
			double Tfo, int particle_idx );

	////////////////////////////////////////////////////////
	//Stuff I have to define in the header file...
	template <typename T>
	vector<size_t> ordered(vector<T> const& values, int lt_or_gt = 0)
	{
		using namespace boost::phoenix;
		using namespace boost::phoenix::arg_names;

		vector<size_t> indices(values.size());
		int i = 0;
		transform(values.begin(), values.end(), indices.begin(), ref(i)++);
		if (lt_or_gt == 0)
			sort(indices.begin(), indices.end(), ref(values)[arg1] < ref(values)[arg2]);
		else
			sort(indices.begin(), indices.end(), ref(values)[arg1] > ref(values)[arg2]);
		return indices;
	}
}

#endif
