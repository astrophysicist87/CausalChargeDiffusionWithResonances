#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<cstdlib>
#include<gsl/gsl_sf_bessel.h>
#include<vector>
#include<algorithm>

#include "readindata.h"

using namespace std;

extern vector<double> pT_pts, pT_wts;
extern vector<double> pphi_pts, pphi_wts;
extern vector<double> Del_pY_pts, Del_pY_wts;

namespace readindata
{
	int n_resonances = -1;
	double * stable_particle_monval;
	int Nstable_particle;
	double ** all_b_j_to_i;
	bool output_mu_and_yield = false;
	const int chosenParticlesMode = 1;

	vector<double> resonance_spectra_RE;
	vector<double> resonance_spectra_IM;

	vector<particle_info> all_particles;

	void load_resonance_info(
			vector<particle_info> & all_particles_in,
			vector<int> & chosen_resonance_indices,
			double Tfo, int particle_idx)
	{
		cout << "Inside load_resonance_info()..." << endl;
		//set and modify particle_info vector
		//while loading resonance information
		read_resonance(all_particles_in);

		//avoid unncessary copies from here on out
		all_particles = all_particles_in;

		calculate_thermal_particle_yield(Tfo);

		compute_total_contribution_percentages(particle_idx);

		// Get chosen particles
		double threshold = 0.0;
		double net_fraction_resonance_contribution = 0.0;

		if (chosenParticlesMode == 0)				// calculate chosen resonances from threshold
		{
			threshold = resonanceThreshold;

			get_important_resonances(particle_idx, chosen_resonance_indices,
										threshold, net_fraction_resonance_contribution);

			get_all_descendants(chosen_resonance_indices);

			//chosen_resonance_indices_ptr.push_back(particle_idx);

			sort_by_mass(chosen_resonance_indices);
		}
		else if (chosenParticlesMode == 1)		// read chosen resonances in from file
		{
			ifstream chosenParticlesStream("./decay/EOS/chosen_particles.dat");
			int len = 0;
			long tmp;
			while (chosenParticlesStream >> tmp)
			{
				chosen_resonance_indices.push_back(lookup_particle_id_from_monval(tmp));
				len++;
			}

			chosenParticlesStream.close();

			sort_by_mass(chosen_resonance_indices);

			//chosen_resonance_indices.pop_back();

			net_fraction_resonance_contribution = 1.0;
		}
		else																// otherwise, you did something wrong!
		{
			cerr << "chosenParticlesMode = " << chosenParticlesMode << " not supported!  Exiting..." << endl;
			exit(1);
		}

		//pass loaded array back to calling routine
		all_particles_in = all_particles;

		return;
	}


	int get_filelength(string filepath)
	{
	   int length=0; 
	   char line[512];
	   ostringstream filepath_stream;
	   filepath_stream << filepath;
	   ifstream infile(filepath_stream.str().c_str());
	   //Determine the length of the file
	   while (!infile.eof ())
	   {
		  infile.getline(line, 512);
		  length++;
	   }
	   length = length-1;
	   infile.close();
	   return(length);
	}

	int get_filewidth(string filepath)
	{
		ostringstream filepath_stream;
		filepath_stream << filepath;
		ifstream infile(filepath_stream.str().c_str());
		string line, temp;
		stringstream ss;
		int ncols=0;
		getline(infile, line);
		ss.clear();
		ss << line;
		
		while (ss >> temp)
			ncols++;

		return ncols;
	}

	void read_resonance(vector<particle_info> & all_particles)
	{
	   int Nparticle=0; 
	   cout << "Reading in particle resonance decay table...";
	   ifstream resofile("./decay/EOS/pdg.dat");
	   int local_i = 0;
	   int dummy_int;
		all_particles.clear();
		all_particles.resize(Maxparticle);
	   while (!resofile.eof())
	   {
		  resofile >> all_particles[local_i].monval;
		  resofile >> all_particles[local_i].name;
		  resofile >> all_particles[local_i].mass;
		  resofile >> all_particles[local_i].width;
		  resofile >> all_particles[local_i].gspin;	      //spin degeneracy
		  resofile >> all_particles[local_i].baryon;
		  resofile >> all_particles[local_i].strange;
		  resofile >> all_particles[local_i].charm;
		  resofile >> all_particles[local_i].bottom;
		  resofile >> all_particles[local_i].gisospin;     //isospin degeneracy
		  resofile >> all_particles[local_i].charge;
		  resofile >> all_particles[local_i].decays;
		  for (int j = 0; j < all_particles[local_i].decays; j++)
		  {
		     resofile >> dummy_int;
		     resofile >> all_particles[local_i].decays_Npart[j];
		     resofile >> all_particles[local_i].decays_branchratio[j];
		     resofile >> all_particles[local_i].decays_part[j][0];
		     resofile >> all_particles[local_i].decays_part[j][1];
		     resofile >> all_particles[local_i].decays_part[j][2];
		     resofile >> all_particles[local_i].decays_part[j][3];
		     resofile >> all_particles[local_i].decays_part[j][4];
		  }
		  
		  //decide whether particle is stable under strong interactions
		  if(all_particles[local_i].decays_Npart[0] == 1) 	      
		     all_particles[local_i].stable = 1;
		  else
		     all_particles[local_i].stable = 0;

		  //add anti-particle entry
		  if(all_particles[local_i].baryon == 1)
		  {
		     local_i++;
		     all_particles[local_i].monval = -all_particles[local_i-1].monval;
		     ostringstream antiname;
		     antiname << "Anti-" << all_particles[local_i-1].name;
		     all_particles[local_i].name = antiname.str();
		     all_particles[local_i].mass = all_particles[local_i-1].mass;
		     all_particles[local_i].width = all_particles[local_i-1].width;
		     all_particles[local_i].gspin = all_particles[local_i-1].gspin;
		     all_particles[local_i].baryon = -all_particles[local_i-1].baryon;
		     all_particles[local_i].strange = -all_particles[local_i-1].strange;
		     all_particles[local_i].charm = -all_particles[local_i-1].charm;
		     all_particles[local_i].bottom = -all_particles[local_i-1].bottom;
		     all_particles[local_i].gisospin = all_particles[local_i-1].gisospin;
		     all_particles[local_i].charge = -all_particles[local_i-1].charge;
		     all_particles[local_i].decays = all_particles[local_i-1].decays;
		     all_particles[local_i].stable = all_particles[local_i-1].stable;
		     for (int j = 0; j < all_particles[local_i].decays; j++)
		     {
		        all_particles[local_i].decays_Npart[j]=all_particles[local_i-1].decays_Npart[j];
		        all_particles[local_i].decays_branchratio[j]=all_particles[local_i-1].decays_branchratio[j];
		        //for (int k=0; k< Maxdecaypart; k++)						//commented out by Chris Plumberg - 06/08/2015
		        //   all_particles[local_i].decays_part[j][k]=all_particles[local_i-1].decays_part[j][k];
		        for (int k=0; k< Maxdecaypart; k++)							//replaced with following k-loop
		        {
		           int idx = 0;  
		           for(int ii=0; ii < local_i; ii++) // find the index for decay particle
		           {
		              if(all_particles[local_i-1].decays_part[j][k] == all_particles[ii].monval)
		              {
		                 idx = ii;
		                 break;
		              }
		           }
		           if(idx == local_i-1 && all_particles[local_i-1].stable == 0)  // check
		           {
		              cout << "Error: can not find decay all_particles index for anti-baryon!" << endl;
		              cout << "all_particles monval : " << all_particles[local_i-1].decays_part[j][k] << endl;
		              exit(1);
		           }
		           if(all_particles[idx].baryon == 0 && all_particles[idx].charge == 0 && all_particles[idx].strange == 0)
		              all_particles[local_i].decays_part[j][k]= all_particles[local_i-1].decays_part[j][k];
		           else
		              all_particles[local_i].decays_part[j][k]= -all_particles[local_i-1].decays_part[j][k];
		        } // end of k-loop
		     }
		   }
		   local_i++;	// Add one to the counting variable "i" for the meson/baryon
	   }
	   resofile.close();
	   Nparticle=local_i-1; //take account the final fake one
	   for(int i=0; i < Nparticle; i++)
	   {
		  if(all_particles[i].baryon==0)
		     all_particles[i].sign=-1;
		  else
		     all_particles[i].sign=1;
	   }
	   cout << "done! Antiparticles are added!" << endl;
		all_particles.resize(Nparticle);
	   return;
	}

	void calculate_thermal_particle_yield(double Temperature)
	{
		int Nparticle = (int)all_particles.size();
		double one_by_Tconv = 1./Temperature;
		//double * all_particle_fugacities = new double [Maxparticle];
		//set_to_zero(all_particle_thermal, Nparticle);
		for (int j = 0; j < Nparticle; j++)
		{
			double yield = 0.0, fugacity = 0.0;
			double gj = all_particles[j].gspin;
			double mj = all_particles[j].mass;
			if (mj == 0)
			{
				fugacity = 0.0;
				all_particles[j].thermal_yield = 0.0;
				continue;
			}
			double pm = -all_particles[j].sign;
			fugacity = exp(one_by_Tconv * all_particles[j].mu);
			for (int k = 1; k <= 10; k++)
				yield += double(pow(pm, k+1))*pow(fugacity, (double)k)*gsl_sf_bessel_Kn(2, double(k)*mj*one_by_Tconv)/double(k);
			all_particles[j].thermal_yield = yield*gj*mj*mj/(2.*M_PI*M_PI);
		}
		//**********************************************************************
		if (output_mu_and_yield)
		{
			ostringstream output_stream;
			output_stream << "check_mu_and_yield.dat";
			ofstream output(output_stream.str().c_str());
			for (int j = 0; j < Nparticle; j++)
				output << all_particles[j].name << "   " << all_particles[j].mu << "   " << all_particles[j].thermal_yield << endl;
			output.close();
		}
		//**********************************************************************
	
		return;
	}

	void compute_total_contribution_percentages(int stable_particle_idx)
	{
		int Nparticle = (int)all_particles.size();
		double denominator = 0.0, temp;
		all_b_j_to_i = new double * [Nparticle];
		for (int i = 0; i < Nparticle; i++)
		{
			all_particles[i].percent_contribution = 0.0;
			all_b_j_to_i[i] = new double [Nparticle];
			for (int j = 0; j < Nparticle; j++)
				all_b_j_to_i[i][j] = 0.0;
		}
		for (int i = 0; i < Nparticle; i++)
		{
			temp = b_j_to_i(i, stable_particle_idx);
			all_particles[i].percent_contribution = temp * all_particles[i].thermal_yield;
			all_particles[i].effective_branchratio = temp;
			denominator += temp * all_particles[i].thermal_yield;
		}
		for (int i = 0; i < Nparticle; i++)
			all_particles[i].percent_contribution /= 0.01 * denominator;	//0.01 factor makes it a percentage
	
		return;
	}

	double b_j_to_i(int j, int i, int verbose_monval /* = 0*/)
	{
		int Nparticle = (int)all_particles.size();
		double result = 0.0;
		particle_info parent = all_particles[j];
		particle_info target = all_particles[i];
		int verbose = 0;
		if (parent.monval == verbose_monval) verbose = 1;
		if (verbose > 0) cout << "Currently looking at decay chains of " << parent.name << " to " << target.name << ":" << endl;
		if ( (parent.decays == 1) && (parent.decays_Npart[0] == 1) )
		{
			if (verbose > 0) cout << "   --> " << parent.name << " is stable, so moving on!" << endl;
			return(0.0);	// if parent is already stable, return 0
		}
		else
		{
			if (verbose > 0) cout << "   --> " << parent.name << " is unstable, so continuing!" << endl;
			//all_particles[j].stable == 0;	//just in case
		}
	
		// added this to try to save some time and simplify debugging output
		if (fabs(all_b_j_to_i[j][i]) > 1.e-12)
		{
			if (verbose > 0) cout << "   --> already calculated!  Recycling stored value and moving on!" << endl;
			if (verbose > 0) cout << "   --> " << parent.name << "-->" << target.name << ": using b_j_to_i = " << all_b_j_to_i[j][i] << endl;
			return all_b_j_to_i[j][i];
		}
	
		for (int k = 0; k < parent.decays; k++)
		{
			int nki = count_targets(parent.decays_part[k], &target);	// number of target decay particles in kth decay channel
			double bk = parent.decays_branchratio[k];			// branching ratio for kth decay channel
			int nks = count_stable(parent.decays_part[k]);			// number of stable decay particles in kth decay channel
			int nktot = abs(parent.decays_Npart[k]);				// total number of decay particles in kth decay channel
			if (verbose > 0) cout << " - " << parent.name << "(monval = " << parent.monval << "): decay channel " << k + 1 << " of " << parent.decays << endl;
			if (verbose > 0) cout << "   --> nki = " << nki << ", nks = " << nks << ", nktot = " << nktot << endl;
			if ((nki == 0) && (nks == nktot))
			{
				if (verbose > 0) cout << "   --> found no " << target.name << "s or unstable particles, so moving on!" << endl;
				continue;			// if kth decay channel contains no target particles or other particles which might decay to some
			}
			if (nki != 0) result += double(nki)*bk;				// if kth decay channel contains target particles
			parent.decays_effective_branchratio[k] = double(nki)*bk;
			all_particles[j].decays_effective_branchratio[k] = double(nki)*bk;
			if (nks != nktot)						// if there are unstable particles
			{
				if (verbose > 0) cout << "   --> found some unstable particles!" << endl;
				for (int ipart = 0; ipart < nktot; ipart++)
				{		// apply this same process recursively to all unstable daughter resonances
					if ( !is_stable( parent.decays_part[k][ipart] ) )
					{
						int decay_particle_idx = lookup_particle_id_from_monval(parent.decays_part[k][ipart]);
						if (verbose > 0)
							cout << "   --> now working on unstable particle ("
									<< all_particles[decay_particle_idx].name << ") with monval = "
									<< parent.decays_part[k][ipart] << endl << endl;
						double temp_bj2i = b_j_to_i(decay_particle_idx, i);
						if (verbose > 0)
							cout << "   --> " << parent.name << "-->" << target.name
									<< ": using b_j_to_i = " << temp_bj2i << endl;
						result += bk * temp_bj2i;
						parent.decays_effective_branchratio[k] += bk * temp_bj2i;
						all_particles[j].decays_effective_branchratio[k] += bk * temp_bj2i;
					}
				}
			}
		}
	
		all_b_j_to_i[j][i] = result;
		if (verbose > 0) cout << "***FINAL***: " << parent.name << "-->" << target.name << ": using b_j_to_i = " << all_b_j_to_i[j][i] << endl;
		return (result);
	}

	int lookup_particle_id_from_monval(int monval)
	{
		int Nparticle = (int)all_particles.size();
		for (int ipart = 0; ipart < Nparticle; ipart++)
		{
			//cout << "lookup(" << monval << "): " << all_particles[ipart].name << "   " << all_particles[ipart].monval << endl;
			if (monval == all_particles[ipart].monval) return ipart;
		}
		cerr << "monval = " << monval << endl;
		cerr << "Only available monvals are:" << endl;
		for (int ipart = 0; ipart < Nparticle; ipart++)
			cerr << all_particles[ipart].name << "   " << all_particles[ipart].monval << endl;
		cerr << "Could not find monval in PDG table!  Aborting..." << endl;
		exit(1);
	}

	int count_targets(int * decay_channel_particles, particle_info * i)
	{
		int count = 0;
		for (int idcp = 0; idcp < Maxdecaypart; idcp++)
			if (i->monval == decay_channel_particles[idcp]) count++;
		return(count);
	}

	int count_stable(int * decay_channel_particles)
	{
		int count = 0;
		for (int idcp = 0; idcp < Maxdecaypart; idcp++)
		{
			if (decay_channel_particles[idcp] == 0) continue;
			if (is_stable(decay_channel_particles[idcp])) count++;
		}
		return(count);
	}

	bool is_stable(int monval)
	{
		int local_idx = lookup_particle_id_from_monval(monval);
		return( all_particles[local_idx].decays_Npart[0] == 1 );
	}

	int set_stable_particle_monval()
	{
		int local_Nstable_particle;
		int Idummy;
		char cdummy[256];
		ifstream particletable("./decay/EOS/EOS_particletable.dat");
		particletable >> local_Nstable_particle;
		stable_particle_monval = new double [local_Nstable_particle];
		for(int i=0; i<local_Nstable_particle; i++)
		{
			particletable >> Idummy >> stable_particle_monval[i];
			particletable.getline(cdummy, 256);
		}
		particletable.close();

		return(local_Nstable_particle);
	}

	void print_particle_stability()
	{
		int Nparticle = (int)all_particles.size();
		for (int ipart = 0; ipart < Nparticle; ipart++)
			cout << all_particles[ipart].name << "   " << all_particles[ipart].monval << "   "
					<< all_particles[ipart].mass << "   " << all_particles[ipart].stable << endl;
		return;
	}

	int get_number_of_decay_channels(vector<int> & chosen_resonances)
	{
		int count = 0, total_number_of_decays = 0;
		for (int icr = 0; icr < (int)chosen_resonances.size(); icr++)
		{
			total_number_of_decays = all_particles[chosen_resonances[icr]].decays;
			count += total_number_of_decays;
		}
		return count;
	}



	void get_important_resonances(
			int chosen_target_particle_idx,
			vector<int> & chosen_resonance_indices,
			double threshold, double &running_total_percentage )
	{
		//**********************************************************************************
		//SELECT RESONANCES TO INCLUDE IN SOURCE VARIANCES CALCULATIONS
		//sort resonances by importance and loop over all resonances needed to achieve a certain minimum percentage of total decay pions
		//double threshold = 0.6;	//include only enough of the most important resonances to account for fixed fraction of total resonance-decay pion(+)s
						//threshold = 1.0 means include all resonance-decay pion(+)s,
						//threshold = 0.0 means include none of them
		vector<double> percent_contributions;
		int Nparticle = (int)all_particles.size();
		for (int i = 0; i < Nparticle; i++)
			percent_contributions.push_back(all_particles[i].percent_contribution);
		vector<size_t> sorted_resonance_indices = ordered(percent_contributions);
		reverse(sorted_resonance_indices.begin(), sorted_resonance_indices.end());
		//vector<int> chosen_resonance_indices;
		//double running_total_percentage = 0.0;
		int count = 0;
		if (threshold < 1e-12)
		{
			cout << "No resonances included." << endl;
		}
		else if (fabs(1. - threshold) < 1e-12)
		{
			count = Nparticle;	//			if (sorted_resonance_indices[ii - 1] == chosen_target_particle_idx)
			for (int ii = 1; ii <= count; ii++)
				chosen_resonance_indices.push_back(sorted_resonance_indices[ii - 1]);
			running_total_percentage = 1.0;
			cout << "All resonances included." << endl;
		}
		else
		{
			while (running_total_percentage <= threshold)
			{
				running_total_percentage += 0.01 * all_particles[sorted_resonance_indices[count]].percent_contribution;
				chosen_resonance_indices.push_back(sorted_resonance_indices[count]);
				count++;
			}
			if (chosen_resonance_indices.size() == 0)
			{
				cout << "No resonances included!  Choose a higher threshold!" << endl;
				exit(1);
			}
			else
			{
				cout << "Including the following " << count << " resonances (accounting for " << 100.*running_total_percentage
					<< "%, threshold = " << 100.*threshold << "%): " << endl;
				for (int ii = 1; ii <= count; ii++)
					cout << "\t" << ii << ": " << all_particles[sorted_resonance_indices[ii - 1]].name << endl;
			}
		}
		//END OF CODE SELECTING INCLUDED RESONANCES
		//**********************************************************************************
	
		return;
	}

	void get_all_descendants(vector<int> & chosen_resonance_indices)
	{
		// This function ensures that no particles are missing from the chosen_resonances vector,
		// even if a truncated (threshold < 1.0) version is used
		int amount_of_output = 0;
		bool no_missing_descendants = false;
		int original_size = (int)chosen_resonance_indices.size();
		while (!no_missing_descendants)
		{
			int old_size = (int)chosen_resonance_indices.size();
			int new_size = old_size;
			for (int icr = 0; icr < old_size; icr++)
			{
				particle_info resonance = all_particles[chosen_resonance_indices[icr]];
				if (amount_of_output > 1) cout << "Currently looking at resonance " << resonance.name << endl;
				if (resonance.stable == 1 && resonance.decays_Npart[0] == 1)
					continue;
				int number_of_decays = resonance.decays;
				for (int k = 0; k < number_of_decays; k++)
				{
					if (amount_of_output > 1) cout << resonance.name << ": currently looking at decay #" << k << endl;
					int nb = abs(resonance.decays_Npart[k]);
					for (int l = 0; l < nb; l++)
					{
						int pid = lookup_particle_id_from_monval(resonance.decays_part[k][l]);
						if (amount_of_output > 1) cout << resonance.name << ", decay #" << k
									<< ": currently looking at daughter " << all_particles[pid].name << endl;
						//if this decay particle is not in the list and has large effective br, include it
						bool decay_particle_is_not_in_list
								= ( find ( chosen_resonance_indices.begin(),
											chosen_resonance_indices.end(), pid )
									==
									chosen_resonance_indices.end() );
						bool br_tot_is_not_small = ( all_particles[pid].effective_branchratio >= 1.e-12 );
						if (decay_particle_is_not_in_list)
							if (amount_of_output > 1) cout << resonance.name << ", decay #" << k
								<< ", daughter " << all_particles[pid].name << ": daughter is not in list!" << endl;
						if (br_tot_is_not_small)
							if (amount_of_output > 1) cout << resonance.name << ", decay #" << k
								<< ", daughter " << all_particles[pid].name << ": effective br. is not small!" << endl;
						if (decay_particle_is_not_in_list && br_tot_is_not_small)
						{
							chosen_resonance_indices.push_back(pid);
							if (amount_of_output > 1) cout << resonance.name << ", decay #" << k
								<< ", daughter " << all_particles[pid].name << ": adding daughter to list!" << endl;
							new_size++;
						}
					}
				}
			}
			no_missing_descendants = ( old_size == new_size );
			if (new_size == original_size)
				if (amount_of_output > 0) cout << "get_all_descendants(): No new particles added!" << endl;
			else
				if (amount_of_output > 0) cout << "get_all_descendants(): " << new_size - original_size << " new particles added!" << endl;
		}
		return;
	}

	void sort_by_mass(vector<int> & chosen_resonance_indices)
	{	//with the heaviest first
		//cout << "sort_by_mass(): Sorting by mass..." << endl;
		int number_of_chosen_particles = (int)chosen_resonance_indices.size();
		for (int m = 0; m < number_of_chosen_particles; m++)
		for (int n = 0; n < number_of_chosen_particles - m - 1; n++)
		if (all_particles[chosen_resonance_indices[n]].mass < all_particles[chosen_resonance_indices[n + 1]].mass)
		{
			// swap them
			int particle_idx = chosen_resonance_indices[n + 1];
			chosen_resonance_indices[n + 1] = chosen_resonance_indices[n];
			chosen_resonance_indices[n] = particle_idx;
		}
		//cout << "sort_by_mass(): ...finished!" << endl;
		return;
	}

}

//End of file
