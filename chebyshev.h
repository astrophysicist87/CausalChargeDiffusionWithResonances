#ifndef CHEBYSHEV_H
#define CHEBYSHEV_H

#include <vector>

#include <gsl/gsl_math.h>
#include <gsl/gsl_chebyshev.h>

#include "chebyshev_library.h"

using namespace std;

namespace cheb_int
{
	inline void reset(vector<double> * v, double val = 0.0)
	{
		fill(v->begin(), v->end(), val);
	}


	void chebyshev_interpolate(
			vector<double> * old_nodes, vector<double> * old_values,
			vector<double> * new_nodes, vector<double> * new_values )
	{
		int npts = (int)old_nodes->size();
		double old_min = old_nodes->at(0);
		double old_max = old_nodes->at(npts-1);
		int Npts = (int)new_nodes->size();
		double new_min = new_nodes->at(0);
		double new_max = new_nodes->at(npts-1);

		reset(new_values);

		vector<double> x_pts(npts);
		for (int k = 0; k < npts; ++k)
			x_pts.at(k) = - cos( M_PI*(2.*(k+1.) - 1.) / (2.*npts) );

		int n_coeffs = npts;
		double coeffs_array[npts];
		double nums[npts*npts];
		double dens[npts];

		for (int j = 0; j < npts; ++j)
		{
			dens[j] = 0.0;
			for (int k = 0; k < npts; ++k)
			{
				double Tjk = csf::Tfun(j, x_pts[k]);
				dens[j] += Tjk*Tjk;
				nums[j*npts+k] = Tjk;
			}
		}

		int icf = 0;
		double * cfs_for_eval = new double [npts];

		//////////////////////////////////
		//separate out 0th coefficient for additional factor of 2.0
		coeffs_array[0] = 0.0;
		for (int k = 0; k < npts; ++k)
			coeffs_array[0] += 2.0*old_values->at(k) * nums[0*npts+k];

		cfs_for_eval[icf++] = coeffs_array[0] / dens[0];
		for (int j = 1; j < npts; ++j)
		{
			coeffs_array[j] = 0.0;
			for (int k = 0; k < npts; ++k)
				coeffs_array[j] += old_values->at(k) * nums[j*npts+k];
			cfs_for_eval[icf++] = coeffs_array[j] / dens[j];
		}

		//////////////////////////////////
		gsl_cheb_series *cs;
		cs = gsl_cheb_alloc (npts - 1);
		cs->a = old_min;
		cs->b = old_max;
		cs->c = cfs_for_eval;
//cout << old_min << "   " << old_max << endl;
//for (int i = 0; i < npts; ++i)
//	cout << cfs_for_eval[i] << endl;

		for (int iN = 0; iN < Npts; ++iN)
		{
	//cout << "iN = " << iN << " of " << new_nodes->size() << " and " << new_values->size() << endl;
			double new_point = new_nodes->at(iN);
			new_values->at(iN) = gsl_cheb_eval (cs, new_point);
		}

		delete [] cfs_for_eval;

		return;
	}
}

#endif
