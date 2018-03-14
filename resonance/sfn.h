#ifndef SFN_H
#define SFN_H

namespace sfn
{
	//function prototypes
	void
		operation_mode_0(
		bool skip_thermal );
	vector<complex<double> >
		operation_mode_1(
		vector<double> Del_pY_pts_in,
		vector<double> k_pts_in,
		bool skip_thermal );
}

#endif
