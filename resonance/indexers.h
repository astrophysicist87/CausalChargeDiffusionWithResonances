#ifndef INDEXERS_H
#define INDEXERS_H

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

inline int res_vector_indexer2(const int ik, const int iRes, const int ipT, const int ipphi, const int ipY)
{
	return (
		( ( ( ik * n_resonances + iRes ) * n_pT_pts + ipT ) * n_pphi_pts + ipphi ) * n_pY_pts + ipY
	);
}

inline int res_FIX_K_vector_indexer(const int iRes, const int ipT, const int ipphi, const int ipY)
{
	return (
		( ( iRes * n_pT_pts + ipT ) * n_pphi_pts + ipphi ) * tmp_n_pY_pts + ipY
	);
}

inline int res_FIX_K_vector_indexer2(const int iRes, const int ipT, const int ipphi, const int ipY)
{
	return (
		( ( iRes * n_pT_pts + ipT ) * n_pphi_pts + ipphi ) * n_pY_pts + ipY
	);
}

inline int scriptFn_vector_indexer(const int ik, const int ipY)
{
	return (
		ik * n_pY_pts + ipY
	);
}

inline int mom_indexer(const int ipT, const int ipphi, const int ipY)
{
	return (
		( ipT * n_pphi_pts + ipphi ) * tmp_n_pY_pts + ipY
	);
}

// End of file

#endif
