#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <complex>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>

using namespace std;

#include "decay/parameters.h"
#include "lib.h"
#include "decay/readindata.h"
#include "gauss/gauss_quadrature.h"
#include "decay/decay.h"
#include "defs.h"
#include "thermal.h"
#include "chebyshev.h"
#include "main.h"
#include "Stopwatch.h"

const int operation_mode = 1;
// 0 - run test of normal spectra (using exact toy model defined in thermal.h)
// 1 - try to get corresponding results for FT'd thermal "spectra" needed for calculation

/////////////////////////////////
int main(int argc, char *argv[])
{
	Stopwatch sw;
	
	sw.Start();
	switch (operation_mode)
	{
		case 0:
			operation_mode_0();
			break;
		case 1:
			operation_mode_1(true);
			break;
		default:
			cout << "Option not supported!" << endl;
			break;
	}
	sw.Stop();
	cout << "Finished everything in " << sw.printTime() << " seconds." << endl;

	return 0;
}

// End of file
