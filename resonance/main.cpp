#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <complex>

using namespace std;

#include "sfn.h"
#include "Stopwatch.h"

const int operation_mode = 1;
// 0 - run test of normal spectra (using exact toy model defined in thermal.h)
// 1 - try to get corresponding results for FT'd thermal "spectra" needed for calculation

/////////////////////////////////
int main(int argc, char *argv[])
{
	Stopwatch sw;
	
	vector<complex<double> > results;

	sw.Start();
	switch (operation_mode)
	{
		case 0:
			operation_mode_0(false);
			break;
		case 1:
			results = operation_mode_1(false);
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
