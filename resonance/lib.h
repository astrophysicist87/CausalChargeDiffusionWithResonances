#ifndef LIB_H
#define LIB_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <complex>

#include <gsl/gsl_integration.h>

using namespace std;

namespace reslib
{

	const complex<double> i(0, 1);

	const double tolerance = 1.e-15;

	inline complex<double> cot(complex<double> x){return (cos(x)/sin(x));}
	inline complex<double> csc(complex<double> x){return (1.0/sin(x));}
	inline complex<double> sec(complex<double> x){return (1.0/cos(x));}

	void linspace(vector<double> & x, double a, double b);

	double integrate_1D(double (*f)(double), double * xpts, double * xwts, int nx);
	double integrate_1D(double (*f)(double, void *), double * xpts, double * xwts, int nx, void * p);
	complex<double> integrate_1D_FT(double (*f)(double), double * xpts, double * xwts, int nx, double k);
	complex<double> integrate_1D_FT(double (*f)(double, void *), double * xpts, double * xwts, int nx, double k, void * p);
	double integrate_2D(double (*f)(double, double), double * xpts, double * ypts, double * xwts, double * ywts, int nx, int ny);

	//main interpolation routine
	double interpolate1D(double * x, double * y, double x0, long size, int kind, bool uniform_spacing,
				bool returnflag = false, double default_return_value = 0.0);

	//subsidiary interpolation routines
	double interpLinearDirect(double * x, double * y, double x0, long size,
					bool returnflag = false, double default_return_value = 0.0);
	double interpLinearNondirect(double * x, double * y, double x0, long size,
					bool returnflag = false, double default_return_value = 0.0);
	double interpCubicDirect(double * x, double * y, double x0, long size,
					bool returnflag = false, double default_return_value = 0.0);
	double interpCubicNonDirect(double * x, double * y, double x0, long size,
					bool returnflag = false, double default_return_value = 0.0);
	long binarySearch(double * A, int length, double value,
					bool skip_out_of_range = true, bool verbose = false);

	//////////////////////////////////
	// special functions
	//////////////////////////////////

	complex<double> Hypergeometric1F1(complex<double> a, complex<double> b, complex<double> z);
	complex<double> Hypergeometric1F2(complex<double> a_in, complex<double> b0_in, complex<double> b1_in, complex<double> z_in);

	//stuff I have to define in a header file
	template <typename T> int sgn(T val)
	{
		return (T(0) < val) - (val < T(0));
	}
}

#endif
