#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <complex>

#include <gsl/gsl_integration.h>

using namespace std;

#include "lib.h"

double integrate_1D(double (*f)(double), double * xpts, double * xwts, int nx)
{
    double sum = 0.0;
    for (int ix = 0; ix < nx; ix++)
        sum += xwts[ix] * (*f)(xpts[ix]);

    return (sum);
}

double integrate_1D(double (*f)(double, void *), double * xpts, double * xwts, int nx, void * p)
{
    double sum = 0.0;
    for (int ix = 0; ix < nx; ix++)
        sum += xwts[ix] * (*f)(xpts[ix], p);

    return (sum);
}

complex<double> integrate_1D_FT(double (*f)(double), double * xpts, double * xwts, int nx, double k)
{
    complex<double> sum = (0,0);
	for (int ix = 0; ix < nx; ix++)
		sum += xwts[ix] * exp(- i * k * xpts[ix]) * (*f)(xpts[ix]);

    return (sum);
}

complex<double> integrate_1D_FT(double (*f)(double, void *), double * xpts, double * xwts, int nx, double k, void * p)
{
    complex<double> sum = (0,0);
	for (int ix = 0; ix < nx; ix++)
	{
		sum += xwts[ix] * exp(- i * k * xpts[ix]) * (*f)(xpts[ix], p);
		//cout << "FT: " << xpts[ix] << "   " << xwts[ix] << "   " << k << "   "
		//		<< exp(- i * k * xpts[ix]).real() << "   " << exp(- i * k * xpts[ix]).imag() << "   "
		//		<< (*f)(xpts[ix], p) << endl;
	}
	//if (1) exit(0);

    return (sum);
}

double integrate_2D(double (*f)(double, double), double * xpts, double * ypts, double * xwts, double * ywts, int nx, int ny)
{
	double sum = 0.0;
	for (int ix = 0; ix < nx; ix++)
	for (int iy = 0; iy < ny; iy++)
		sum += xwts[ix] * ywts[iy] * (*f)(xpts[ix], ypts[iy]);

	return (sum);
}

/*double integrate_1D(double (*f)(double), double Lx, int nx)
{
    double sum = 0.0;
    for (int ix = 1; ix < nx; ix++)
	{
		double xp = (double)ix * M_PI / (double)nx;
		double s1 = sin(xp);
        sum += (Lx * M_PI) * (*f)(Lx * cot(xp)) / (nx*s1*s1);
	}

    return (sum);
}


double integrate_2D(double (*f)(double, double), double Lx, double Ly, int nx, int ny)
{
    double sum = 0.0;
    for (int ix = 1; ix < nx; ix++)
    for (int iy = 1; iy < ny; iy++)
	{
		double xp = (double)ix * M_PI / (double)nx;
		double yp = (double)iy * M_PI / (double)ny;
		double s1 = sin(xp);
		double s2 = sin(yp);
        sum += (Lx * Ly * M_PI * M_PI) * (*f)(Lx * cot(xp), Ly * cot(yp)) / (nx*ny*s1*s1*s2*s2);
	}

    return (sum);
}*/

// some interpolation routines here
double interpolate1D(double * x, double * y, double x0, long size, int kind, bool uniform_spacing, bool returnflag /*= false*/, double default_return_value /* = 0*/)
{
// kind == 0: linear interpolation
// kind == 1: cubic interpolation
	switch (kind)
	{
		case 0:
		{
			if (uniform_spacing)
				return interpLinearDirect(x, y, x0, size, returnflag, default_return_value);
			else
				return interpLinearNondirect(x, y, x0, size, returnflag, default_return_value);
			break;
		}
		case 1:
		{
			if (uniform_spacing)
				return interpCubicDirect(x, y, x0, size, returnflag, default_return_value);
			else
				return interpCubicNonDirect(x, y, x0, size, returnflag, default_return_value);
			break;
		}
		default:
		{
			cerr << "Error (interpolate1D): interpolation kind not supported!" << endl;
			exit(1);
			break;
		}
	}
	return default_return_value;
}

//**********************************************************************
double interpLinearDirect(double * x, double * y, double x0, long size, bool returnflag /*= false*/, double default_return_value /* = 0*/)
// Returns the interpreted value of y=y(x) at x=x0 using linear interpolation method.
// -- x,y: the independent and dependent tables; x is assumed to be equal spaced and increasing
// -- x0: where the interpolation should be performed
{
	if (size==1) {cout<<"interpLinearDirect warning: table size = 1"<<endl; return y[0];}
	double dx = x[1]-x[0]; // increment in x

	// if close to left end:
	if (abs(x0-x[0])<dx*1e-30) return y[0];

	// find x's integer index
	long idx = floor((x0-x[0])/dx);

	if (idx<0 || idx>=size-1)
	{
		if (!returnflag)	//i.e., if returnflag is false, exit
		{
			cout    << "interpLinearDirect: x0 out of bounds." << endl
				<< "x ranges from " << x[0] << " to " << x[size-1] << ", "
				<< "x0=" << x0 << ", " << "dx=" << dx << ", " << "idx=" << idx << endl;
			exit(1);
		}
		else
		{
			idx = (idx<0) ? 0 : size-2;	//uses extrapolation
		}
		//else return default_return_value;
	}

  return y[idx] + (y[idx+1]-y[idx])/dx*(x0-x[idx]);
}

//**********************************************************************
double interpLinearNondirect(double * x, double * y, double x0, long size, bool returnflag /*= false*/, double default_return_value /* = 0*/)
// Returns the interpreted value of y=y(x) at x=x0 using linear interpolation method.
// -- x,y: the independent and dependent tables; x is assumed to be increasing but not equal spaced
// -- x0: where the interpolation should be performed
{
	if (size==1) {cout<<"interpLinearNondirect warning: table size = 1"<<endl; return y[0];}
	double dx = x[1]-x[0]; // increment in x

	// if close to left end:
	if (abs(x0-x[0])<dx*1e-30) return y[0];

	// find x's integer index
	long idx = binarySearch(x, size, x0, true);
	if (idx<0 || idx>=size-1)
	{
		if (!returnflag)	//i.e., if returnflag is false, exit
		{
			cout    << "interpLinearNondirect: x0 out of bounds." << endl
				<< "x ranges from " << x[0] << " to " << x[size-1] << ", "
				<< "x0=" << x0 << ", " << "dx=" << dx << ", " << "idx=" << idx << endl;
			exit(1);
		}
		else
		{
			idx = (idx<0) ? 0 : size-2;	//uses extrapolation
			//cout //<< "interpLinearNondirect: x0 out of bounds." << endl
			//		<< "interpLinearNondirect(" << (int)default_return_value << "): x ranges from " << x[0] << " to " << x[size-1] << ", "
			//		<< "x0=" << x0 << ", " << "dx=" << dx << ", " << "idx=" << idx << endl
			//		<< "interpLinearNondirect: data = " << x[idx] << "   " << x[idx+1] << "   " << y[idx] << "   " << y[idx+1] << endl
			//		<< "interpLinearNondirect: result = " << y[idx] + (y[idx+1]-y[idx])/(x[idx+1]-x[idx])*(x0-x[idx]) << endl;
			//idx = (idx<0) ? 0 : size-2;	//uses extrapolation
		}
		//else return default_return_value;
	}

	return y[idx] + (y[idx+1]-y[idx])/(x[idx+1]-x[idx])*(x0-x[idx]);
}

//**********************************************************************
double interpCubicDirect(double * x, double * y, double x0, long size, bool returnflag /*= false*/, double default_return_value /* = 0*/)
// Returns the interpreted value of y=y(x) at x=x0 using cubic polynomial interpolation method.
// -- x,y: the independent and dependent tables; x is assumed to be equal spaced and increasing
// -- x0: where the interpolation should be performed
{
  if (size==1) {cout<<"interpCubicDirect warning: table size = 1"; return y[0];}
  double dx = x[1]-x[0]; // increment in x

  // if close to left end:
  if (abs(x0-x[0])<dx*1e-30) return y[0];

  // find x's integer index
  long idx = floor((x0-x[0])/dx);

	if (idx < 0 || idx >= size-1)
	{
		if (!returnflag)	//i.e., if returnflag is false, exit
		{
			cerr << "interpCubicDirect(): index out of range!  Aborting!" << endl
				<< "interpCubicDirect(): size = " << size << ", x0 = " << x0 << ", " << "dx=" << dx << ", " << "idx=" << idx << endl;
			exit(1);
		}
		else
		{
			idx = (idx<0) ? 0 : size-2;	//uses linear extrapolation
			return y[idx] + (y[idx+1]-y[idx])/(x[idx+1]-x[idx])*(x0-x[idx]);
		}
		//else return default_return_value;
	}

  if (idx==0)
  {
    // use quadratic interpolation at left end
    double A0 = y[0], A1 = y[1], A2 = y[2], deltaX = x0 - x[0]; // deltaX is the increment of x0 compared to the closest lattice point
    return (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX - (3.0*A0-4.0*A1+A2)/(2.0*dx)*deltaX + A0;
  }
  else if (idx==size-2)
  {
    // use quadratic interpolation at right end
    double A0 = y[size-3], A1 = y[size-2], A2 = y[size-1], deltaX = x0 - (x[0] + (idx-1)*dx);
    return (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX - (3.0*A0-4.0*A1+A2)/(2.0*dx)*deltaX + A0;
  }
  else
  {
    // use cubic interpolation
    double A0 = y[idx-1], A1 = y[idx], A2 = y[idx+1], A3 = y[idx+2], deltaX = x0 - (x[0] + idx*dx);
    //cout << A0 << "  " << A1 << "  " << A2 << "  " << A3 << endl;
    return (-A0+3.0*A1-3.0*A2+A3)/(6.0*dx*dx*dx)*deltaX*deltaX*deltaX
            + (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX
            - (2.0*A0+3.0*A1-6.0*A2+A3)/(6.0*dx)*deltaX
            + A1;
  }

}

//**********************************************************************
double interpCubicNonDirect(double * x, double * y, double xi, long size, bool returnflag /*= false*/, double default_return_value /* = 0*/)
// Returns the interpreted value of y=y(x) at x=x0 using cubic polynomial interpolation method.
// -- x,y: the independent and dependent double x0ables; x is *NOT* assumed to be equal spaced but it has to be increasing
// -- xi: where the interpolation should be performed
{
  if (size==1) {cout<<"interpCubicNondirect warning: table size = 1"<<endl; return y[0];}

  // if close to left end:
  if (abs(xi-x[0])<(x[1]-x[0])*1e-30) return y[0];

  // find x's integer index
  long idx = binarySearch(x, size, xi, true);

	if (idx < 0 || idx >= size-1)
	{
		if (!returnflag)	//i.e., if returnflag is false, exit
		{
			cerr << "interpCubicNonDirect(): index out of range!  Aborting!" << endl
				<< "interpCubicNonDirect(): size = " << size << ", x0 = " << xi << ", " << "idx=" << idx << endl;
			exit(1);
		}
		else
		{
			idx = (idx<0) ? 0 : size-2;	//uses linear extrapolation
		}
		//else return default_return_value;
	}

  if (idx==0)
  {
    // use linear interpolation at the left end
    return y[0] + (y[1]-y[0])/(x[1]-x[0])*(xi-x[0]);
  }
  else if (idx==size-2)
  {
    // use linear interpolation at the right end
    return y[size-2] + (y[size-1]-y[size-2] )/(x[size-1]-x[size-2] )*(xi-x[size-2]);
  }
  else
  {
    // use cubic interpolation
    long double y0 = y[idx-1], y1 = y[idx], y2 = y[idx+1], y3 = y[idx+2];
    long double y01=y0-y1, y02=y0-y2, y03=y0-y3, y12=y1-y2, y13=y1-y3, y23=y2-y3;
    long double x0 = x[idx-1], x1 = x[idx], x2 = x[idx+1], x3 = x[idx+2];
    long double x01=x0-x1, x02=x0-x2, x03=x0-x3, x12=x1-x2, x13=x1-x3, x23=x2-x3;
    long double x0s=x0*x0, x1s=x1*x1, x2s=x2*x2, x3s=x3*x3;
    long double denominator = x01*x02*x12*x03*x13*x23;
    long double C0, C1, C2, C3;
    C0 = (x0*x02*x2*x03*x23*x3*y1
          + x1*x1s*(x0*x03*x3*y2+x2s*(-x3*y0+x0*y3)+x2*(x3s*y0-x0s*y3))
          + x1*(x0s*x03*x3s*y2+x2*x2s*(-x3s*y0+x0s*y3)+x2s*(x3*x3s*y0-x0*x0s*y3))
          + x1s*(x0*x3*(-x0s+x3s)*y2+x2*x2s*(x3*y0-x0*y3)+x2*(-x3*x3s*y0+x0*x0s*y3))
          )/denominator;
    C1 = (x0s*x03*x3s*y12
          + x2*x2s*(x3s*y01+x0s*y13)
          + x1s*(x3*x3s*y02+x0*x0s*y23-x2*x2s*y03)
          + x2s*(-x3*x3s*y01-x0*x0s*y13)
          + x1*x1s*(-x3s*y02+x2s*y03-x0s*y23)
          )/denominator;
    C2 = (-x0*x3*(x0s-x3s)*y12
          + x2*(x3*x3s*y01+x0*x0s*y13)
          + x1*x1s*(x3*y02+x0*y23-x2*y03)
          + x2*x2s*(-x3*y01-x0*y13)
          + x1*(-x3*x3s*y02+x2*x2s*y03-x0*x0s*y23)
          )/denominator;
    C3 = (x0*x03*x3*y12
          + x2s*(x3*y01+x0*y13)
          + x1*(x3s*y02+x0s*y23-x2s*y03)
          + x2*(-x3s*y01-x0s*y13)
          + x1s*(-x3*y02+x2*y03-x0*y23)
          )/denominator;
    return C0 + C1*xi + C2*xi*xi + C3*xi*xi*xi;
  }
}

//**********************************************************************
long binarySearch(double * A, int length, double value, bool skip_out_of_range /*== true*/, bool verbose /*== false*/)
// Return the index of the largest number less than value in the list A
// using binary search. Index starts with 0.
// If skip_out_of_range is set to true, then it will return -1 for those
// samples that are out of the table range (default is true).
{
   int idx_i, idx_f, idx;
   idx_i = 0;
   idx_f = length-1;

   if(value > A[idx_f])
   {
      if (verbose)
		cerr << "binarySearch: desired value is too large, exceeding the end of the table: value = "
				<< value << " and A[idx_f] = " << A[idx_f] << endl;
      if (skip_out_of_range) return length;
      exit(1);
   }
   if(value < A[idx_i])
   {
      if (verbose)
		cerr << "binarySearch: desired value is too small, exceeding the beginning of table: value = "
				<< value << " and A[idx_i] = " << A[idx_i] << endl;
      if (skip_out_of_range) return -1;
      exit(1);
   }
   idx = (int) floor((idx_f+idx_i)/2.);
   while((idx_f-idx_i) > 1)
   {
     if(A[idx] < value)
        idx_i = idx;
     else
        idx_f = idx;
     idx = (int) floor((idx_f+idx_i)/2.);
   }
   return(idx_i);
}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// some special functions below this line
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/*complex<double> Hypergeometric0F1(complex<double> a_in, complex<double> z_in)
{
	complex<double> a = a_in, z = z_in;

	complex<double> sum = 0.0, oldsum = 0.0;
	complex<double> term = 1.0;
	long int i = 0;
	//for (int i = 0; i <= max_terms; ++i)
	do
	{
		oldsum = sum;
		sum += term;
		term *= z/(a*double(i+1));
		a += 1.0;
		++i;
	} while (2.0*abs(sum-oldsum)/abs(sum+oldsum) >= tolerance);
	//cout << setw(15) << setprecision(12) << "Hypergeometric0F1: " << a_in.real() << "   " << a_in.imag()
	//		<< "   " << z_in.real() << "   " << z_in.imag() << "   " << sum.real() << "   " << sum.imag() << endl;
	return (sum);
}*/

complex<double> Hypergeometric1F1(complex<double> a_in, complex<double> b_in, complex<double> z_in)
{
	complex<double> a = a_in, b = b_in, z = z_in;

	complex<double> sum = 0.0, oldsum = 0.0;
	complex<double> term = 1.0;
	long int i = 0;
	//for (int i = 0; i <= max_terms; ++i)
	do
	{
		oldsum = sum;
		sum += term;
		term *= a*z/(b*double(i+1));
		a += 1.0;
		b += 1.0;
		//cout << i << "   " << oldsum << "   " << sum << endl;
		++i;
	} while (2.0*abs(sum-oldsum)/abs(sum+oldsum) >= tolerance);
	//cout << "Convergence reached at counter = " << i << "!" << endl;
	//cout << setw(15) << setprecision(12) << "Hypergeometric1F1: " << a_in.real() << "   " << a_in.imag() << "   " << b_in.real() << "   " << b_in.imag()
	//		<< "   " << z_in.real() << "   " << z_in.imag() << "   " << sum.real() << "   " << sum.imag() << endl;
	return (sum);
}

complex<double> Hypergeometric1F2(complex<double> a_in, complex<double> b0_in, complex<double> b1_in, complex<double> z_in)
{
	complex<double> a = a_in, b0 = b0_in, b1 = b1_in, z = z_in;

	complex<double> sum = 0.0, oldsum = 0.0;
	complex<double> term = 1.0;
	long int i = 0;
	do
	{
		oldsum = sum;
		sum += term;
		term *= a*z/(b0*b1*double(i+1));
		a += 1.0;
		b0 += 1.0;
		b1 += 1.0;
		++i;
	} while (2.0*abs(sum-oldsum)/abs(sum+oldsum) >= tolerance);
	return (sum);
}

/*
complex<double> SFL_Hypergeometric1F1(complex<double> a_in, complex<double> b_in, complex<double> z_in)
{
	struct SFL::complex a, b, z, ans;
	int top = -1;
	a.x = a_in.real();
	a.y = a_in.imag();
	b.x = b_in.real();
	b.y = b_in.imag();
	z.x = z_in.real();
	z.y = z_in.imag();

	SFL::c1f1(&a, &b, &z, top, &ans);

	complex<double> ans_out(ans.x, ans.y);

	//cout << setw(15) << setprecision(12) << "SFL_Hypergeometric1F1: " << a_in.real() << "   " << a_in.imag() << "   " << b_in.real() << "   " << b_in.imag()
	//		<< "   " << z_in.real() << "   " << z_in.imag() << "   " << ans_out.real() << "   " << ans_out.imag() << endl;

	return (ans_out);
}
*/

//End of file
