#ifndef ASYMPTOTICS_H
#define ASYMPTOTICS_H

#include <complex>
#include <iostream>
#include <cmath>
#include <vector>

#include <gsl/gsl_sf_airy.h>

using namespace std;

namespace asymptotics
{
	//small number to avoid divisions by zero, etc.
	const double eps = 1.0e-10;

	//some useful numerical constants used below
	const complex<double> one = 1.0;
	const complex<double> i(0, 1);
	const complex<double> m1_to_1_3 = 0.5 + 0.86602540378443864676 * i;
	const complex<double> three_by_two_to_two_by_three = 1.3103706971044483036;

	inline void get_Bi_nu_and_Bi_nu_prime(complex<double> nu, complex<double> z, complex<double> & Bi_nu, complex<double> & Bi_nu_prime)
	{
		double z_re = z.real();

		double Ai = gsl_sf_airy_Ai (z_re, GSL_PREC_DOUBLE);
		double Ai_prime = gsl_sf_airy_Ai_deriv (z_re, GSL_PREC_DOUBLE);
		double Bi = gsl_sf_airy_Bi (z_re, GSL_PREC_DOUBLE);
		double Bi_prime = gsl_sf_airy_Bi_deriv (z_re, GSL_PREC_DOUBLE);

		if (0)
		{
			if (z_re < 0.0)
			{
				double sz = sqrt(-z_re);
				double lz = 2.0*sz*sz*sz/3.0;
				Ai = cos(lz - 0.25*M_PI) / (sqrt(M_PI)*pow(z_re,0.25));
				Ai_prime = pow(z_re,0.25)*sin(lz - 0.25*M_PI) / sqrt(M_PI);
				Bi = -sin(lz - 0.25*M_PI) / (sqrt(M_PI)*pow(z_re,0.25));
				Bi_prime = pow(z_re,0.25)*cos(lz - 0.25*M_PI) / sqrt(M_PI);
			}
			else
			{
				double sz = sqrt(z_re);
				double lz = 2.0*sz*sz*sz/3.0;
				Ai = exp(-lz) / (2.0*sqrt(M_PI)*pow(z_re,0.25));
				Ai_prime = -pow(z_re,0.25)*exp(-lz) / (2.0*sqrt(M_PI));
				Bi = exp(lz) / (sqrt(M_PI)*pow(z_re,0.25));
				Bi_prime = pow(z_re,0.25)*exp(lz) / (sqrt(M_PI));
			}
		}

		Bi_nu = Bi - i * tanh(M_PI*nu) * Ai;
		Bi_nu_prime = Bi_prime - i * tanh(M_PI*nu) * Ai_prime;

		return;
	}

	inline double zeta(double x)
	{
		complex<double> phase = complex<double>( x <= 1.0 ) - m1_to_1_3 * complex<double>( x > 1.0 );
		phase *= three_by_two_to_two_by_three;
		complex<double> root = sqrt(one - x*x);

		return ( ( phase * pow(
								log( ( one + root ) / x ) - root,
								2.0/3.0
							) ).real() );
	}

	inline double zeta_prime(double x)
	{
		complex<double> phase = complex<double>( x <= 1.0 ) - m1_to_1_3 * complex<double>( x > 1.0 );
		phase *= pow(2.0/3.0,1.0/3.0);
		complex<double> root = sqrt(one - x*x);

		return ( ( -phase * root * pow(
								log( ( one + root ) / x ) - root,
								-1.0/3.0
							) / x ).real() );
	}

	//not necessary
	inline double A_0(double x, double local_zeta)
	{
		return (1.0);
	}

	inline double B_0(double x, double local_zeta)
	{
		return (
				-5.0 / (48.0 * local_zeta * local_zeta)
				+ sqrt(4.0 * local_zeta / ( 1.0 - x * x ))
					* ( 5.0 / ( 1.0 - x * x ) - 3.0 )
					/ (48.0 * local_zeta)
				);
	}

	inline double B_0_prime(double x, double lz, double lzp)
	{
		double n1 = -6.0*x*(4.0+x*x)*lz*lz*lz;
		double n2 = (2.0 + x*x - 3.0*x*x*x*x)*lz*lz*lzp;
		double y = x*x-1.0;
		double n3 = 10.0*y*y*y*lzp*sqrt(-lz/y);
		double d = 48.0*y*y*y*lz*lz*lz*sqrt(-lz/y);
		return (
				( n1 + n2 + n3 ) / d
				);
	}


	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//CURRENTLY ASSUMING NU IS PURE IMAGINARY AND Z IS REAL!!!
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	inline complex<double> I(complex<double> nu_in, double z)
	{
		if ( abs( nu_in.real() ) > eps )
		{
			//return (0.0);
			cerr << "Can't handle real/complex nu=" << nu_in
					<< " in this asymptotic function!!!  Exiting..." << endl;
			exit(1);
			//double nu = nu_in.real();
			//return ( exp(z)/sqrt(2.0*M_PI*z)*( 1.0 - 0.125 * (4.0*nu*nu-1.0) / z ) );
		}

		double nu = nu_in.imag();
		if (nu_in.imag() <= 0.0) nu *= -1.0;
		double r = z/nu;
		double nu_to_1_3 = pow(nu, 1.0/3.0);
		double nu_to_2_3 = nu_to_1_3 * nu_to_1_3;
		double nu_to_4_3 = nu_to_2_3 * nu_to_2_3;
		double zeta_at_r = zeta(r);

		//assume z != 1.0
		double prefactor = 0.5 * exp(0.5*nu*M_PI) * pow( 4.0*zeta_at_r / (1.0-r*r), 0.25 ) / nu_to_1_3;
		complex<double> Bi_nu(0,0), Bi_nu_prime(0,0);
		double A0 = A_0(r, zeta_at_r);	//not necessary
		double B0 = B_0(r, zeta_at_r);

		//cout << "I(): " << setprecision(20) << nu << "   " << z << "   " << r
		//		<< "   " << zeta_at_r << "   " << -nu_to_2_3*zeta_at_r << endl;
		get_Bi_nu_and_Bi_nu_prime(nu, -nu_to_2_3*zeta_at_r, Bi_nu, Bi_nu_prime);
		//cout << "I(): " << setprecision(20) << nu << "   " << z << "   " << -nu_to_2_3*zeta_at_r 
		//		<< "   " << Bi_nu << "   " << Bi_nu_prime << endl;

		complex<double> pre_result = prefactor * ( Bi_nu + B0 * Bi_nu_prime / nu_to_4_3 );
		complex<double> final_result = 0.0;
		if (nu_in.imag() <= 0.0)
			final_result = conj( pre_result );
		else
			final_result = pre_result;

		return ( final_result );
	}

	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//CURRENTLY ASSUMING NU IS PURE IMAGINARY AND Z IS REAL!!!
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	inline complex<double> Iprime(complex<double> nu_in, double z)
	{
		if ( abs( nu_in.real() ) > eps )
		{
			//return (0.0);
			cerr << "Can't handle real/complex nu=" << nu_in
					<< " in this asymptotic function!!!  Exiting..." << endl;
			exit(1);
			//double nu = nu_in.real();
			//return ( exp(z)/sqrt(2.0*M_PI*z)*( 1.0 - 0.125 * (4.0*nu*nu+3.0) / z ) );
		}

		double nu = nu_in.imag();
		if (nu_in.imag() <= 0.0) nu *= -1.0;
		double r = z/nu;
		double nu_to_1_3 = pow(nu, 1.0/3.0);
		double nu_to_2_3 = nu_to_1_3 * nu_to_1_3;
		double nu_to_4_3 = nu_to_2_3 * nu_to_2_3;
		double nu_to_8_3 = nu_to_4_3 * nu_to_4_3;
		double zeta_at_r = zeta(r);
		double y = 1.0-r*r;

		//assume z != 1.0
		double prefactor = exp(0.5*nu*M_PI)
							/ ( 4.0*sqrt(2.0)*nu_to_8_3*y*y
								* pow( zeta_at_r / y, 0.75 ) );
		complex<double> Bi_nu(0,0), Bi_nu_prime(0,0);
		double A0 = A_0(r, zeta_at_r);	//not necessary
		double B0 = B_0(r, zeta_at_r);
		double zeta_prime_at_r = zeta_prime(r);
		double B0prime = B_0_prime(r, zeta_at_r, zeta_prime_at_r);

		//cout << "Iprime(): " << setprecision(20) << nu << "   " << z
		//		 << "   " << zeta_at_r << "   " << -nu_to_2_3*zeta_at_r << endl;
		get_Bi_nu_and_Bi_nu_prime(nu, -nu_to_2_3*zeta_at_r, Bi_nu, Bi_nu_prime);
		//cout << "Iprime(): " << setprecision(20) << nu << "   " << z << "   " << -nu_to_2_3*zeta_at_r << "   "
		//		<< Bi_nu << "   " << Bi_nu_prime << endl;

		double factor1 = B0*(y*zeta_prime_at_r + 2.0*r*zeta_at_r)
							+ 4.0*y*zeta_at_r*(B0prime - nu*nu*zeta_prime_at_r);
		double factor2 = 2.0*r*zeta_at_r
							+ y*zeta_prime_at_r*(4.0*B0*zeta_at_r*zeta_at_r + 1.0);


		complex<double> pre_result = prefactor * ( factor1 * Bi_nu_prime
													+ nu_to_4_3 * factor2 * Bi_nu );
		complex<double> final_result = 0.0;
		if (nu_in.imag() <= 0.0)
			final_result = conj( pre_result );
		else
			final_result = pre_result;

		return ( final_result );
	}
}

#endif
