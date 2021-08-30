#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <numeric>
#include <tuple>
#include <fstream>
#include <ctime>
#include <chrono>
#include <stdio.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <omp.h>

// Defining special types
typedef std::complex<double> dcomp;
typedef std::vector<double> dvec;
typedef std::vector< std::vector<double> > ddvec;
typedef std::vector<int> ivec;
typedef std::vector< std::vector<int> > iivec;
typedef std::vector< std::complex<double> > dcvec;
typedef std::vector< std::vector<std::complex<double> > > ddcvec;

/* --------------------------------------------------------------------

		FUNCTIONS

-----------------------------------------------------------------------*/

inline const double pi()
{
	// Set teh value fo pi
	const double pi = std::atan(1) * 4;
	return pi;
}

/* --------------------------------------------------------------------
		FROM MICROSECONDS TO HH:MM:SS:MS
-----------------------------------------------------------------------*/
std::tuple<long long, long long, long long, long long> ms_to_time(long long ms)
{
	// Defining the variables
	long long s, m, h;

	// Calculating the time in hour:minutes:seconds:microseconds
	s = ms / 1000000;
	ms = ms % 1000000;
	m = s / 60;
	s = s % 60;
	h = m / 60;
	m = m % 60;

	return { ms,s,m,h };
}
/* --------------------------------------------------------------------
		PRINT in HH:MM:SS:MS FORMAT
-----------------------------------------------------------------------*/
void time_print(long long ms, long long s, long long m, long long h)
{

	// Print the time
	if (h < 10)
	{
		std::cout << "0" << h << ":";
	}
	else
	{
		std::cout << "" << h << ":";
	}
	if (m < 10)
	{
		std::cout << "0" << m << ":";
	}
	else
	{
		std::cout << "" << m << ":";
	}
	if (s < 10)
	{
		std::cout << "0" << s << ":" << ms;
	}
	else
	{
		std::cout << "" << s << ":" << ms;
	}

	return;
}
/* --------------------------------------------------------------------
		CHI FROM Z
-----------------------------------------------------------------------*/
inline dcomp chi_from_z(dcomp z, dcomp z1, dcomp z2, double L, double mu)
{
	// Defining the variables
	dcomp z0, Z, chi;

	// Calculating the chi from a z value
	z0 = 0.5 * (z1 + z2);
	Z = exp(dcomp(0, -1) * mu) * 2.0 * (z - z0) / L;
	chi = Z + sqrt(Z - 1.0) * sqrt(Z + 1.0);

	return chi;
}
/* --------------------------------------------------------------------
		Z FROM CHI
-----------------------------------------------------------------------*/
inline dcomp z_from_chi(dcomp chi, dcomp z1, dcomp z2, double L, double mu)
{
	// Defining the variables
	dcomp z, z0, Z;

	// Calculating the z from a chi value
	z0 = 0.5 * (z1 + z2);
	Z = 0.5 * (chi + 1.0 / chi);
	z = 0.5 * L * Z * exp(dcomp(0, 1) * mu) + z0;

	return z;
}
/* --------------------------------------------------------------------
		ANGEL CHANGE
-----------------------------------------------------------------------*/
inline double angel_change(dcomp z, dcomp z1, dcomp z2)
{
	// Defining the variables
	double u, v, eta;

	// Calculating the angels
	u = arg(z - z1);
	v = arg(z1 - z2);

	// Correct to angel of 0 to 2pi
	if (u < 0)
	{
		u += 2 * pi();
	}
	if (v < 0)
	{
		v += 2 * pi();
	}

	// Calculate angel between the vectors
	eta = abs(u - v);

	// Correct to angel of 0 to pi
	if (eta > pi())
	{
		eta -= 2 * pi();
	}

	eta = abs(eta);

	return eta;
}
/* --------------------------------------------------------------------
		FIND INTERSECTION POINT
-----------------------------------------------------------------------*/
inline std::tuple<dcomp, int> intersection_point(dcomp z1, dcomp z2, dcomp z3, dcomp z4)
{
	// Defining the variables
	dcomp zint, Z1, Z2;
	int int_check = 0;

	if (z1 != z3 || z2 != z4)
	{
		// Calculating the intersection point
		zint = ((conj(z2 - z1) * z1 - (z2 - z1) * conj(z1)) * (z4 - z3) - (conj(z4 - z3) * z3 - (z4 - z3) * conj(z3)) * (z2 - z1)) / ((z4 - z3) * conj(z2 - z1) - (z2 - z1) * conj(z4 - z3));

		// Check if intersection is on the lines
		Z1 = (zint - .5 * (z1 + z2)) / (.5 * (z1 - z2));
		Z2 = (zint - .5 * (z3 + z4)) / (.5 * (z3 - z4));

		if (real(Z1 * conj(Z1)) < 1.0)
		{
			if (real(Z2 * conj(Z2)) < 1.0)
			{
				int_check = 1;
			}
		}
		
	}

	return { zint, int_check };
}
/* --------------------------------------------------------------------
		GET DELTA
-----------------------------------------------------------------------*/
inline dcomp get_delta(ddcvec ab, int nab, int m, ivec pos_z1, ivec pos_z2)
{
	// Defintion of variables
	dcomp delta;

	delta = 0;
	for (int ii = 0; ii < nab; ii++)
	{
		// Adding the jump for the z1 elements
		if (pos_z1[ii] == 1)
		{
			for (int jj = 0; jj < m; jj++)
			{
				delta += ab[ii][jj]* pow(-1, jj);
			}
		}
		// Adding the jump for the z2 elements
		if (pos_z2[ii] == 1)
		{
			for (int jj = 0; jj < m; jj++)
			{
				delta -= ab[ii][jj];
			}
		}

	}

	return delta;
}
/* --------------------------------------------------------------------
		TAU UNIFORM STRESS
-----------------------------------------------------------------------*/
inline std::tuple<dcomp, dcomp>  tau_uni(double sigma_11inf)
{
	// Defining the variables
	dcomp tau_11, tau_12, phi, psi;

	// Get the phi and psi
	phi = -0.5 * sigma_11inf;
	psi = -0.5 * sigma_11inf;

	// calculating the tau
	tau_11 = -phi - psi;
	tau_12 = -phi - phi;

	return { tau_11, tau_12 };
}
/* --------------------------------------------------------------------
		TAU CRACK
-----------------------------------------------------------------------*/
inline std::tuple<dcomp, dcomp> tau_crack(dcomp z, dcomp z1, dcomp z2, double L, double mu, int m, dcvec a)
{
	// Defining the variables
	dcomp chi, chi_bar, Z, chi_pow;
	dcomp dphi, dphi_bar, ddphi, dpsi;
	dcomp tau_11, tau_12, S1, L_frac;

	// Getting the chi - and Z - coordinates
	chi = chi_from_z(z, z1, z2, L, mu);
	chi_bar = conj(chi);
	Z = exp(dcomp(0, -1) * mu) * 2.0 * (z - 0.5 * (z1 + z2)) / L;

	// Calculating the series
	dphi = 0;
	dphi_bar = 0;
	ddphi = 0;
	dpsi = 0;
	chi_pow = chi * chi - 1.0;
	for (int ii = 0; ii < m; ii++)
	{
		double n = ii + 1.0;
		dcomp beta_n = a[ii] * n;
		dcomp chipow = pow(chi, (1.0 - n)) / chi_pow;
		dphi += conj(beta_n) * chipow;
		dphi_bar += beta_n * conj(chipow);
		ddphi -= conj(beta_n) * (pow(chi, (2.0 - n)) / (chi_pow * chi_pow * chi_pow)) * ((n + 1.0) * chi * chi - n + 1.0);
		dpsi -= beta_n * chipow;
	}

	// Multiplying the constants
	L_frac = (4.0 / L) * exp(dcomp(0, -1) * mu);
	dphi *= L_frac;
	dphi_bar *= conj(L_frac);
	ddphi *= (16.0 / (L * L)) * exp(dcomp(0, -2) * mu);
	dpsi *= L_frac;

	// Calcualting tau
	tau_11 = -0.5 * L * (Z - conj(Z)) * ddphi - exp(dcomp(0, -1) * mu) * (dphi + dpsi);
	tau_12 = -exp(dcomp(0, 1) * mu) * dphi - exp(dcomp(0, -1) * mu) * dphi_bar;

	return { tau_11, tau_12 };
}
/* --------------------------------------------------------------------
		TAU TOTAL
-----------------------------------------------------------------------*/
inline std::tuple<dcomp, dcomp>  tau_total(dcomp z, double sigma_11inf, dcvec z1, dcvec z2, dvec L, dvec mu, int ma, int na, ddcvec a, int m_not_a)
{
	// Defining the variables
	dcomp tau_11, tau_12;

	// Add the unifrom stress field
	std::tie(tau_11, tau_12) = tau_uni(sigma_11inf);

	// Add the analytic element for a crack
	if (na > 0)
	{
		for (int ii = 0; ii < na; ii++)
		{
			if (ii != m_not_a)
			{
				dcomp tau_11c, tau_12c;
				std::tie(tau_11c, tau_12c) = tau_crack(z, z1[ii], z2[ii], L[ii], mu[ii], ma, a[ii]);
				tau_11 += tau_11c;
				tau_12 += tau_12c;
			}
		}
	}

	return { tau_11, tau_12 };
}
/* --------------------------------------------------------------------
		T TOTAL
-----------------------------------------------------------------------*/
inline dcomp T_total(dcomp z, double sigma_11inf, dcvec z1a, dcvec z2a, dvec La, dvec mua, int ma, int na, ddcvec a, int m_not_a, int m_is_a)
{
	// Defining the variables
	dcomp tau_11, tau_12, T;

	// Get teh taus
	std::tie(tau_11, tau_12) = tau_total(z, sigma_11inf, z1a, z2a, La, mua, ma, na, a, m_not_a);

	// Calculate the tractions
	T = dcomp(0, -.5) * (tau_11 * exp(dcomp(0, 2) * mua[m_is_a]) - tau_12);

	return { T };
}
/* --------------------------------------------------------------------
		SOLVE CRACK ANALYTIC ELEMENT A (BETA IN PAPER)
-----------------------------------------------------------------------*/
inline dcvec AE_crack_solver(Eigen::VectorXd T_s, Eigen::VectorXd T_n, Eigen::MatrixXd A, int ma, int Na, dvec pa, double sigma_11inf, dcvec z1a, dcvec z2a, dvec La, dvec mua, int na, ddcvec a_in, dcvec zint, ivec int_check, int m_is_a)
{
	// Defining the variables
	dcvec a(ma);
	Eigen::VectorXd b1(ma);
	Eigen::VectorXd b2(ma);

	// Count number of intersections
	int int_count = std::accumulate(int_check.begin(), int_check.end(), 0.0);

	// Check if the crack has any intersection
	if (int_count == 0)
	{
		// Solving the linear system (without intersections)
		b1 = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(T_s);
		b2 = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(T_n);
	}
	else
	{
		// Define new matrices for A and T and assign the pre-calculated values
		Eigen::MatrixXd A_int_s(Na + int_count * 3, ma);
		Eigen::MatrixXd A_int_n(Na + int_count, ma);
		Eigen::VectorXd T_s_int(Na + int_count * 3);
		Eigen::VectorXd T_n_int(Na + int_count);
		#pragma omp parallel for default(none) shared(A_int_s, A_int_n, A, T_n_int, T_n, T_s_int, T_s)
		for (int ii = 0; ii < Na; ii++)
		{
			for (int jj = 0; jj < ma; jj++)
			{
				/*A_int_re(ii, jj) = A(ii, jj);
				A_int_im(ii, jj) = A(ii, jj);*/
				A_int_s(ii, jj) = A(ii, jj);
				A_int_n(ii, jj) = A(ii, jj);
			}
			T_s_int(ii) = T_s(ii);
			T_n_int(ii) = T_n(ii);
		}
		int cnt_s = 0;
		int cnt_n = 0;
		for (int ii = 0; ii < na; ii++)
		{
			// Check if element ii has any intersections
			if (int_check[ii] == 1)
			{
				// Define variables
				dcomp chia, tau11, chi_pow, L_frac, T_temp;
				double thetaa;

				// Compute the intersection angel in the chi-plane
				chia = chi_from_z(zint[ii], z1a[m_is_a], z2a[m_is_a], La[m_is_a], mua[m_is_a]);
				thetaa = imag(log(chia));

				// Pre-calcualte terms
				chi_pow = chia /( chia * chia - 1.0 );
				L_frac = (dcomp(0,8.0) / La[m_is_a]) * exp(dcomp(0, -2) * mua[m_is_a]);

				// Add the intersection to the A matrix
				for (int jj = 0; jj < ma; jj++)
				{
					double n = jj + 1.0;
					// Add a control point at the intersection - tau^11 condition
					A_int_s(Na + cnt_s, jj) = real( chi_pow * L_frac * n * pow(chia, -n) );
					A_int_s(Na + (cnt_s + 1), jj) = imag( chi_pow * L_frac * n * pow(chia, -n) );
					
					// Add a control point at the intersection - traction condition
					A_int_s(Na + (cnt_s + 2), jj) = n * sin(n * thetaa);
					A_int_n(Na + cnt_n, jj) = n * sin(n * thetaa);
				}

				// Calculate the taus at the intersection
				dcomp tau_11, tau_12;
				std::tie(tau_11, tau_12) = tau_total(zint[ii], sigma_11inf, z1a, z2a, La, mua, ma, na, a_in, m_is_a);

				// Add the intersection to the T_s vector
				T_s_int(Na + cnt_s) = -real(tau_11);
				T_s_int(Na + (cnt_s + 1)) = -imag(tau_11);

				// Calculate the traction at the control point at the intersection
				T_temp = T_total(zint[ii], sigma_11inf, z1a, z2a, La, mua, ma, na, a_in, m_is_a, m_is_a);

				// Add the intersection to the T vectors
				T_s_int(Na + (cnt_s + 2)) = real(T_temp) * 0.5 * La[m_is_a] * sin(thetaa);
				T_n_int(Na + cnt_n) = (pa[ii] + imag(T_temp)) * -0.5 * La[m_is_a] * sin(thetaa);

				// Add step to loop count
				cnt_s += 3;
				cnt_n += 1;
			}
			
		}

		// Solving the linear system (with intersections)
		b1 = A_int_s.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(T_s_int);
		b2 = A_int_n.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(T_n_int);
	}
	// Assign to a
	for (int ii = 0; ii < ma; ii++)
	{
		a[ii] = dcomp(b2[ii], b1[ii]);
	}


	return { a };
}
/* --------------------------------------------------------------------
		ITERATOR
-----------------------------------------------------------------------*/
inline ddcvec iterator(double cond, int ITR, int Na, dvec pa, double sigma_11inf, dcvec z1a, dcvec z2a, dvec La, dvec mua, int ma, int na)
{
	// Defining the variables
	dcomp tau_11, tau_12, T;
	Eigen::MatrixXd A(Na, ma);
	dvec theta_a(Na);
	ddvec term(na, dvec(Na));
	ddcvec za(na, dcvec(Na));
	ddcvec a(na, dcvec(ma));
	ddcvec zint(na, dcvec(na));
	iivec int_check(na, ivec(na));
	ivec int_count(na);
	double theta_0a, delthetaa;

	// Print the conditions
	std::cout << std::scientific;
	std::cout << "Solver for " << na << " analytical element cracks with " << ma << " coefficients at " << Na << " integration points." << std::endl;
	std::cout << "Iterations break after: error < " << cond << " or iterations > " << ITR << "." << std::endl << std::endl;
	std::cout << "	Error:" << "		Iteration:" << std::endl;


	// Assigning variables
	theta_0a = pi() / (Na);
	delthetaa = (pi() - 2.5*theta_0a) / (Na - 1.0);

	if (na > 0)
	{
		// Find the intersection points
		#pragma omp parallel for default(none) shared(zint, int_check, int_count, z1a, z2a)
		for (int ii = 0; ii < na; ii++)
		{
			for (int jj = 0; jj < na; jj++)
			{
				std::tie(zint[ii][jj], int_check[ii][jj]) = intersection_point(z1a[ii], z2a[ii], z1a[jj], z2a[jj]);
			}
			int_count[ii] = std::accumulate(int_check[ii].begin(), int_check[ii].end(), 0.0);
		}

		// Calculating the A, a theta and term matrices
		#pragma omp parallel for default(none) shared(theta_a, theta_0a, delthetaa)
		for (int ii = 0; ii < Na; ii++)
		{
			theta_a[ii] = theta_0a + ii * delthetaa;
		}

		//#pragma omp parallel for default(none) shared(A, za, term, z1a, z2a, La, mua, theta_a)
		for (int ii = 0; ii < Na; ii++)
		{
			for (int mm = 0; mm < ma; mm++)
			{
				A(ii, mm) = (mm + 1.0) * sin((mm + 1.0) * theta_a[ii]);
			}
			for (int jj = 0; jj < na; jj++)
			{
				dcomp chi;
				chi = exp(dcomp(0, 1) * theta_a[ii]);
				za[jj][ii] = z_from_chi(chi, z1a[jj], z2a[jj], La[jj], mua[jj]);
				term[jj][ii] = -0.5 * La[jj] * sin(theta_a[ii]);
			}
		}
		
	}

	// Sovle the a (betas in paper)
	double error = 1;
	ddcvec a_old(na, dcvec(ma, 0));
	double error_a;
	int NIT = 0;
	while ( error > cond && NIT < ITR )
	{
				for (int ii = 0; ii < na; ii++)
				{
					Eigen::VectorXd T_s(Na);
					Eigen::VectorXd T_n(Na);
					#pragma omp parallel for default(none) shared(T_s, T_n, za, sigma_11inf, z1a, z2a, La, mua, ma, na, a)
					for (int jj = 0; jj < Na; jj++)
					{
						dcomp T;
						T = T_total(za[ii][jj], sigma_11inf, z1a, z2a, La, mua, ma, na, a, ii, ii);
						T_s(jj) = -real(T) * term[ii][jj];
						T_n(jj) = (pa[ii] + imag(T)) * term[ii][jj];
					}
					dcomp tau_11, tau_12;
					a[ii] = AE_crack_solver(T_s, T_n, A, ma, Na, pa, sigma_11inf, z1a, z2a, La, mua, na, a, zint[ii], int_check[ii], ii);
					std::tie(tau_11, tau_12) = tau_total(zint[0][1], sigma_11inf, z1a, z2a, La, mua, ma, na, a, -1);
				}

		// Calcualte the error
		error_a = 0;
		if (na > 0)
		{
			dvec dela(na);
			#pragma omp parallel for default(none) shared(dela, a, error_a)
			for (int ii = 0; ii < na; ii++)
			{
				dvec dela_temp(ma);
				for (int jj = 0; jj < ma; jj++)
				{
					dela_temp[jj] = abs(a[ii][jj] - a_old[ii][jj]);
				}
				dela[ii] = *max_element(dela_temp.begin(), dela_temp.end());
			}
			error_a = *max_element(dela.begin(), dela.end());
		}
		dvec error_vec{ error_a };
		error = *max_element(error_vec.begin(), error_vec.end());

		NIT += 1;
		a_old = a;

		std::cout << std::scientific;
		std::cout << "	" << error << "	" << NIT << std::endl;

	}

	return { a };
}
/* --------------------------------------------------------------------
		STRESS TOTAL
-----------------------------------------------------------------------*/
inline std::tuple<double, double, double> sigma_total(dcomp z, double sigma_11inf, dcvec z1a, dcvec z2a, dvec La, dvec mua, int ma, int na, ddcvec a)
{
	// Defining the variables
	dcomp S1, S2, tau_11, tau_12;
	double sigma_11, sigma_22, sigma_12;

	// Calculating the tau
	std::tie(tau_11, tau_12) = tau_total(z, sigma_11inf, z1a, z2a, La, mua, ma, na, a, -1);

	// Calculate the sigmas
	S1 = .5 * (tau_11 + tau_12);
	S2 = .5 * (-tau_11 + tau_12);
	sigma_11 = real(S1);
	sigma_22 = real(S2);
	sigma_12 = -imag(S1);

	return { sigma_11, sigma_22, sigma_12 };
}
/* --------------------------------------------------------------------
		STRESS FIELD
-----------------------------------------------------------------------*/
std::tuple<dvec, dvec, ddvec, ddvec, ddvec> stress_field(double xfrom, double xto, double yfrom, double yto, int Nx, int Ny, double sigma_11inf, dcvec z1a, dcvec z2a, dvec La, dvec mua, int ma, int na, ddcvec a)
{
	// Defining the variables
	ddvec grid_11(Nx, dvec(Ny));
	ddvec grid_22(Nx, dvec(Ny));
	ddvec grid_12(Nx, dvec(Ny));
	double dx;
	double dy;
	dvec x_vec(Nx);
	dvec y_vec(Ny);

	// Calcualte the sigma grids
	dx = (xto - xfrom) / (Nx - 1.0);
	dy = (yto - yfrom) / (Ny - 1.0);
	#pragma omp parallel for default(none) shared(grid_11, grid_22, grid_12, x_vec, y_vec, sigma_11inf, z1a, z2a, La, mua, ma, na, a)
	for (int ii = 0; ii < Nx; ii++)
	{
		for (int jj = Ny; jj--;)
		{
			x_vec[ii] = xfrom + ii * dx;
			y_vec[jj] = yfrom + jj * dy;
			std::tie(grid_11[ii][jj], grid_22[ii][jj], grid_12[ii][jj]) = sigma_total(dcomp(x_vec[ii], y_vec[jj]), sigma_11inf, z1a, z2a, La, mua, ma, na, a);
		}
	}
	return { x_vec, y_vec, grid_11, grid_22, grid_12 };
}
/* --------------------------------------------------------------------
		PRINCIPAL STRESSES
-----------------------------------------------------------------------*/
std::tuple<double, double, double> principal_sigma(dcomp z, double sigma_11inf, dcvec z1a, dcvec z2a, dvec La, dvec mua, int ma, int na, ddcvec a)
{
	// Defining the variables
	double sigma_1;
	double sigma_2;
	double theta_p;
	double frac1, frac2, sqrt1;
	dcomp S1, S2, tau_11, tau_12;
	double sigma_11, sigma_22, sigma_12;

	// Calculating the tau
	std::tie(tau_11, tau_12) = tau_total(z, sigma_11inf, z1a, z2a, La, mua, ma, na, a, -1);

	// Calculate the sigmas
	S1 = .5 * (tau_11 + tau_12);
	S2 = .5 * (-tau_11 + tau_12);
	sigma_11 = real(S1);
	sigma_22 = real(S2);
	sigma_12 = -imag(S1);

	// Calculating the terms
	frac1 = (sigma_11 + sigma_22) / 2.0;
	frac2 = (sigma_11 - sigma_22) / 2.0;
	sqrt1 = sqrt(frac2 * frac2 + sigma_12 * sigma_12);

	// Calculating the principal stresses and the angel of sigma
	sigma_1 = frac1 + sqrt1;
	sigma_2 = frac1 - sqrt1;

	// Calcuating the angel theta_p
	theta_p = -0.5 * imag(log(tau_11));

	// Changing the absolut sigma is sigma_2 > sigma_1 (TURNED OFF)
	/*if (abs(sigma_2) > abs(sigma_1))
	{
		sigma_1 = frac1 - sqrt1;
		sigma_2 = frac1 + sqrt1;
		theta_p = -0.5*imag(log(tau_11)) + pi() / 2;
	}
	*/



	return { sigma_1, sigma_2, theta_p };
}
/* --------------------------------------------------------------------
		PRINCIPAL STRESSES PLOT
-----------------------------------------------------------------------*/
std::tuple<double, double, double> principal_sigma_plt(double sigma_11, double sigma_22, double sigma_12)
{
	// Defining the variables
	double sigma_1;
	double sigma_2;
	double theta_p;
	double frac1, frac2, sqrt1;
	dcomp tau_11;

	// Calculating the tau
	tau_11 = sigma_11 - sigma_22 - dcomp(0, 2) * sigma_12;

	// Calculating the terms
	frac1 = (sigma_11 + sigma_22) / 2.0;
	frac2 = (sigma_11 - sigma_22) / 2.0;
	sqrt1 = sqrt(frac2 * frac2 + sigma_12 * sigma_12);

	// Calculating the principal stresses and the angel of sigma
	sigma_1 = frac1 + sqrt1;
	sigma_2 = frac1 - sqrt1;

	// Calcuating the angel theta_p
	theta_p = -0.5 * imag(log(tau_11));

	// Changing the absolut sigma is sigma_2 > sigma_1 (TURNED OFF)
	/*if (abs(sigma_2) > abs(sigma_1))
	{
		sigma_1 = frac1 - sqrt1;
		sigma_2 = frac1 + sqrt1;
		theta_p = -0.5*imag(log(tau_11)) + pi() / 2;
	}
	*/



	return { sigma_1, sigma_2, theta_p };
}
/* --------------------------------------------------------------------
		PRINCIPAL STRESS FIELDS
-----------------------------------------------------------------------*/
std::tuple<ddvec, ddvec, ddvec> principal_stress_field(int Nx, int Ny, ddvec grid_11, ddvec grid_22, ddvec grid_12)
{
	// Defining the variables
	ddvec grid_1(Nx, dvec(Ny));
	ddvec grid_2(Nx, dvec(Ny));
	ddvec grid_tp(Nx, dvec(Ny));

	// Calculate teh principal stresses from the stress feilds
	#pragma omp parallel for default(none) shared(grid_1, grid_2, grid_tp, grid_11, grid_22, grid_12)
	for (int ii = 0; ii < Nx; ii++)
	{
		for (int jj = Ny; jj--;)
		{
			std::tie(grid_1[ii][jj], grid_2[ii][jj], grid_tp[ii][jj]) = principal_sigma_plt(grid_11[ii][jj], grid_22[ii][jj], grid_12[ii][jj]);
		}
	}
	return { grid_1, grid_2, grid_tp };
}
/* --------------------------------------------------------------------
		PRINCIPAL STRESS TRAJECTORIES
-----------------------------------------------------------------------*/
std::tuple<ddcvec, ddcvec> principal_stress_trajectories(double xfrom, double xto, double yfrom, double yto, dcvec xtraj, dcvec ytraj, int Ntraj, int lvs_traj, double sigma_11inf, dcvec z1a, dcvec z2a, dvec La, dvec mua, int ma, int na, ddcvec a)
{
	// Defining the variables
	ddcvec traj_1(lvs_traj * 2, dcvec(Ntraj));
	ddcvec traj_2(lvs_traj * 2, dcvec(Ntraj));
	double dx, dy, dx_lvsre, dx_lvsim, dy_lvsre, dy_lvsim, pi_val;
	double cond;
	int NIT;

	// Getting the starting points
	dx = (xto - xfrom) / (Ntraj * 0.8);
	dy = (yto - yfrom) / (Ntraj * 0.8);
	dx_lvsre = (real(xtraj[1]) - real(xtraj[0])) / (lvs_traj - 1.0);
	dx_lvsim = (imag(xtraj[1]) - imag(xtraj[0])) / (lvs_traj - 1.0);
	dy_lvsre = (real(ytraj[1]) - real(ytraj[0])) / (lvs_traj - 1.0);
	dy_lvsim = (imag(ytraj[1]) - imag(ytraj[0])) / (lvs_traj - 1.0);
	pi_val = 0.5 * pi();
	cond = 1e-6;
	NIT = 10;

	// SIGMA 1
	#pragma omp parallel for default(none) shared(traj_1, sigma_11inf, z1a, z2a, La, mua, ma, na, a)
	for (int ii = 0; ii < lvs_traj; ii++)
	{
		traj_1[ii][0] = dcomp(real(xtraj[0]) + ii * dx_lvsre, imag(xtraj[0]) + ii * dx_lvsim);
		dcomp z = traj_1[ii][0];
		dcomp z_old = z;
		double dx1 = dx;
		for (int jj = 1; jj < Ntraj; jj++)
		{
			dcomp zt, z11, z_oldt, sigma, sigmat;
			double sigma_1, theta_p;
			double ee, eta;
			std::tie(sigma_1, std::ignore, theta_p) = principal_sigma(z, sigma_11inf, z1a, z2a, La, mua, ma, na, a);
			sigma = abs(sigma_1) * exp(dcomp(0, 1) * theta_p);
			z11 = z + sigma / abs(sigma) * dx1;
			eta = angel_change(z11, z, z_old);
			if (eta > pi_val && jj > 1) {
				dx1 = -dx1;
			}
			zt = z + sigma / abs(sigma) * dx1;

			ee = 1;
			for (int rr = NIT; rr--;)
			{
				z_oldt = zt;
				std::tie(sigma_1, std::ignore, theta_p) = principal_sigma(zt, sigma_11inf, z1a, z2a, La, mua, ma, na, a);
				sigmat = abs(sigma_1) * exp(dcomp(0, 1) * theta_p);
				z11 = z + (sigma + sigmat) / abs(sigma + sigmat) * dx1;
				eta = angel_change(z11, z, z_old);
				if (eta > pi_val && jj > 1) {
					dx1 = -dx1;
				}
				zt = z + (sigma + sigmat) / abs(sigma + sigmat) * dx1;
				ee = std::norm(z_oldt - zt);
				if (ee < cond)
				{
					break;
				}
			}

			traj_1[ii][jj] = zt;
			z_old = z;
			z = zt;
		}

		int kk = ii + lvs_traj;
		traj_1[kk][0] = traj_1[ii][0];
		z = traj_1[kk][0];
		z_old = z;
		dx1 = -dx;
		for (int jj = 1; jj < Ntraj; jj++)
		{
			dcomp zt, z11, z_oldt, sigma, sigmat;
			double sigma_1, theta_p;
			double ee, eta;
			std::tie(sigma_1, std::ignore, theta_p) = principal_sigma(z, sigma_11inf, z1a, z2a, La, mua, ma, na, a);
			sigma = abs(sigma_1) * exp(dcomp(0, 1) * theta_p);
			z11 = z + sigma / abs(sigma) * dx1;
			eta = angel_change(z11, z, z_old);
			if (eta > pi_val && jj > 1) {
				dx1 = -dx1;
			}
			zt = z + sigma / abs(sigma) * dx1;

			ee = 1;
			for (int rr = NIT; rr--;)
			{
				z_oldt = zt;
				std::tie(sigma_1, std::ignore, theta_p) = principal_sigma(zt, sigma_11inf, z1a, z2a, La, mua, ma, na, a);
				sigmat = abs(sigma_1) * exp(dcomp(0, 1) * theta_p);
				z11 = z + (sigma + sigmat) / abs(sigma + sigmat) * dx1;
				eta = angel_change(z11, z, z_old);
				if (eta > pi_val && jj > 1) {
					dx1 = -dx1;
				}
				zt = z + (sigma + sigmat) / abs(sigma + sigmat) * dx1;
				ee = std::norm(z_oldt - zt);
				if (ee < cond)
				{
					break;
				}
			}

			traj_1[kk][jj] = zt;
			z_old = z;
			z = zt;
		}
	}

	// SIGMA 2
	#pragma omp parallel for default(none) shared(traj_2, sigma_11inf, z1a, z2a, La, mua, ma, na, a)
	for (int ii = 0; ii < lvs_traj; ii++)
	{
		traj_2[ii][0] = dcomp(real(ytraj[0]) + ii * dy_lvsre, imag(ytraj[0]) + ii * dy_lvsim);
		double dy1 = dy;
		dcomp z = traj_2[ii][0];
		dcomp z_old = z;
		for (int jj = 1; jj < Ntraj; jj++)
		{
			dcomp zt, z11, z_oldt, sigma, sigmat;
			double sigma_2, theta_p;
			double ee, eta;
			std::tie(std::ignore, sigma_2, theta_p) = principal_sigma(z, sigma_11inf, z1a, z2a, La, mua, ma, na, a);
			sigma = abs(sigma_2) * exp(dcomp(0, 1) * (theta_p + pi_val));
			z11 = z + sigma / abs(sigma) * dy1;
			eta = angel_change(z11, z, z_old);
			if (eta > pi_val && jj > 1) {
				dy1 = -dy1;
			}
			zt = z + sigma / abs(sigma) * dy1;

			ee = 1;
			for (int rr = NIT; rr--;)
			{
				z_oldt = zt;
				std::tie(std::ignore, sigma_2, theta_p) = principal_sigma(zt, sigma_11inf, z1a, z2a, La, mua, ma, na, a);
				sigmat = abs(sigma_2) * exp(dcomp(0, 1) * (theta_p + pi_val));
				z11 = z + (sigma + sigmat) / abs(sigma + sigmat) * dy1;
				eta = angel_change(z11, z, z_old);
				if (eta > pi_val && jj > 1) {
					dy1 = -dy1;
				}
				zt = z + (sigma + sigmat) / abs(sigma + sigmat) * dy1;
				ee = std::norm(z_oldt - zt);
				if (ee < cond)
				{
					break;
				}
			}

			traj_2[ii][jj] = zt;
			z_old = traj_2[ii][(int)(jj - 1)];
			z = zt;
		}

		int kk = ii + lvs_traj;
		traj_2[kk][0] = traj_2[ii][0];
		dy1 = -dy;
		z = traj_2[kk][0];
		z_old = z;
		for (int jj = 1; jj < Ntraj; jj++)
		{
			dcomp zt, z11, z_oldt, sigma, sigmat;
			double sigma_2, theta_p;
			double ee, eta;
			std::tie(std::ignore, sigma_2, theta_p) = principal_sigma(z, sigma_11inf, z1a, z2a, La, mua, ma, na, a);
			sigma = abs(sigma_2) * exp(dcomp(0, 1) * (theta_p + pi_val));
			z11 = z + sigma / abs(sigma) * dy1;
			eta = angel_change(z11, z, z_old);
			if (eta > pi_val && jj > 1) {
				dy1 = -dy1;
			}
			zt = z + sigma / abs(sigma) * dy1;

			ee = 1;
			for (int rr = NIT; rr--;)
			{
				z_oldt = zt;
				std::tie(std::ignore, sigma_2, theta_p) = principal_sigma(zt, sigma_11inf, z1a, z2a, La, mua, ma, na, a);
				sigmat = abs(sigma_2) * exp(dcomp(0, 1) * (theta_p + pi_val));
				z11 = z + (sigma + sigmat) / abs(sigma) * dy1;
				eta = angel_change(z11, z, z_old);
				if (eta > pi_val && jj > 1) {
					dy1 = -dy1;
				}
				zt = z + (sigma + sigmat) / abs(sigma + sigmat) * dy1;
				ee = std::norm(z_oldt - zt);
				if (ee < cond)
				{
					break;
				}
			}

			traj_2[kk][jj] = zt;
			z_old = z;
			z = zt;
		}
	}
	return { traj_1, traj_2 };
}
/* --------------------------------------------------------------------
		W UNIFORM STRESS
-----------------------------------------------------------------------*/
inline dcomp  w_uni(dcomp z, double kappa, double G, double sigma_11inf)
{
	// Defining the variables
	dcomp phi_bar, dphi, psi, w;

	// calculating the veriables
	phi_bar = -0.5 * sigma_11inf * conj(z);
	dphi = -0.5 * sigma_11inf;
	psi = -0.5 * sigma_11inf * z;

	// Calculating w
	w = 1 / (4 * G) * ((z - conj(z)) * dphi + kappa * phi_bar + psi);

	return { w };
}
/* --------------------------------------------------------------------
		W CRACK - ANALYTIC ELEMENT
-----------------------------------------------------------------------*/
inline dcomp w_crack(dcomp z, double kappa, double G, dcomp z1, dcomp z2, double L, double mu, int m, dcvec a)
{
	// Defining the variables
	dcomp chi, chi_bar, Z, chi_pow;
	dcomp dphi, phi_bar, psi;
	dcomp w, L_frac;
	double n;

	// Getting the chi - and Z - coordinates
	chi = chi_from_z(z, z1, z2, L, mu);
	chi_bar = conj(chi);
	Z = exp(dcomp(0, -1) * mu) * 2.0 * (z - 0.5 * (z1 + z2)) / L;

	// Calculating the series
	phi_bar = 0;
	dphi = 0;
	psi = 0;
	n = 0;
	chi_pow = chi * chi - 1.0;
	for (int ii = 0; ii < m; ii++)
	{
		dcomp a_n;
		n += 1;
		a_n = a[ii] * n;
		dphi += conj(a_n) * pow(chi, (1.0 - n)) / chi_pow;
		phi_bar -= a[ii] * pow(chi_bar, -n);
		psi += a[ii] * pow(chi, -n);
	}

	// Multiplying the constants
	L_frac = (4.0 / L) * exp(dcomp(0, -1) * mu);
	dphi *= L_frac;

	// Calcualting w
	w = 1 / (4 * G) * (0.5 * L * (Z - conj(Z)) * exp(dcomp(0, 1) * mu) * dphi + kappa * phi_bar + psi);

	return { w };
}
/* --------------------------------------------------------------------
		W TOTAL
-----------------------------------------------------------------------*/
inline dcomp  w_total(dcomp z, double kappa, double G, double sigma_11inf, dcvec z1a, dcvec z2a, dvec La, dvec mua, int ma, int na, ddcvec a)
{
	// Defining the variables
	dcomp w, wg;

	// Add the unfirm stress field
	w = w_uni(z, kappa, G, sigma_11inf);

	// Add the anlytic element for a crack
	if (na > 0)
	{
		for (int ii = 0; ii < na; ii++)
		{
			dcomp wc;
			wc = w_crack(z, kappa, G, z1a[ii], z2a[ii], La[ii], mua[ii], ma, a[ii]);
			w += wc * exp(dcomp(0, -1) * mua[ii]);
		}
	}

	return { w };
}
/* --------------------------------------------------------------------
		DISPLACEMENT FIELD
-----------------------------------------------------------------------*/
std::tuple<dvec, dvec, ddcvec> w_field(double xfrom, double xto, double yfrom, double yto, int Nw, double kappa, double G, double sigma_11inf, dcvec z1a, dcvec z2a, dvec La, dvec mua, int ma, int na, ddcvec a)

{
	// Defining the variables
	ddcvec grid_w(Nw, dcvec(Nw));
	dvec x_vecw(Nw), y_vecw(Nw);
	double dx;
	double dy;

	// Calcualte the displacement grid
	dx = (xto - xfrom) / (Nw - 1.0);
	dy = (yto - yfrom) / (Nw - 1.0);
#pragma omp parallel for default(none) shared(grid_w, x_vecw, y_vecw, kappa, G, sigma_11inf, z1a, z2a, La, mua, ma, na, a)
	for (int ii = 0; ii < Nw; ii++)
	{
		for (int jj = 0; jj < Nw; jj++)
		{

			x_vecw[ii] = xfrom + ii * dx;
			y_vecw[jj] = yfrom + jj * dy;
			grid_w[ii][jj] = w_total(dcomp(x_vecw[ii], y_vecw[jj]), kappa, G, sigma_11inf, z1a, z2a, La, mua, ma, na, a);
		}
	}
	return { x_vecw, y_vecw, grid_w };
}
/*---------------------------------------------------------------------
		DISPLACEMENT TRAJECTORIES
-----------------------------------------------------------------------*/
ddcvec w_trajectories(double xfrom, double xto, double yfrom, double yto, int Ntraj, int Nw, double kappa, double G, double sigma_11inf, dcvec z1a, dcvec z2a, dvec La, dvec mua, int ma, int na, ddcvec a)
{
	// Defining the variables
	ddcvec traj_w((int)(Nw * 2), dcvec(Ntraj));
	double dx, dy_lvs;
	double cond, pi_val;
	int NIT;

	// Getting the starting points
	dx = (xto - xfrom) / (Ntraj);
	dy_lvs = (yto - yfrom) / (Nw - 1.0);
	pi_val = 0.5 * pi();
	cond = 1e-6;
	NIT = 10;

	// w trajectories
#pragma omp parallel for default(none) shared(traj_w, kappa, G, sigma_11inf, z1a, z2a, La, mua, ma, na, a)
	for (int ii = 0; ii < Nw; ii++)
	{
		traj_w[ii][0] = dcomp(xfrom, yfrom + ii * dy_lvs);
		dcomp z = traj_w[ii][0];
		dcomp z_old = z;
		double dx1 = dx;
		for (int jj = 1; jj < Ntraj; jj++)
		{
			dcomp zt, z_oldt, w, w1;
			double ee;
			w = w_total(z, kappa, G, sigma_11inf, z1a, z2a, La, mua, ma, na, a);
			zt = z + conj(w) / abs(w) * dx1;

			ee = 1;
			for (int rr = NIT; rr--;)
			{
				z_oldt = zt;
				w1 = w_total(zt, kappa, G, sigma_11inf, z1a, z2a, La, mua, ma, na, a);
				zt = z + conj(w + w1) / abs(w + w1) * dx1;
				ee = std::norm(z_oldt - zt);
				if (ee < cond)
				{
					break;
				}
			}

			traj_w[ii][jj] = zt;
			z_old = z;
			z = zt;
		}

		int kk = ii + Nw;
		traj_w[kk][0] = dcomp(xto, yfrom + ii * dy_lvs);
		z = traj_w[kk][0];
		z_old = z;
		dx1 = dx;
		for (int jj = 1; jj < Ntraj; jj++)
		{
			dcomp zt, z_oldt, w, w1;
			double ee;
			w = w_total(z, kappa, G, sigma_11inf, z1a, z2a, La, mua, ma, na, a);
			zt = z + conj(w) / abs(w) * dx1;

			ee = 1;
			for (int rr = NIT; rr--;)
			{
				z_oldt = zt;
				w1 = w_total(zt, kappa, G, sigma_11inf, z1a, z2a, La, mua, ma, na, a);
				zt = z + conj(w + w1) / abs(w + w1) * dx1;
				ee = std::norm(z_oldt - zt);
				if (ee < cond)
				{
					break;
				}
			}

			traj_w[kk][jj] = zt;
			z_old = z;
			z = zt;

		}
	}
	return { traj_w };
}

/* --------------------------------------------------------------------

		MAIN SCRIPT

-----------------------------------------------------------------------*/
int main()
{
	/* --------------------------------------------------------------------
			Defining variables and importing data
	-----------------------------------------------------------------------*/

	// Header in console window
	auto start = std::chrono::high_resolution_clock::now(); // Start the clock

	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "		ANALYTIC ELEMENT LINEAR ELASTIC SOLVER	" << std::endl << std::endl;
	std::cout << "=================================================================" << std::endl << std::endl;

	auto date_time1 = std::chrono::system_clock::now();
	std::time_t start_time = std::chrono::system_clock::to_time_t(date_time1);
	char str_time1[26];
	ctime_s(str_time1, sizeof str_time1, &start_time);
	std::cout << "Program started: " << str_time1 << std::endl;

	// Setting the data types
	dcomp z;
	double kappa, sigma_11inf, G, cond;
	int na, ma, Na, NIT;

	// Read the input data from binary file PART 1
	std::ifstream input_file("geometry_data.bin", std::ios::in | std::ios::binary | std::ios::ate);
	std::streampos size = input_file.tellg();
	char* memblock = new char[size];
	input_file.seekg(0, std::ios::beg);
	input_file.read(memblock, size);
	double* fin = (double*)memblock;//reinterpret as doubles

	sigma_11inf = fin[0];
	kappa = fin[1];
	G = fin[2];
	na = (int)fin[3];
	ma = (int)fin[4];
	Na = (int)fin[5];
	cond = (double)fin[6];
	NIT = (int)fin[7];

	// Declaring the vecotrs
	dcvec z1a(na), z2a(na);
	dvec pa(na), La(na), mua(na);
	ddcvec a(na, dcvec(ma));

	int pos = 7 + 1;
	if (na > 0)
	{
		for (int ii = 0; ii < na; ii++)
		{
			int re = pos + ii;
			int im = pos + na + ii;
			z1a[ii] = dcomp(fin[re], fin[im]);
		}
		pos += 2 * na;
		for (int ii = 0; ii < na; ii++)
		{
			int re = pos + ii;
			int im = pos + na + ii;
			z2a[ii] = dcomp(fin[re], fin[im]);
		}
		pos += 2 * na;
		for (int ii = 0; ii < na; ii++)
		{
			pa[ii] = fin[pos + ii];
		}
		pos += na;
		for (int ii = 0; ii < na; ii++)
		{
			La[ii] = fin[pos + ii];
		}
		pos += na;
		for (int ii = 0; ii < na; ii++)
		{
			mua[ii] = fin[pos + ii];
		}
		pos += na;
	}
	else
	{
		z1a = { dcomp(0, 0) };
		z2a = { dcomp(0, 0) };
		La = { 0 };
		mua = { 0 };
		a = { { dcomp(0,0) } };
		ma = 1;
		Na = 1;
	}

	// Setting the data types
	double xfrom, xto, yfrom, yto;
	int Nx, Ny, Nw, Ntraj, lvs_traj;
	ddvec grid_11, grid_22, grid_12, grid_1, grid_2, theta_p;
	ddcvec traj_1, traj_2, grid_w, traj_w;
	dvec x_vec, y_vec, x_vecw, y_vecw;
	dcvec xtraj, ytraj;

	// Read the plot data from binary file
	std::ifstream plot_file("plot_data.bin", std::ios::in | std::ios::binary | std::ios::ate);
	std::streampos size2 = plot_file.tellg();
	char* memblock2 = new char[size2];
	plot_file.seekg(0, std::ios::beg);
	plot_file.read(memblock2, size2);
	double* fplot = (double*)memblock2;//reinterpret as doubles

	xfrom = fplot[0];
	xto = fplot[1];
	yfrom = fplot[2];
	yto = fplot[3];
	Nx = (int)fplot[4];
	Ny = (int)fplot[5];
	Nw = (int)fplot[6];
	Ntraj = (int)fplot[7];
	lvs_traj = (int)fplot[8];
	xtraj = { fplot[9] + dcomp(0,1) * fplot[10], fplot[11] + dcomp(0,1) * fplot[12] };
	ytraj = { fplot[13] + dcomp(0,1) * fplot[14], fplot[15] + dcomp(0,1) * fplot[16] };

	// Displying the plot data in the command window
	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "		THE GEOMETRY DATA	" << std::endl << std::endl;
	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "This is the retrived geometry data:\n";
	std::cout << "sigma_11inf = " << sigma_11inf << std::endl;
	std::cout << "      kappa = " << kappa << std::endl;
	std::cout << "          G = " << G << std::endl;
	std::cout << "         na = " << na << std::endl;
	std::cout << "         ma = " << ma << std::endl;
	std::cout << "         Na = " << Na << std::endl;
	// Displying the plot data in the command window
	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "		THE READ PLOT DATA	" << std::endl << std::endl;
	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "This is the retrived plot data:\n";
	std::cout << "x from: " << xfrom << " to " << xto << std::endl;
	std::cout << "y from: " << yfrom << " to " << yto << std::endl;
	std::cout << "x resolution: " << Nx << std::endl;
	std::cout << "y resolution: " << Ny << std::endl;
	std::cout << "Total number of points: " << Nx * Ny << std::endl;
	std::cout << "Number of steps in trajectories: " << Ntraj << std::endl;
	std::cout << "Number of trajectory levels: " << lvs_traj << std::endl;
	std::cout << "Total number of trajectory points: " << Ntraj * lvs_traj * 4 << std::endl;
	std::cout << "sigma_1 starting line from: " << xtraj[0] << " to " << xtraj[1] << std::endl;
	std::cout << "sigma_2 starting line from: " << ytraj[0] << " to " << ytraj[1] << std::endl;
	std::cout << "x and y quiver resolution: " << Nw << std::endl << std::endl;
	std::cout << std::endl;

	/* --------------------------------------------------------------------
			Solve the system
	-----------------------------------------------------------------------*/
	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "		INITIALIZING THE SOVLER	" << std::endl << std::endl;
	std::cout << "=================================================================" << std::endl << std::endl;

	auto start_solv = std::chrono::high_resolution_clock::now();
	a = iterator(cond, NIT, Na, pa, sigma_11inf, z1a, z2a, La, mua, ma, na);
	auto stop_solv = std::chrono::high_resolution_clock::now();
	auto duration_solv = std::chrono::duration_cast<std::chrono::microseconds>(stop_solv - start_solv);

	// Displaying the computation time
	std::cout << std::endl;
	long long ms_solv = duration_solv.count();
	long long s_solv, m_solv, h_solv;
	std::tie(ms_solv, s_solv, m_solv, h_solv) = ms_to_time(ms_solv);
	std::cout << "Computations finnished after ";
	time_print(ms_solv, s_solv, m_solv, h_solv);

	std::cout << std::endl;

	/* --------------------------------------------------------------------
			Checking the error
	-----------------------------------------------------------------------*/
	double error_med_a_re{}, error_mean_a_re{}, error_max_a_re{}, error_med_a_im{}, error_mean_a_im{}, error_max_a_im{}, error_med_int{}, error_mean_int{}, error_max_int{};
	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "		ERRORS	" << std::endl << std::endl;
	std::cout << "=================================================================" << std::endl << std::endl;
	int numa;
	numa = (int)round(Na);
	dcvec T_check(numa* na);
	if (na > 0)
	{
		// Calculate the error along the cracks, i.e. T = t_s + i*t_n = 0 + i*p
		dvec T_check_a_re(numa * na), T_check_a_im(numa * na);
		double theta_0a = pi()/numa;
		double delthetaa = (pi() - 2.0*theta_0a) / (numa-1.0);
		for (int ii = 0; ii < na; ii++)
		{
			for (int jj = 0; jj < numa; jj++)
			{
				dcomp z, tau_11, tau_12;
				double theta;
				theta = theta_0a + (jj) * delthetaa;
				z = z_from_chi(exp(dcomp(0, 1) * theta), z1a[ii], z2a[ii], La[ii], mua[ii]);
				T_check[(ii * numa) + jj] = T_total(z, sigma_11inf, z1a, z2a, La, mua, ma, na, a, -1, ii) + dcomp(0, pa[ii]);
				T_check_a_re[(ii * numa) + jj] = abs(real(T_check[(ii * numa) + jj]));
				T_check_a_im[(ii * numa) + jj] = abs(imag(T_check[(ii * numa) + jj]));
			}
		}
		// Median
		size_t n1 = T_check_a_re.size() / 2;
		nth_element(T_check_a_re.begin(), T_check_a_re.begin() + n1, T_check_a_re.end());
		error_med_a_re = T_check_a_re[n1];
		nth_element(T_check_a_im.begin(), T_check_a_im.begin() + n1, T_check_a_im.end());
		error_med_a_im = T_check_a_im[n1];
		// Mean
		error_mean_a_re = 1.0 * std::accumulate(T_check_a_re.begin(), T_check_a_re.end(), 0.0) / (n1 * 2);
		error_mean_a_im = 1.0 * std::accumulate(T_check_a_im.begin(), T_check_a_im.end(), 0.0) / (n1 * 2);
		// Max
		error_max_a_re = *max_element(T_check_a_re.begin(), T_check_a_re.end());
		error_max_a_im = *max_element(T_check_a_im.begin(), T_check_a_im.end());
		// Print the results
		std::cout << "     Analytic Element for a Crack" << std::endl;
		std::cout << "     Difference for ts" << std::endl;
		std::cout << "     Maximum = " << error_max_a_re << std::endl;
		std::cout << "     Mean    = " << error_mean_a_re << std::endl;
		std::cout << "     Median  = " << error_med_a_re << std::endl;
		std::cout << "     Difference for tn" << std::endl;
		std::cout << "     Maximum = " << error_max_a_im << std::endl;
		std::cout << "     Mean    = " << error_mean_a_im << std::endl;
		std::cout << "     Median  = " << error_med_a_im << std::endl;

		// Error at the intersection
		ddcvec zint(na, dcvec(na));
		iivec int_check(na, ivec(na));
		ivec int_count(na);

		// Find the intersection points
		#pragma omp parallel for default(none) shared(zint, int_check, int_count, z1a, z2a)
		for (int ii = 0; ii < na; ii++)
		{
			for (int jj = 0; jj < na; jj++)
			{
				std::tie(zint[ii][jj], int_check[ii][jj]) = intersection_point(z1a[ii], z2a[ii], z1a[jj], z2a[jj]);
			}
			int_count[ii] = std::accumulate(int_check[ii].begin(), int_check[ii].end(), 0.0);
		}
		int int_sum = std::accumulate(int_count.begin(), int_count.end(), 0.0);
		dvec int_error(int_sum);
		int cnt_int = 0;

		// Calcualte the error at the intersections, i.e. real(tau^11) = 0
		if (int_sum > 0)
		{
			for (int ii = 0; ii < na; ii++)
			{
				if (int_count[ii] > 0)
				{
					int cnt_er = 0;
					for (int jj = 0; jj < na; jj++)
					{
						if (int_check[ii][jj] == 1)
						{
							dcomp tau_11, tau_12, T;
							std::tie(tau_11, tau_12) = tau_total(zint[ii][jj], sigma_11inf, z1a, z2a, La, mua, ma, na, a, -1);
							T = T_total(z, sigma_11inf, z1a, z2a, La, mua, ma, na, a, -1, ii);
							std::cout << "z = " << zint[ii][jj] << " tau_11 = " << tau_11 << " tau_12 = " << tau_12 << " T = " << T << std::endl;
							int_error[cnt_int + cnt_er] = abs(real(tau_11));
							cnt_er += 1;
						}
					}
				}
				cnt_int += int_count[ii];
			}

			// Mean
			error_mean_int = std::accumulate(int_error.begin(), int_error.end(), 0.0) / (int_sum);
			// Max
			error_max_int = *max_element(int_error.begin(), int_error.end());
			// Median
			size_t n2 = int_error.size() / 2;
			nth_element(int_error.begin(), int_error.begin() + n2, int_error.end());
			error_med_int = int_error[n2];
			// Print the results
			std::cout << "     Intersection/s:" << std::endl;
			std::cout << "     Maximum = " << error_max_int << std::endl;
			std::cout << "     Mean    = " << error_mean_int << std::endl;
			std::cout << "     Median  = " << error_med_int << std::endl;
		}
	}

	// Create/open the two output files
	std::ofstream outfile_coef("input_data.bin", std::ios::out | std::ios::binary);

	// Save the plot properties
	dvec prop = { sigma_11inf,  kappa, G, 1.0 * na, 1.0 * ma};
	const char* pointerprop = reinterpret_cast<const char*>(&prop[0]);
	std::size_t bytesprop = prop.size() * sizeof(prop[0]);
	outfile_coef.write(pointerprop, bytesprop);


	// saving the coordinates of the cracks
	if (na > 0)
	{
		dvec fz1_re(na);
		dvec fz1_im(na);
		for (int jj = 0; jj < na; jj++)
		{
			fz1_re[jj] = real(z1a[jj]);
			fz1_im[jj] = imag(z1a[jj]);
		}
		const char* pointerz1_re = reinterpret_cast<const char*>(&fz1_re[0]);
		std::size_t bytesz1_re = fz1_re.size() * sizeof(fz1_re[0]);
		outfile_coef.write(pointerz1_re, bytesz1_re);
		const char* pointerz1_im = reinterpret_cast<const char*>(&fz1_im[0]);
		std::size_t bytesz1_im = fz1_im.size() * sizeof(fz1_im[0]);
		outfile_coef.write(pointerz1_im, bytesz1_im);

		dvec fz2_re(na);
		dvec fz2_im(na);
		for (int jj = 0; jj < na; jj++)
		{
			fz2_re[jj] = real(z2a[jj]);
			fz2_im[jj] = imag(z2a[jj]);
		}
		const char* pointerz2_re = reinterpret_cast<const char*>(&fz2_re[0]);
		std::size_t bytesz2_re = fz2_re.size() * sizeof(fz2_re[0]);
		outfile_coef.write(pointerz2_re, bytesz2_re);
		const char* pointerz2_im = reinterpret_cast<const char*>(&fz2_im[0]);
		std::size_t bytesz2_im = fz2_im.size() * sizeof(fz2_im[0]);
		outfile_coef.write(pointerz2_im, bytesz2_im);

		dvec fL = La;
		const char* pointerL = reinterpret_cast<const char*>(&fL[0]);
		std::size_t bytesL = fL.size() * sizeof(fL[0]);
		outfile_coef.write(pointerL, bytesL);

		dvec fmu = mua;
		const char* pointermu = reinterpret_cast<const char*>(&fmu[0]);
		std::size_t bytesmu = fmu.size() * sizeof(fmu[0]);
		outfile_coef.write(pointermu, bytesmu);

		// Save the a
		for (int ii = 0; ii < na; ii++)
		{
			dvec fbeta_re(ma);
			dvec fbeta_im(ma);

			for (int jj = 0; jj < ma; jj++)
			{
				fbeta_re[jj] = real(a[ii][jj]);
				fbeta_im[jj] = imag(a[ii][jj]);
			}
			const char* pointerbeta_re = reinterpret_cast<const char*>(&fbeta_re[0]);
			std::size_t bytesbeta_re = fbeta_re.size() * sizeof(fbeta_re[0]);
			outfile_coef.write(pointerbeta_re, bytesbeta_re);
			const char* pointerbeta_im = reinterpret_cast<const char*>(&fbeta_im[0]);
			std::size_t bytesbeta_im = fbeta_im.size() * sizeof(fbeta_im[0]);
			outfile_coef.write(pointerbeta_im, bytesbeta_im);
		}
	}

	// Estimate the time of the program
	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "		ESTIMATED CALCULATION TIME	" << std::endl << std::endl;
	std::cout << "=================================================================" << std::endl << std::endl;
	int Ne = (Nx + Ny + Ntraj) / 100;
	if (Ne < 8)
	{
		Ne = 8;
	}
	auto start0 = std::chrono::high_resolution_clock::now();
	stress_field(xfrom, xto, yfrom, yto, Ne, Ne, sigma_11inf, z1a, z2a, La, mua, ma, na, a);
	auto stop0 = std::chrono::high_resolution_clock::now();
	auto start01 = std::chrono::high_resolution_clock::now();
	principal_stress_trajectories(xfrom, xto, yfrom, yto, xtraj, ytraj, Ne, Ne, sigma_11inf, z1a, z2a, La, mua, ma, na, a);
	auto stop01 = std::chrono::high_resolution_clock::now();
	auto duration0 = std::chrono::duration_cast<std::chrono::microseconds>(stop0 - start0);
	long long ms0 = duration0.count();
	auto duration01 = std::chrono::duration_cast<std::chrono::microseconds>(stop01 - start01);
	long long ms01 = duration01.count();
	long long calcs0 = (int)(2 * Nx * Ny) + (int)(Nw * Nw);
	long long calcs01 = (int)(Ntraj * lvs_traj) + (int)(Ntraj * Nw);
	ms0 = ((ms0 / ((int)Ne * Ne) * calcs0 + ms01 / ((int)Ne * Ne) * calcs01))*0.7;
	long long s0, m0, h0;
	std::tie(ms0, s0, m0, h0) = ms_to_time(ms0);
	std::cout << "Estimated calculation time:  ";
	time_print(ms0, s0, m0, h0);
	std::cout << std::endl << std::endl;
	auto date_time2 = std::chrono::system_clock::now();
	std::time_t start_time2 = std::chrono::system_clock::to_time_t(date_time2);
	char str_time2[26];
	ctime_s(str_time2, sizeof str_time2, &start_time2);
	std::cout << "Plotting started: " << str_time2 << std::endl << std::endl;

	/* --------------------------------------------------------------------
			Calculating the plots
	-----------------------------------------------------------------------*/

	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "		COMPUTING THE PLOTS	" << std::endl << std::endl;
	std::cout << "=================================================================" << std::endl << std::endl;

	// Get the Cartisian stress field
	std::cout << "Initiating the stress field calculation\n";
	auto start1 = std::chrono::high_resolution_clock::now();
	std::tie(x_vec, y_vec, grid_11, grid_22, grid_12) = stress_field(xfrom, xto, yfrom, yto, Nx, Ny, sigma_11inf, z1a, z2a, La, mua, ma, na, a);
	auto stop1 = std::chrono::high_resolution_clock::now();
	auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(stop1 - start1);
	std::cout << "Completed, time taken by function: stress_field = ";
	long long ms1, s1, m1, h1;
	ms1 = duration1.count();
	std::tie(ms1, s1, m1, h1) = ms_to_time(ms1);
	time_print(ms1, s1, m1, h1);
	std::cout << std::endl << std::endl;

	// Get the principal stress field
	std::cout << "Initiating the principal stress field calculation\n";
	auto start2 = std::chrono::high_resolution_clock::now();
	std::tie(grid_1, grid_2, theta_p) = principal_stress_field(Nx, Ny, grid_11, grid_22, grid_12);
	auto stop2 = std::chrono::high_resolution_clock::now();
	auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2 - start2);
	std::cout << "Completed, time taken by function: principal_stress_field = ";
	ms1 = duration2.count();
	std::tie(ms1, s1, m1, h1) = ms_to_time(ms1);
	time_print(ms1, s1, m1, h1);
	std::cout << std::endl << std::endl;

	// Get the stress trajectories
	std::cout << "Initiating the principal stress trajectories calculation\n";
	auto start3 = std::chrono::high_resolution_clock::now();
	std::tie(traj_1, traj_2) = principal_stress_trajectories(xfrom, xto, yfrom, yto, xtraj, ytraj, Ntraj, lvs_traj, sigma_11inf, z1a, z2a, La, mua, ma, na, a);
	auto stop3 = std::chrono::high_resolution_clock::now();
	auto duration3 = std::chrono::duration_cast<std::chrono::microseconds>(stop3 - start3);
	std::cout << "Completed, time taken by function: principal_stress_trajectories = ";
	ms1 = duration3.count();
	std::tie(ms1, s1, m1, h1) = ms_to_time(ms1);
	time_print(ms1, s1, m1, h1);
	std::cout << std::endl << std::endl;

	// Get the displacement field
	std::cout << "Initiating the displacement field calculation\n";
	auto start4 = std::chrono::high_resolution_clock::now();
	std::tie(x_vecw, y_vecw, grid_w) = w_field(xfrom, xto, yfrom, yto, Nw, kappa, G, sigma_11inf, z1a, z2a, La, mua, ma, na, a);
	auto stop4 = std::chrono::high_resolution_clock::now();
	auto duration4 = std::chrono::duration_cast<std::chrono::microseconds>(stop4 - start4);
	std::cout << "Completed, time taken by function: w_field = ";
	ms1 = duration4.count();
	std::tie(ms1, s1, m1, h1) = ms_to_time(ms1);
	time_print(ms1, s1, m1, h1);
	std::cout << std::endl << std::endl;


	// Get the displacement trajectories
	std::cout << "Initiating the displacement trajectories calculation\n";
	auto start5 = std::chrono::high_resolution_clock::now();
	traj_w = w_trajectories(xfrom, xto, yfrom, yto, Ntraj, Nw, kappa, G, sigma_11inf, z1a, z2a, La, mua, ma, na, a);
	auto stop5 = std::chrono::high_resolution_clock::now();
	auto duration5 = std::chrono::duration_cast<std::chrono::microseconds>(stop5 - start5);
	std::cout << "Completed, time taken by function: w_trajectories = ";
	ms1 = duration5.count();
	std::tie(ms1, s1, m1, h1) = ms_to_time(ms1);
	time_print(ms1, s1, m1, h1);
	std::cout << std::endl << std::endl;

	// Displaying the computation time
	std::cout << "=================================================================" << std::endl << std::endl;
	long long mstime = duration1.count() + duration2.count() + duration3.count() + duration4.count() + duration5.count();
	long long stime, mtime, htime;
	std::tie(mstime, stime, mtime, htime) = ms_to_time(mstime);
	std::cout << "Total calculation time: ";
	time_print(mstime, stime, mtime, htime);
	std::cout << std::endl << std::endl;

	/* --------------------------------------------------------------------
			Saving the data as binary files
	-----------------------------------------------------------------------*/

	// Create/open the two output files
	std::ofstream outfile("data.bin", std::ios::out | std::ios::binary);
	std::ofstream outfiledim("dim_data.bin", std::ios::out | std::ios::binary);

	// Save the x and y vectors
	dvec fx = x_vec;
	const char* pointerx = reinterpret_cast<const char*>(&fx[0]);
	std::size_t bytesx = fx.size() * sizeof(fx[0]);
	outfile.write(pointerx, bytesx);

	dvec fy = y_vec;
	const char* pointery = reinterpret_cast<const char*>(&fy[0]);
	std::size_t bytesy = fy.size() * sizeof(fy[0]);
	outfile.write(pointery, bytesy);

	dvec fxw = x_vecw;
	const char* pointerxw = reinterpret_cast<const char*>(&fxw[0]);
	std::size_t bytesxw = fxw.size() * sizeof(fxw[0]);
	outfile.write(pointerxw, bytesxw);

	dvec fyw = y_vecw;
	const char* pointeryw = reinterpret_cast<const char*>(&fyw[0]);
	std::size_t bytesyw = fyw.size() * sizeof(fyw[0]);
	outfile.write(pointeryw, bytesyw);

	// Save the grids
	for (size_t ii = 0; ii < grid_11.size(); ii++)
	{
		dvec fg11 = grid_11[ii];
		const char* pointerg11 = reinterpret_cast<const char*>(&fg11[0]);
		std::size_t bytesg11 = fg11.size() * sizeof(fg11[0]);
		outfile.write(pointerg11, bytesg11);
	}
	for (size_t ii = 0; ii < grid_22.size(); ii++)
	{
		dvec fg22 = grid_22[ii];
		const char* pointerg22 = reinterpret_cast<const char*>(&fg22[0]);
		std::size_t bytesg22 = fg22.size() * sizeof(fg22[0]);
		outfile.write(pointerg22, bytesg22);
	}
	for (size_t ii = 0; ii < grid_12.size(); ii++)
	{
		dvec fg12 = grid_12[ii];
		const char* pointerg12 = reinterpret_cast<const char*>(&fg12[0]);
		std::size_t bytesg12 = fg12.size() * sizeof(fg12[0]);
		outfile.write(pointerg12, bytesg12);
	}
	for (size_t ii = 0; ii < grid_1.size(); ii++)
	{
		dvec fg1 = grid_1[ii];
		const char* pointerg1 = reinterpret_cast<const char*>(&fg1[0]);
		std::size_t bytesg1 = fg1.size() * sizeof(fg1[0]);
		outfile.write(pointerg1, bytesg1);
	}
	for (size_t ii = 0; ii < grid_2.size(); ii++)
	{
		dvec fg2 = grid_2[ii];
		const char* pointerg2 = reinterpret_cast<const char*>(&fg2[0]);
		std::size_t bytesg2 = fg2.size() * sizeof(fg2[0]);
		outfile.write(pointerg2, bytesg2);
	}
	for (size_t ii = 0; ii < theta_p.size(); ii++)
	{
		dvec fgtp = theta_p[ii];
		const char* pointergtp = reinterpret_cast<const char*>(&fgtp[0]);
		std::size_t bytesgtp = fgtp.size() * sizeof(fgtp[0]);
		outfile.write(pointergtp, bytesgtp);
	}
	for (size_t ii = 0; ii < traj_1.size(); ii++)
	{
		dvec fgt1_re(Ntraj);
		dvec fgt1_im(Ntraj);

		for (size_t jj = 0; jj < traj_1[0].size(); jj++)
		{
			fgt1_re[jj] = real(traj_1[ii][jj]);
			fgt1_im[jj] = imag(traj_1[ii][jj]);
		}
		const char* pointergt1_re = reinterpret_cast<const char*>(&fgt1_re[0]);
		std::size_t bytesgt1_re = fgt1_re.size() * sizeof(fgt1_re[0]);
		outfile.write(pointergt1_re, bytesgt1_re);
		const char* pointergt1_im = reinterpret_cast<const char*>(&fgt1_im[0]);
		std::size_t bytesgt1_im = fgt1_im.size() * sizeof(fgt1_im[0]);
		outfile.write(pointergt1_im, bytesgt1_im);
	}
	for (size_t ii = 0; ii < traj_2.size(); ii++)
	{
		dvec fgt2_re(Ntraj);
		dvec fgt2_im(Ntraj);
		for (size_t jj = 0; jj < traj_2[0].size(); jj++)
		{
			fgt2_re[jj] = real(traj_2[ii][jj]);
			fgt2_im[jj] = imag(traj_2[ii][jj]);
		}
		const char* pointergt2_re = reinterpret_cast<const char*>(&fgt2_re[0]);
		std::size_t bytesgt2_re = fgt2_re.size() * sizeof(fgt2_re[0]);
		outfile.write(pointergt2_re, bytesgt2_re);
		const char* pointergt2_im = reinterpret_cast<const char*>(&fgt2_im[0]);
		std::size_t bytesgt2_im = fgt2_im.size() * sizeof(fgt2_im[0]);
		outfile.write(pointergt2_im, bytesgt2_im);
	}
	for (size_t ii = 0; ii < grid_w.size(); ii++)
	{
		dvec fgw_re(Nw);
		dvec fgw_im(Nw);
		for (size_t jj = 0; jj < grid_w[0].size(); jj++)
		{
			fgw_re[jj] = real(grid_w[ii][jj]);
			fgw_im[jj] = imag(grid_w[ii][jj]);
		}
		const char* pointergw_re = reinterpret_cast<const char*>(&fgw_re[0]);
		std::size_t bytesgw_re = fgw_re.size() * sizeof(fgw_re[0]);
		outfile.write(pointergw_re, bytesgw_re);
		const char* pointergw_im = reinterpret_cast<const char*>(&fgw_im[0]);
		std::size_t bytesgw_im = fgw_im.size() * sizeof(fgw_im[0]);
		outfile.write(pointergw_im, bytesgw_im);
	}
	for (size_t ii = 0; ii < traj_w.size(); ii++)
	{
		dvec fgtw_re(Ntraj);
		dvec fgtw_im(Ntraj);
		for (size_t jj = 0; jj < traj_w[0].size(); jj++)
		{
			fgtw_re[jj] = real(traj_w[ii][jj]);
			fgtw_im[jj] = imag(traj_w[ii][jj]);
		}
		const char* pointergtw_re = reinterpret_cast<const char*>(&fgtw_re[0]);
		std::size_t bytesgtw_re = fgtw_re.size() * sizeof(fgtw_re[0]);
		outfile.write(pointergtw_re, bytesgtw_re);
		const char* pointergtw_im = reinterpret_cast<const char*>(&fgtw_im[0]);
		std::size_t bytesgtw_im = fgtw_im.size() * sizeof(fgtw_im[0]);
		outfile.write(pointergtw_im, bytesgtw_im);
	}
	if (na > 0)
	{
		dvec fT_res(numa * na);
		dvec fT_ims(numa * na);
		for (int jj = 0; jj < numa * na; jj++)
		{
			fT_res[jj] = real(T_check[jj]);
			fT_ims[jj] = imag(T_check[jj]);
		}
		const char* pointerT_res = reinterpret_cast<const char*>(&fT_res[0]);
		std::size_t bytesT_res = fT_res.size() * sizeof(fT_res[0]);
		outfile.write(pointerT_res, bytesT_res);
		const char* pointerT_ims = reinterpret_cast<const char*>(&fT_ims[0]);
		std::size_t bytesT_ims = fT_ims.size() * sizeof(fT_ims[0]);
		outfile.write(pointerT_ims, bytesT_ims);
	}
	

	// Save the plot properties
	dvec dim = { 1.0 * Nx, 1.0 * Ny, 1.0 * Nw, 1.0 * Ntraj, 1.0 * lvs_traj, 1.0 * na, error_med_a_re, error_mean_a_re, error_max_a_re, error_med_a_im, error_mean_a_im, error_max_a_im, error_med_int, error_mean_int, error_max_int};
	const char* pointerdim = reinterpret_cast<const char*>(&dim[0]);
	std::size_t bytesdim = dim.size() * sizeof(dim[0]);
	outfiledim.write(pointerdim, bytesdim);

	// saving the coordinates of the cracks
	if (na > 0)
	{
		dvec fz1_res(na);
		dvec fz1_ims(na);
		for (int jj = 0; jj < na; jj++)
		{
			fz1_res[jj] = real(z1a[jj]);
			fz1_ims[jj] = imag(z1a[jj]);
		}
		const char* pointerz1_res = reinterpret_cast<const char*>(&fz1_res[0]);
		std::size_t bytesz1_res = fz1_res.size() * sizeof(fz1_res[0]);
		outfiledim.write(pointerz1_res, bytesz1_res);
		const char* pointerz1_ims = reinterpret_cast<const char*>(&fz1_ims[0]);
		std::size_t bytesz1_ims = fz1_ims.size() * sizeof(fz1_ims[0]);
		outfiledim.write(pointerz1_ims, bytesz1_ims);

		dvec fz2_res(na);
		dvec fz2_ims(na);
		for (int jj = 0; jj < na; jj++)
		{
			fz2_res[jj] = real(z2a[jj]);
			fz2_ims[jj] = imag(z2a[jj]);
		}
		const char* pointerz2_re = reinterpret_cast<const char*>(&fz2_res[0]);
		std::size_t bytesz2_re = fz2_res.size() * sizeof(fz2_res[0]);
		outfiledim.write(pointerz2_re, bytesz2_re);
		const char* pointerz2_im = reinterpret_cast<const char*>(&fz2_ims[0]);
		std::size_t bytesz2_im = fz2_ims.size() * sizeof(fz2_ims[0]);
		outfiledim.write(pointerz2_im, bytesz2_im);

	}

	// Close the output files
	outfile.close();
	outfiledim.close();

	// Get the date and execution time
	auto end = std::chrono::high_resolution_clock::now();
	auto date_time3 = std::chrono::system_clock::now();
	auto elapsed_seconds = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
	long long mseconds = elapsed_seconds.count();
	long long seconds, hours, minutes;
	std::tie(mseconds, seconds, minutes, hours) = ms_to_time(mseconds);
	std::time_t end_time = std::chrono::system_clock::to_time_t(date_time3);
	char str_time[26];
	ctime_s(str_time, sizeof str_time, &end_time);

	std::cout << "Program finnished after ";
	time_print(mseconds, seconds, minutes, hours);
	std::cout << std::endl;
	std::cout << "Output data saved to binary files: data.bin and dim_data.bin" << std::endl;

	return 0;
}