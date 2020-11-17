#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <numeric>
#include <tuple>
#include <fstream>
#include <ctime>
#include <chrono>
#include <omp.h>
#include <stdio.h>
#include <Eigen/Dense>

// Defining special types
typedef std::complex<double> dcomp;
typedef std::vector<double> dvec;
typedef std::vector< std::vector<double> > ddvec;
typedef std::vector< std::complex<double> > dcvec;
typedef std::vector< std::vector<std::complex<double> > > ddcvec;

/* --------------------------------------------------------------------

		FUNCTIONS

-----------------------------------------------------------------------*/

inline const double pi()
{
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
	z0 = 0.5*(z1 + z2);
	Z = exp(dcomp(0, -1) * mu) * 2.0 * (z - z0) / L;
	chi = Z + sqrt(Z - 1.0)*sqrt(Z + 1.0);

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
	z0 = 0.5*(z1 + z2);
	Z = 0.5*(chi + 1.0 / chi);
	z = 0.5*L*Z*exp(dcomp(0, 1)*mu) + z0;

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
		TAU UNIFORM STRESS
-----------------------------------------------------------------------*/
inline std::tuple<dcomp, dcomp>  tau_uni(double sigma_11inf)
{
	// Defining the variables
	dcomp tau_11, tau_12, phi, psi;

	// Get the phi and psi
	phi = -0.5*sigma_11inf;
	psi = -0.5*sigma_11inf;

	// calculating the tau
	tau_11 = -phi - psi;
	tau_12 = -phi - phi;

	return { tau_11, tau_12 };
}

/* --------------------------------------------------------------------
		TAU CRACK
-----------------------------------------------------------------------*/
inline std::tuple<dcomp, dcomp> tau_crack(dcomp z, dcomp z1, dcomp z2, double L, double mu, int m, dcvec beta)
{
	// Defining the variables
	dcomp chi, chi_bar, Z, chi_pow;
	dcomp dphi, dphi_bar, ddphi, dpsi;
	dcomp tau_11, tau_12, S1, L_frac;
	double n;

	// Getting the chi - and Z - coordinates
	chi = chi_from_z(z, z1, z2, L, mu);
	chi_bar = conj(chi);
	Z = exp(dcomp(0, -1)*mu) * 2.0 * (z - 0.5*(z1 + z2)) / L;

	// Calculating the series
	dphi = 0;
	dphi_bar = 0;
	ddphi = 0;
	dpsi = 0;
	n = 0;
	chi_pow = chi * chi - 1.0;
	for (int ii = 0; ii < m; ii++)
	{
		dcomp beta_n;
		n += 1;
		beta_n = beta[ii] * n;
		dphi += conj(beta_n) * pow(chi, (1.0 - n)) / chi_pow;
		dphi_bar += beta_n * pow(chi_bar, (1.0 - n)) / (chi_bar * chi_bar - 1.0);
		ddphi -= conj(beta_n) * (pow(chi, (2.0 - n)) / (chi_pow*chi_pow*chi_pow))*((n + 1.0)*chi*chi - n + 1.0);
		dpsi -= beta_n * (pow(chi, (1.0 - n))) / chi_pow;
	}

	// Multiplying the constants
	L_frac = (4.0 / L) * exp(dcomp(0, -1)*mu);
	dphi *= L_frac;
	dphi_bar *= conj(L_frac);
	ddphi *= (16.0 / (L*L))*exp(dcomp(0, -2)*mu);
	dpsi *= L_frac;

	// Calcualting tau
	tau_11 = -0.5*L*(Z - conj(Z))*ddphi - exp(dcomp(0, -1)*mu)*(dphi + dpsi);
	tau_12 = -exp(dcomp(0, 1)*mu)*dphi - exp(dcomp(0, -1)*mu)*dphi_bar;

	return { tau_11, tau_12 };
}

/* --------------------------------------------------------------------
		TAU TOTAL
-----------------------------------------------------------------------*/
inline std::tuple<dcomp, dcomp>  tau_total(dcomp z, double sigma_11inf, dcvec z1, dcvec z2, dvec L, dvec mu, int m, int nc, ddcvec beta, int m_not)
{
	// Defining the variables
	dcomp tau_11, tau_12;

	std::tie(tau_11, tau_12) = tau_uni(sigma_11inf);
	if (m > 0)
	{
		for (int ii = 0; ii < nc; ii++)
		{
			if (ii != m_not)
			{
				dcomp tau_11c, tau_12c;
				std::tie(tau_11c, tau_12c) = tau_crack(z, z1[ii], z2[ii], L[ii], mu[ii], m, beta[ii]);
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
inline dcomp T_total(dcomp z, double sigma_11inf, dcvec z1, dcvec z2, dvec L, dvec mu, int m, int nc, ddcvec beta, int m_not, int m_is)
{
	// Defining the variables
	dcomp tau_11, tau_12, T;

	std::tie(tau_11, tau_12) = tau_total(z, sigma_11inf, z1, z2, L, mu, m, nc, beta, m_not);

	T = -0.5 * dcomp(0, 1)*(tau_11 * exp(dcomp(0, 2)*mu[m_is]) - tau_12);


	return { T };
}

/* --------------------------------------------------------------------
		SOLVE CRACK BETA
-----------------------------------------------------------------------*/
inline dcvec AE_crack_solver(Eigen::VectorXd T_s, Eigen::VectorXd T_n, dvec term, Eigen::MatrixXd A, int N, int m)
{
	// Defining the variables
	dcvec beta(m);

	// Solving the linear system
	Eigen::VectorXd b1 = A.colPivHouseholderQr().solve(T_s);
	Eigen::VectorXd b2 = A.colPivHouseholderQr().solve(T_n);

	// Assigning to beta
	#pragma omp parallel for default(none) shared(beta)
	for (int ii = 0; ii < m; ii++)
	{
		beta[ii] = dcomp(b2[ii], b1[ii]);
	}


	return { beta };
}

/* --------------------------------------------------------------------
		ITERATOR
-----------------------------------------------------------------------*/
inline ddcvec iterator(double cond, int ITR, int N, dvec p, double sigma_11inf, dcvec z1, dcvec z2, dvec L, dvec mu, int m, int nc)
{
	// Defining the variables
	dcomp tau_11, tau_12, T;
	Eigen::MatrixXd A(N, m);
	dvec theta(N);
	ddvec term(nc, dvec(N));
	ddcvec z(nc, dcvec(N));
	ddcvec beta(nc, dcvec(m));
	double theta_0, deltheta;

	// Print the conditions
	std::cout << std::scientific;
	std::cout << "Solver for " << nc << " cracks with " << m << " coefficients at " << N << " integration points." << std::endl;
	std::cout << "Iterations break after: error < " << cond << " or iterations > " << ITR << "." << std::endl << std::endl;
	std::cout << "	Error:" << "		Iteration:" << std::endl;


	// Assigning variables
	theta_0 = 0.1;
	deltheta = 2 * pi() / N;

	// Calculating the A, a theta and term matrices
#pragma omp parallel for default(none) shared(theta)
	for (int ii = 0; ii < N; ii++)
	{
		theta[ii] = theta_0 + (ii + 1) * deltheta;
	}

#pragma omp parallel for default(none) shared(A, z, term)
	for (int ii = 0; ii < N; ii++)
	{
		for (int mm = 0; mm < m; mm++)
		{
			A(ii, mm) = (mm + 1) * sin((mm + 1)*theta[ii]);
		}
		for (int jj = 0; jj < nc; jj++)
		{
			dcomp chi;
			chi = exp(dcomp(0, 1)*theta[ii]);
			z[jj][ii] = z_from_chi(chi, z1[jj], z2[jj], L[jj], mu[jj]);
			term[jj][ii] = -0.5*L[jj] * sin(theta[ii]);
		}
	}

	// Sovle the betas
	double error = 1;
	ddcvec error_b(nc, dcvec(m, 0));
	int NIT = 0;
	while (error > cond && NIT < ITR)
	{
		for (int ii = 0; ii < nc; ii++)
		{
			Eigen::VectorXd T_s(N);
			Eigen::VectorXd T_n(N);
#pragma omp parallel for default(none) shared(T_s, T_n)
			for (int jj = 0; jj < N; jj++)
			{
				dcomp T;
				T = T_total(z[ii][jj], sigma_11inf, z1, z2, L, mu, m, nc, beta, ii, ii);
				T_s(jj) = -real(T)*term[ii][jj];
				T_n(jj) = (p[ii] + imag(T))*term[ii][jj];
			}
			beta[ii] = AE_crack_solver(T_s, T_n, term[ii], A, N, m);

		}

		// Calcualte the error
		dvec delbeta(nc);
#pragma omp parallel for default(none) shared(delbeta)
		for (int ii = 0; ii < nc; ii++)
		{
			dvec delb_temp(m);
			for (int jj = 0; jj < m; jj++)
			{
				delb_temp[jj] = abs(beta[ii][jj] - error_b[ii][jj]);
			}
			delbeta[ii] = *max_element(delb_temp.begin(), delb_temp.end());
		}
		error = *max_element(delbeta.begin(), delbeta.end());

		NIT += 1;
		error_b = beta;

		std::cout << std::scientific;
		std::cout << "	" << error << "	" << NIT << std::endl;
	}



	return { beta };
}

/* --------------------------------------------------------------------
		STRESS TOTAL
-----------------------------------------------------------------------*/
inline std::tuple<double, double, double> sigma_total(dcomp z, double sigma_11inf, dcvec z1, dcvec z2, dvec L, dvec mu, int m, int nc, ddcvec beta, int m_not)
{
	// Defining the variables
	dcomp S1, S2, tau_11, tau_12;
	double sigma_11, sigma_22, sigma_12;

	// Calculating the tau
	std::tie(tau_11, tau_12) = tau_total(z, sigma_11inf, z1, z2, L, mu, m, nc, beta, m_not);

	// Calculate the sigmas
	S1 = .5*(tau_11 + tau_12);
	S2 = .5*(-tau_11 + tau_12);
	sigma_11 = real(S1);
	sigma_22 = real(S2);
	sigma_12 = -imag(S1);

	return { sigma_11, sigma_22, sigma_12 };
}

/* --------------------------------------------------------------------
		STRESS FIELD
-----------------------------------------------------------------------*/
std::tuple<dvec, dvec, ddvec, ddvec, ddvec> stress_field(double xfrom, double xto, double yfrom, double yto, int Nx, int Ny, double sigma_11inf, dcvec z1, dcvec z2, dvec L, dvec mu, int m, int nc, ddcvec beta, int m_not)
{
	// Defining the variables
	ddvec grid_11(Nx, dvec(Ny));
	ddvec grid_22(Nx, dvec(Ny));
	ddvec grid_12(Nx, dvec(Ny));
	double dx;
	double dy;
	dvec x_vec(Nx);
	dvec y_vec(Ny);

	// Retriving the terms from the sigma function for z
	dx = (xto - xfrom) / (Nx - 1);
	dy = (yto - yfrom) / (Ny - 1);
#pragma omp parallel for default(none) shared(grid_11, grid_22, grid_12, x_vec, y_vec)
	for (int ii = 0; ii < Nx; ii++)
	{
		for (int jj = Ny; jj--;)
		{
			x_vec[ii] = xfrom + ii * dx;
			y_vec[jj] = yfrom + jj * dy;
			std::tie(grid_11[ii][jj], grid_22[ii][jj], grid_12[ii][jj]) = sigma_total(dcomp(x_vec[ii], y_vec[jj]), sigma_11inf, z1, z2, L, mu, m, nc, beta, m_not);
		}
	}
	return { x_vec, y_vec, grid_11, grid_22, grid_12 };
}

/* --------------------------------------------------------------------
		PRINCIPAL STRESSES
-----------------------------------------------------------------------*/
std::tuple<double, double, double> principal_sigma(dcomp z, double sigma_11inf, dcvec z1, dcvec z2, dvec L, dvec mu, int m, int nc, ddcvec beta, int m_not)
{
	// Defining the variables
	double sigma_1;
	double sigma_2;
	double theta_p;
	double frac1, frac2, sqrt1;
	dcomp S1, S2, tau_11, tau_12;
	double sigma_11, sigma_22, sigma_12;

	// Calculating the tau
	std::tie(tau_11, tau_12) = tau_total(z, sigma_11inf, z1, z2, L, mu, m, nc, beta, m_not);

	// Calculate the sigmas
	S1 = .5*(tau_11 + tau_12);
	S2 = .5*(-tau_11 + tau_12);
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
	theta_p = -0.5*imag(log(tau_11));

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
std::tuple<ddvec, ddvec, ddvec> principal_stress_field(double xfrom, double xto, double yfrom, double yto, int Nx, int Ny, double sigma_11inf, dcvec z1, dcvec z2, dvec L, dvec mu, int m, int nc, ddcvec beta, int m_not)
{
	// Defining the variables
	ddvec grid_1(Nx, dvec(Ny));
	ddvec grid_2(Nx, dvec(Ny));
	ddvec grid_tp(Nx, dvec(Ny));
	double dx;
	double dy;
	dvec x_vec(Nx);
	dvec y_vec(Ny);

	// Retriving the terms from the sigma function for z
	dx = (xto - xfrom) / (Nx - 1);
	dy = (yto - yfrom) / (Ny - 1);
	#pragma omp parallel for default(none) shared(grid_1, grid_2, grid_tp, x_vec, y_vec)
	for (int ii = 0; ii < Nx; ii++)
	{
		for (int jj = Ny; jj--;)
		{
			x_vec[ii] = xfrom + ii * dx;
			y_vec[jj] = yfrom + jj * dy;
			std::tie(grid_1[ii][jj], grid_2[ii][jj], grid_tp[ii][jj]) = principal_sigma(dcomp(x_vec[ii], y_vec[jj]), sigma_11inf, z1, z2, L, mu, m, nc, beta, m_not);
		}
	}
	return { grid_1, grid_2, grid_tp };
}

/* --------------------------------------------------------------------
		PRINCIPAL STRESS TRAJECTORIES
-----------------------------------------------------------------------*/
std::tuple<ddcvec, ddcvec> principal_stress_trajectories(double xfrom, double xto, double yfrom, double yto, dcvec xtraj, dcvec ytraj, int Ntraj, int lvs_traj, double sigma_11inf, dcvec z1, dcvec z2, dvec L, dvec mu, int m, int nc, ddcvec beta, int m_not)
{
	// Defining the variables
	ddcvec traj_1(lvs_traj * 2, dcvec(Ntraj));
	ddcvec traj_2(lvs_traj * 2, dcvec(Ntraj));
	double dx, dy, dx_lvsre, dx_lvsim, dy_lvsre, dy_lvsim, pi_val;
	double cond;
	int NIT;

	// Getting the starting points
	dx = (xto - xfrom) / (Ntraj*0.8);
	dy = (yto - yfrom) / (Ntraj*0.8);
	dx_lvsre = (real(xtraj[1]) - real(xtraj[0])) / (lvs_traj - 1);
	dx_lvsim = (imag(xtraj[1]) - imag(xtraj[0])) / (lvs_traj - 1);
	dy_lvsre = (real(ytraj[1]) - real(ytraj[0])) / (lvs_traj - 1);
	dy_lvsim = (imag(ytraj[1]) - imag(ytraj[0])) / (lvs_traj - 1);
	pi_val = 0.5*pi();
	cond = 1e-6;
	NIT = 10;

	// SIGMA 1
	#pragma omp parallel for default(none) shared(traj_1)
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
			std::tie(sigma_1, std::ignore, theta_p) = principal_sigma(z, sigma_11inf, z1, z2, L, mu, m, nc, beta, m_not);
			sigma = abs(sigma_1)*exp(dcomp(0, 1)*theta_p);
			z11 = z + sigma / abs(sigma)*dx1;
			eta = angel_change(z11, z, z_old);
			if (eta > pi_val && jj > 1) {
				dx1 = -dx1;
			}
			zt = z + sigma / abs(sigma)*dx1;

			ee = 1;
			for (int rr = NIT; rr--;)
			{
				z_oldt = zt;
				std::tie(sigma_1, std::ignore, theta_p) = principal_sigma(zt, sigma_11inf, z1, z2, L, mu, m, nc, beta, m_not);
				sigmat = abs(sigma_1)*exp(dcomp(0, 1)*theta_p);
				z11 = z + (sigma + sigmat) / abs(sigma + sigmat)*dx1;
				eta = angel_change(z11, z, z_old);
				if (eta > pi_val && jj > 1) {
					dx1 = -dx1;
				}
				zt = z + (sigma + sigmat) / abs(sigma + sigmat)*dx1;
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
			std::tie(sigma_1, std::ignore, theta_p) = principal_sigma(z, sigma_11inf, z1, z2, L, mu, m, nc, beta, m_not);
			sigma = abs(sigma_1)*exp(dcomp(0, 1)*theta_p);
			z11 = z + sigma / abs(sigma)*dx1;
			eta = angel_change(z11, z, z_old);
			if (eta > pi_val && jj > 1) {
				dx1 = -dx1;
			}
			zt = z + sigma / abs(sigma)*dx1;

			ee = 1;
			for (int rr = NIT; rr--;)
			{
				z_oldt = zt;
				std::tie(sigma_1, std::ignore, theta_p) = principal_sigma(zt, sigma_11inf, z1, z2, L, mu, m, nc, beta, m_not);
				sigmat = abs(sigma_1)*exp(dcomp(0, 1)*theta_p);
				z11 = z + (sigma + sigmat) / abs(sigma + sigmat)*dx1;
				eta = angel_change(z11, z, z_old);
				if (eta > pi_val && jj > 1) {
					dx1 = -dx1;
				}
				zt = z + (sigma + sigmat) / abs(sigma + sigmat)*dx1;
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
	#pragma omp parallel for default(none) shared(traj_2)
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
			std::tie(std::ignore, sigma_2, theta_p) = principal_sigma(z, sigma_11inf, z1, z2, L, mu, m, nc, beta, m_not);
			sigma = abs(sigma_2)*exp(dcomp(0, 1)*(theta_p + pi_val));
			z11 = z + sigma / abs(sigma)*dy1;
			eta = angel_change(z11, z, z_old);
			if (eta > pi_val && jj > 1) {
				dy1 = -dy1;
			}
			zt = z + sigma / abs(sigma)*dy1;

			ee = 1;
			for (int rr = NIT; rr--;)
			{
				z_oldt = zt;
				std::tie(std::ignore, sigma_2, theta_p) = principal_sigma(zt, sigma_11inf, z1, z2, L, mu, m, nc, beta, m_not);
				sigmat = abs(sigma_2)*exp(dcomp(0, 1)*(theta_p + pi_val));
				z11 = z + (sigma + sigmat) / abs(sigma + sigmat)*dy1;
				eta = angel_change(z11, z, z_old);
				if (eta > pi_val && jj > 1) {
					dy1 = -dy1;
				}
				zt = z + (sigma + sigmat) / abs(sigma + sigmat)*dy1;
				ee = std::norm(z_oldt - zt);
				if (ee < cond)
				{
					break;
				}
			}

			traj_2[ii][jj] = zt;
			z_old = traj_2[ii][jj - 1];
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
			std::tie(std::ignore, sigma_2, theta_p) = principal_sigma(z, sigma_11inf, z1, z2, L, mu, m, nc, beta, m_not);
			sigma = abs(sigma_2)*exp(dcomp(0, 1)*(theta_p + pi_val));
			z11 = z + sigma / abs(sigma)*dy1;
			eta = angel_change(z11, z, z_old);
			if (eta > pi_val && jj > 1) {
				dy1 = -dy1;
			}
			zt = z + sigma / abs(sigma)*dy1;

			ee = 1;
			for (int rr = NIT; rr--;)
			{
				z_oldt = zt;
				std::tie(std::ignore, sigma_2, theta_p) = principal_sigma(zt, sigma_11inf, z1, z2, L, mu, m, nc, beta, m_not);
				sigmat = abs(sigma_2)*exp(dcomp(0, 1)*(theta_p + pi_val));
				z11 = z + (sigma + sigmat) / abs(sigma)*dy1;
				eta = angel_change(z11, z, z_old);
				if (eta > pi_val && jj > 1) {
					dy1 = -dy1;
				}
				zt = z + (sigma + sigmat) / abs(sigma + sigmat)*dy1;
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
	phi_bar = -0.5*sigma_11inf*conj(z);
	dphi = -0.5*sigma_11inf;
	psi = -0.5*sigma_11inf*z;

	// Calculating w
	w = 1 / (4 * G)*((z - conj(z))*dphi + kappa * phi_bar + psi);

	return { w };
}

/* --------------------------------------------------------------------
		W CRACK
-----------------------------------------------------------------------*/
inline dcomp w_crack(dcomp z, double kappa, double G, dcomp z1, dcomp z2, double L, double mu, int m, dcvec beta)
{
	// Defining the variables
	dcomp chi, chi_bar, Z, chi_pow;
	dcomp dphi, phi_bar, psi;
	dcomp w, L_frac;
	double n;

	// Getting the chi - and Z - coordinates
	chi = chi_from_z(z, z1, z2, L, mu);
	chi_bar = conj(chi);
	Z = exp(dcomp(0, -1)*mu) * 2.0 * (z - 0.5*(z1 + z2)) / L;

	// Calculating the series
	phi_bar = 0;
	dphi = 0;
	psi = 0;
	n = 0;
	chi_pow = chi * chi - 1.0;
	for (int ii = 0; ii < m; ii++)
	{
		dcomp beta_n;
		n += 1;
		beta_n = beta[ii] * n;
		dphi += conj(beta_n) * pow(chi, (1.0 - n)) / chi_pow;
		phi_bar -= beta[ii] * pow(chi_bar, -n);
		psi += beta[ii] * pow(chi, -n);
	}

	// Multiplying the constants
	L_frac = (4.0 / L) * exp(dcomp(0, -1)*mu);
	dphi *= L_frac;

	// Calcualting w
	w = 1 / (4 * G)*(0.5*L*(Z - conj(Z))*dphi + exp(dcomp(0, -1)*mu)*kappa*phi_bar + exp(dcomp(0, -1)*mu)*psi);

	return { w };
}

/* --------------------------------------------------------------------
		W TOTAL
-----------------------------------------------------------------------*/
inline dcomp  w_total(dcomp z, double kappa, double G, double sigma_11inf, dcvec z1, dcvec z2, dvec L, dvec mu, int m, int nc, ddcvec beta)
{
	// Defining the variables
	dcomp w, wg;

	w = w_uni(z, kappa, G, sigma_11inf);
	if (m > 0)
	{
		for (int ii = 0; ii < nc; ii++)
		{
			dcomp wc;
			wc = w_crack(z, kappa, G, z1[ii], z2[ii], L[ii], mu[ii], m, beta[ii]);
			w += wc;
		}
	}

	return { w };
}

/* --------------------------------------------------------------------
		DISPLACEMENT FIELD
-----------------------------------------------------------------------*/
std::tuple<dvec, dvec, ddcvec> w_field(double xfrom, double xto, double yfrom, double yto, int Nw, double kappa, double G, double sigma_11inf, dcvec z1, dcvec z2, dvec L, dvec mu, int m, int nc, ddcvec beta)

{
	// Defining the variables
	ddcvec grid_w(Nw, dcvec(Nw));
	dvec x_vecw(Nw), y_vecw(Nw);
	double dx;
	double dy;

	// Retriving the terms from the sigma function for z
	dx = (xto - xfrom) / (Nw - 1.0);
	dy = (yto - yfrom) / (Nw - 1.0);
	#pragma omp parallel for default(none) shared(grid_w, x_vecw, y_vecw)
	for (int ii = 0; ii < Nw; ii++)
	{
		for (int jj = 0; jj < Nw; jj++)
		{

			x_vecw[ii] = xfrom + ii * dx;
			y_vecw[jj] = yfrom + jj * dy;
			grid_w[ii][jj] = w_total(dcomp(x_vecw[ii], y_vecw[jj]), kappa, G, sigma_11inf, z1, z2, L, mu, m, nc, beta);
		}
	}
	return { x_vecw, y_vecw, grid_w };
}

/* --------------------------------------------------------------------
		DISPLACEMENT TRAJECTORIES
-----------------------------------------------------------------------*/
ddcvec w_trajectories(double xfrom, double xto, double yfrom, double yto, int Ntraj, int Nw, double kappa, double G, double sigma_11inf, dcvec z1, dcvec z2, dvec L, dvec mu, int m, int nc, ddcvec beta)
{
	// Defining the variables
	ddcvec traj_w(Nw * 2, dcvec(Ntraj));
	double dx, dy_lvs;
	double cond, pi_val;
	int NIT;

	// Getting the starting points
	dx = (xto - xfrom) / (Ntraj);
	dy_lvs = (yto - yfrom) / (Nw - 1.0);
	pi_val = 0.5*pi();
	cond = 1e-6;
	NIT = 10;

	// w trajectories
	#pragma omp parallel for default(none) shared(traj_w)
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
			w = w_total(z, kappa, G, sigma_11inf, z1, z2, L, mu, m, nc, beta);
			zt = z + conj(w) / abs(w)*dx1;

			ee = 1;
			for (int rr = NIT; rr--;)
			{
				z_oldt = zt;
				w1 = w_total(zt, kappa, G, sigma_11inf, z1, z2, L, mu, m, nc, beta);
				zt = z + conj(w + w1) / abs(w + w1)*dx1;
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
			w = w_total(z, kappa, G, sigma_11inf, z1, z2, L, mu, m, nc, beta);
			zt = z + conj(w) / abs(w)*dx1;

			ee = 1;
			for (int rr = NIT; rr--;)
			{
				z_oldt = zt;
				w1 = w_total(zt, kappa, G, sigma_11inf, z1, z2, L, mu, m, nc, beta);
				zt = z + conj(w + w1) / abs(w + w1)*dx1;
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
	std::cout << "		AEM CRACK SOLVER	" << std::endl << std::endl;
	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "The documentation for this program can be found on:\nhttps://github.com/eriktoller/AE_LE_crack \n";
	std::cout << "Written by: Erik Toller, erik.toller@geo.uu.se.\n\n";

	auto date_time1 = std::chrono::system_clock::now();
	std::time_t start_time = std::chrono::system_clock::to_time_t(date_time1);
	char str_time1[26];
	ctime_s(str_time1, sizeof str_time1, &start_time);
	std::cout << "Program started: " << str_time1 << std::endl;
	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "		LOADING THE DATA	" << std::endl << std::endl;
	std::cout << "=================================================================" << std::endl << std::endl;


	// Setting the data types
	dcomp z;
	double kappa, sigma_11inf, G;
	int nc, m, N;

	// Read the input data from binary file PART 1
	std::cout << "Loading input data" << std::endl;
	double fin[6000];
	std::ifstream input_file("geometry_data.bin", std::ios::in | std::ios::binary);
	input_file.read((char *)&fin, sizeof fin);
	sigma_11inf = fin[0];
	kappa = fin[1];
	G = fin[2];
	nc = (int)fin[3];
	m = (int)fin[4];
	N = (int)fin[5];

	// Declaring the vecotrs
	dcvec z1(nc), z2(nc);
	dvec p(nc), L(nc), mu(nc);
	ddcvec beta(nc, dcvec(m));

	int pos = 5 + 1;
	if (nc > 0)
	{
		for (int ii = 0; ii < nc; ii++)
		{
			int re = pos + ii;
			int im = pos + nc + ii;
			z1[ii] = dcomp(fin[re], fin[im]);
		}
		pos += 2 * nc;
		for (int ii = 0; ii < nc; ii++)
		{
			int re = pos + ii;
			int im = pos + nc + ii;
			z2[ii] = dcomp(fin[re], fin[im]);
		}
		pos += 2 * nc;
		for (int ii = 0; ii < nc; ii++)
		{
			p[ii] = fin[pos + ii];
		}
		pos += nc;
		for (int ii = 0; ii < nc; ii++)
		{
			L[ii] = fin[pos + ii];
		}
		pos += nc;
		for (int ii = 0; ii < nc; ii++)
		{
			mu[ii] = fin[pos + ii];
		}
		pos += nc;
	}
	else
	{
		z1 = { dcomp(0, 0) };
		z2 = { dcomp(0, 0) };
		L = { 0 };
		mu = { 0 };
		beta = { { dcomp(0,0) } };
	}

	// Setting the data types
	double xfrom, xto, yfrom, yto;
	int Nx, Ny, Nw, Ntraj, lvs_traj;
	ddvec grid_11, grid_22, grid_12, grid_1, grid_2, theta_p;
	ddcvec traj_1, traj_2, grid_w, traj_w;
	dvec x_vec, y_vec, x_vecw, y_vecw;
	dcvec xtraj, ytraj;

	// Read the plot data from binary file
	double fplot[800];
	std::ifstream plot_file("plot_data.bin", std::ios::in | std::ios::binary);
	plot_file.read((char *)&fplot, sizeof fplot);
	xfrom = fplot[0];
	xto = fplot[1];
	yfrom = fplot[2];
	yto = fplot[3];
	Nx = (int)fplot[4];
	Ny = (int)fplot[5];
	Nw = (int)fplot[6];
	Ntraj = (int)fplot[7];
	lvs_traj = (int)fplot[8];
	xtraj = { fplot[9] + dcomp(0,1)*fplot[10], fplot[11] + dcomp(0,1)*fplot[12] };
	ytraj = { fplot[13] + dcomp(0,1)*fplot[14], fplot[15] + dcomp(0,1)*fplot[16] };

	std::cout << " -> Complete" << std::endl << std::endl;

	// Displying the plot data in the command window
	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "		THE GEOMETRY DATA	" << std::endl << std::endl;
	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "This is the retrived geometry data:\n";
	std::cout << "sigma_11inf = " << sigma_11inf << std::endl;
	std::cout << "      kappa = " << kappa << std::endl;
	std::cout << "          G = " << G << std::endl;
	std::cout << "         nc = " << nc << std::endl;
	std::cout << "          m = " << m << std::endl;
	std::cout << "          N = " << N << std::endl;
	if (nc > 0)
	{
		std::cout << "         z1 = [";
		for (int ii = 0; ii < nc; ii++)
		{
			if (ii == (nc - 1))
			{
				std::cout << z1[ii];
			}
			else
			{
				std::cout << z1[ii] << ", ";
			}
			if ((ii + 1) % 5 == 0) {
				std::cout << std::endl << "               ";
			}
		}
		std::cout << "]" << std::endl;
		std::cout << "         z2 = [";
		for (int ii = 0; ii < nc; ii++)
		{
			if (ii == (nc - 1))
			{
				std::cout << z2[ii];
			}
			else
			{
				std::cout << z2[ii] << ", ";
			}
			if ((ii + 1) % 5 == 0) {
				std::cout << std::endl << "               ";
			}
		}
		std::cout << "]" << std::endl;
		std::cout << "          p = [";
		for (int ii = 0; ii < nc; ii++)
		{
			if (ii == (nc - 1))
			{
				std::cout << p[ii];
			}
			else
			{
				std::cout << p[ii] << ", ";
			}
			if ((ii + 1) % 5 == 0) {
				std::cout << std::endl << "               ";
			}
		}
		std::cout << "]" << std::endl;
		std::cout << "          L = [";
		for (int ii = 0; ii < nc; ii++)
		{
			if (ii == (nc - 1))
			{
				std::cout << L[ii];
			}
			else
			{
				std::cout << L[ii] << ", ";
			}
			if ((ii + 1) % 5 == 0) {
				std::cout << std::endl << "               ";
			}
		}
		std::cout << "]" << std::endl;
		std::cout << "         mu = [";
		for (int ii = 0; ii < nc; ii++)
		{
			if (ii == (nc - 1))
			{
				std::cout << mu[ii];
			}
			else
			{
				std::cout << mu[ii] << ", ";
			}
			if ((ii + 1) % 5 == 0) {
				std::cout << std::endl << "               ";
			}
		}
		std::cout << "]" << std::endl;
	}
	std::cout << std::endl;

	/* --------------------------------------------------------------------
			Solve the system
	-----------------------------------------------------------------------*/

	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "		INITIALIZING THE SOVLER	" << std::endl << std::endl;
	std::cout << "=================================================================" << std::endl << std::endl;
	auto start_solv = std::chrono::high_resolution_clock::now();
	beta = iterator(1e-6, 300, N, p, sigma_11inf, z1, z2, L, mu, m, nc);
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
	int num;
	if (m < 100)
	{
		num = 201;
	}
	else {
		num = N * 2 + 1;
	}
	dvec T_check_re(num*nc), T_check_im(num*nc);
	double theta_0 = 0.1;
	double deltheta = 2 * pi() / num;
	for (int ii = 0; ii < nc; ii++)
	{
		for (int jj = 0; jj < num; jj++)
		{
			dcomp z, T_check;
			double theta;
			theta = theta_0 + (jj + 1)*deltheta;
			z = z_from_chi(exp(dcomp(0, 1)*theta), z1[ii], z2[ii], L[ii], mu[ii]);
			T_check = T_total(z, sigma_11inf, z1, z2, L, mu, m, nc, beta, -1, ii) + dcomp(0, p[ii]);
			T_check_re[(ii*num) + jj] = abs(real(T_check));
			T_check_im[(ii*num) + jj] = abs(imag(T_check));
		}
	}
	// Median
	size_t n = T_check_re.size() / 2;
	nth_element(T_check_re.begin(), T_check_re.begin() + n, T_check_re.end());
	double error_med_re = T_check_re[n];
	nth_element(T_check_im.begin(), T_check_im.begin() + n, T_check_im.end());
	double error_med_im = T_check_im[n];
	// Mean
	double error_mean_re = 1.0 * std::accumulate(T_check_re.begin(), T_check_re.end(), 0.0) / (n * 2);
	double error_mean_im = 1.0 * std::accumulate(T_check_im.begin(), T_check_im.end(), 0.0) / (n * 2);
	// Max
	double error_max_re = *max_element(T_check_re.begin(), T_check_re.end());
	double error_max_im = *max_element(T_check_im.begin(), T_check_im.end());
	// Print the results
	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "		ERRORS	" << std::endl << std::endl;
	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "     Difference for ts" << std::endl;
	std::cout << "     Maximum = " << error_max_re << std::endl;
	std::cout << "     Mean    = " << error_mean_re << std::endl;
	std::cout << "     Median  = " << error_med_re << std::endl;
	std::cout << "     Difference for tn" << std::endl;
	std::cout << "     Maximum = " << error_max_im << std::endl;
	std::cout << "     Mean    = " << error_mean_im << std::endl;
	std::cout << "     Median  = " << error_med_im << std::endl;

	std::cout << "=================================================================" << std::endl;
	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "		STRESS AND DISPLACEMENT CALCULATOR	" << std::endl << std::endl;
	std::cout << "=================================================================" << std::endl;
	std::cout << "=================================================================" << std::endl << std::endl;


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

	// Estimate the time of the program
	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "		ESTIMATED CALCULATION TIME	" << std::endl << std::endl;
	std::cout << "=================================================================" << std::endl << std::endl;
	int Ne = (Nx + Ny + Ntraj) / 50;
	if (Ne < 8)
	{
		Ne = 8;
	}
	auto start0 = std::chrono::high_resolution_clock::now();
	stress_field(xfrom, xto, yfrom, yto, Ne, Ne, sigma_11inf, z1, z2, L, mu, m, nc, beta, -1);
	auto stop0 = std::chrono::high_resolution_clock::now();
	auto start01 = std::chrono::high_resolution_clock::now();
	principal_stress_trajectories(xfrom, xto, yfrom, yto, xtraj, ytraj, Ne, Ne, sigma_11inf, z1, z2, L, mu, m, nc, beta, -1);
	auto stop01 = std::chrono::high_resolution_clock::now();
	auto duration0 = std::chrono::duration_cast<std::chrono::microseconds>(stop0 - start0);
	long long ms0 = duration0.count();
	auto duration01 = std::chrono::duration_cast<std::chrono::microseconds>(stop01 - start01);
	long long ms01 = duration01.count();
	long long calcs0 = (Nx*Ny + Nx * Ny + Nw * Nw);
	long long calcs01 = (Ntraj * lvs_traj + Ntraj * Nw/2);
	ms0 = (ms0 / (Ne*Ne) * calcs0 + ms01 / (Ne*Ne) * calcs01);
	long long s0, m0, h0;
	std::tie(ms0, s0, m0, h0) = ms_to_time(ms0);
	std::cout << "Estimated calculation time:  ";
	time_print(ms0, s0, m0, h0);
	std::cout << std::endl << std::endl;

	/* --------------------------------------------------------------------
			Calcuating the stress field
	-----------------------------------------------------------------------*/

	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "		COMPUTING THE PLOTS	" << std::endl << std::endl;
	std::cout << "=================================================================" << std::endl << std::endl;

	// Get the Cartisian stress field
	std::cout << "Initiating the stress field calculation\n";
	auto start1 = std::chrono::high_resolution_clock::now();
	std::tie(x_vec, y_vec, grid_11, grid_22, grid_12) = stress_field(xfrom, xto, yfrom, yto, Nx, Ny, sigma_11inf, z1, z2, L, mu, m, nc, beta, -1);
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
	std::tie(grid_1, grid_2, theta_p) = principal_stress_field(xfrom, xto, yfrom, yto, Nx, Ny, sigma_11inf, z1, z2, L, mu, m, nc, beta, -1);
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
	std::tie(traj_1, traj_2) = principal_stress_trajectories(xfrom, xto, yfrom, yto, xtraj, ytraj, Ntraj, lvs_traj, sigma_11inf, z1, z2, L, mu, m, nc, beta, -1);
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
	std::tie(x_vecw, y_vecw, grid_w) = w_field(xfrom, xto, yfrom, yto, Nw, kappa, G, sigma_11inf, z1, z2, L, mu, m, nc, beta);
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
	traj_w = w_trajectories(xfrom, xto, yfrom, yto, Ntraj, Nw, kappa, G, sigma_11inf, z1, z2, L, mu, m, nc, beta);
	auto stop5 = std::chrono::high_resolution_clock::now();
	auto duration5 = std::chrono::duration_cast<std::chrono::microseconds>(stop5 - start5);
	std::cout << "Completed, time taken by function: w_trajectories = ";
	ms1 = duration5.count();
	std::tie(ms1, s1, m1, h1) = ms_to_time(ms1);
	time_print(ms1, s1, m1, h1);
	std::cout << std::endl << std::endl;


	// Displaying the computation time
	std::cout <<  std::endl << std::endl;
	long long mstime = duration1.count() + duration2.count() + duration3.count() + duration4.count() + duration5.count();
	long long stime, mtime, htime;
	std::tie(mstime, stime, mtime, htime) = ms_to_time(mstime);
	std::cout << "Total calculation time: ";
	time_print(mstime, stime, mtime, htime);
	std::cout << std::endl << std::endl;

	/* --------------------------------------------------------------------
			Saving the data as binary files
	-----------------------------------------------------------------------*/
	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "		SAVING THE OUTPUT DATA	" << std::endl << std::endl;
	std::cout << "=================================================================" << std::endl << std::endl;

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

	// Save the plot properties
	dvec dim = { 1.0 * Nx, 1.0 * Ny, 1.0*Nw, 1.0*Ntraj, 1.0*lvs_traj, 1.0*nc};
	const char* pointerdim = reinterpret_cast<const char*>(&dim[0]);
	std::size_t bytesdim = dim.size() * sizeof(dim[0]);
	outfiledim.write(pointerdim, bytesdim);

	// saving the coordinates of the cracks
	dvec fz1_re(nc);
	dvec fz1_im(nc);
	for (int jj = 0; jj < nc; jj++)
	{
		fz1_re[jj] = real(z1[jj]);
		fz1_im[jj] = imag(z1[jj]);
	}
	const char* pointerz1_re = reinterpret_cast<const char*>(&fz1_re[0]);
	std::size_t bytesz1_re = fz1_re.size() * sizeof(fz1_re[0]);
	outfiledim.write(pointerz1_re, bytesz1_re);
	const char* pointerz1_im = reinterpret_cast<const char*>(&fz1_im[0]);
	std::size_t bytesz1_im = fz1_im.size() * sizeof(fz1_im[0]);
	outfiledim.write(pointerz1_im, bytesz1_im);

	dvec fz2_re(nc);
	dvec fz2_im(nc);
	for (int jj = 0; jj < nc; jj++)
	{
		fz2_re[jj] = real(z2[jj]);
		fz2_im[jj] = imag(z2[jj]);
	}
	const char* pointerz2_re = reinterpret_cast<const char*>(&fz2_re[0]);
	std::size_t bytesz2_re = fz2_re.size() * sizeof(fz2_re[0]);
	outfiledim.write(pointerz2_re, bytesz2_re);
	const char* pointerz2_im = reinterpret_cast<const char*>(&fz2_im[0]);
	std::size_t bytesz2_im = fz2_im.size() * sizeof(fz2_im[0]);
	outfiledim.write(pointerz2_im, bytesz2_im);

	// Close the output files
	outfile.close();
	outfiledim.close();

	std::cout << " -> Complete" << std::endl << std::endl;

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

	// Create the log file
	std::ofstream logfile;
	logfile.open("log.txt");
	logfile << "=================================================================" << std::endl << std::endl;
	logfile << "    STRESS AND DISPLACEMENT CALCULATOR - LOG FILE                 " << std::endl << std::endl;
	logfile << "=================================================================" << std::endl;
	logfile << "This is the log file for computation completed on: " << str_time << std::endl;
	logfile << "The program took: ";
	if (hours < 10 && minutes < 10 && seconds < 10)
	{
		logfile << "0" << hours << ":0" << minutes << ":0" << seconds << ":" << mseconds << std::endl;
	}
	else if (hours < 10 && minutes < 10)
	{
		logfile << "0" << hours << ":0" << minutes << ":" << seconds << ":" << mseconds << std::endl;
	}
	else if (hours < 10)
	{
		logfile << "0" << hours << ":" << minutes << ":" << seconds << ":" << mseconds << std::endl;
	}
	else
	{
		logfile << "" << hours << ":" << minutes << ":" << seconds << ":" << mseconds << std::endl;
	}
	logfile << std::endl << std::endl;
	logfile << "=================================================================" << std::endl;
	logfile << "    DATA FILES								           " << std::endl;
	logfile << "=================================================================" << std::endl;
	logfile << "'data.bin'     [x_vec, y_vec, x_vecw, y_vecw, grid_11, grid_22, grid_12, grid_1, grid_2, theta_p, (traj_1.real, traj_1.imag), (traj_2.real, traj_2.imag), (grid_w.real, grid_w.imag), (traj_w.real, traj_w.imag)]\n";
	logfile << "'dim_data.bin' [Nx, Ny, Nw, Ntraj, lvs_traj, nc, z1, z2]\n";
	logfile << std::endl << std::endl;
	logfile << "=================================================================" << std::endl;
	logfile << "    PLOT DATA                                          " << std::endl;
	logfile << "=================================================================" << std::endl;
	logfile << "x from: " << xfrom << " to " << xto << std::endl;
	logfile << "y from: " << yfrom << " to " << yto << std::endl;
	logfile << "x resolution: " << Nx << std::endl;
	logfile << "y resolution: " << Ny << std::endl;
	logfile << "Total number of points: " << Nx * Ny << std::endl;
	logfile << "Number of trajectory levels: " << lvs_traj << std::endl;
	logfile << "Number of steps in trajectories: " << Ntraj;
	logfile << "sigma_1 starting line from: " << xtraj[0] << " to " << xtraj[1] << std::endl;
	logfile << "sigma_2 starting line from: " << ytraj[0] << " to " << ytraj[1] << std::endl;
	logfile << "x and y quiver resolution: " << Nw;
	logfile << std::endl << std::endl;
	logfile << "=================================================================" << std::endl;
	logfile << "    VARIABLES										   " << std::endl;
	logfile << "=================================================================" << std::endl;
	logfile << "Elasticity:\n";
	logfile << "      kappa = " << kappa << std::endl;
	logfile << std::endl;
	logfile << "Unifrom stress state:\n";
	logfile << "sigma_11inf = " << sigma_11inf << std::endl;
	logfile << std::endl;
	logfile << std::endl;
	if (nc > 0)
	{
		logfile << "Cracks:\n";
		logfile << "         nc = " << nc << std::endl;
		logfile << "          m = " << m << std::endl;
		logfile << "         z1 = [";
		for (int ii = 0; ii < nc; ii++)
		{
			if (ii == (nc - 1))
			{
				logfile << z1[ii];
			}
			else
			{
				logfile << z1[ii] << ", ";
			}
		}
		logfile << "]" << std::endl;
		logfile << "         z2 = [";
		for (int ii = 0; ii < nc; ii++)
		{
			if (ii == (nc - 1))
			{
				logfile << z2[ii];
			}
			else
			{
				logfile << z2[ii] << ", ";
			}
		}
		logfile << "]" << std::endl;
		logfile << "          L = [";
		for (int ii = 0; ii < nc; ii++)
		{
			if (ii == (nc - 1))
			{
				logfile << L[ii];
			}
			else
			{
				logfile << L[ii] << ", ";
			}
		}
		logfile << "]" << std::endl;
		logfile << "         mu = [";
		for (int ii = 0; ii < nc; ii++)
		{
			if (ii == (nc - 1))
			{
				logfile << mu[ii];
			}
			else
			{
				logfile << mu[ii] << ", ";
			}
		}
		logfile << "]" << std::endl;

		for (int ii = 0; ii < nc; ii++)
		{
			logfile << "    beta[" << ii << "] = [";
			for (int jj = 0; jj < m; jj++)
			{
				if (jj == (m - 1))
				{
					logfile << beta[ii][jj];
				}
				else
				{
					logfile << beta[ii][jj] << ", ";
				}
			}
			logfile << "]" << std::endl;
		}
		logfile << std::endl << std::endl;
	}

	logfile.close();
	std::cout << "Program finnished after ";
	time_print(mseconds, seconds, minutes, hours);
	std::cout << std::endl;
	std::cout << "Output data saved to binary files: data.bin and dim_data.bin" << std::endl;
	std::cout << "For log data see log.txt" << std::endl;

	return 0;

	return 0;
}