/* Computes synchroton emission from an arbitrary CR electron spectrum,
 * set in synchro_factors (Donnert & Brunetti 2014)
 * */

#include "../proto.h"
#include "../globals.h"
#include "effects.h"

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_synchrotron.h>
#include <gsl/gsl_spline.h>


#define NSTEP 128 // <1% rel error with power-law N(E)
#define TABLE_SIZE NSTEP

#define X_MAX 70.0 // -> F(x) < 2e-6
#define X_MIN 1e-20// -> F(x) < 5e-6

#define THETA_MAX 0.5*pi

double gsl_sf_synchrotron_1(double);	// GSL synchro kernel func
double gsl_sf_synchrotron_2(double);
static inline double synchro_kernel(double, double *, double *);
static void prepare_kernel();
static double stokes2angle(double Q, double U);
extern void add_Ipol_chi_Pi();

// BP_REAL_CRs
// static inline double bp_energy_integral( double, double, double, float, float );
// static inline double bp_density_integral(double, double, double, float, float );
#ifdef BP_REAL_CRs
#define CNST_ME  9.10953e-28
#define CNST_MP  1.6726e-24
#define CNST_C   2.9979e10
#define CR_DSlope 1.0e-6
#define  ELECTRONCHARGE  4.8032e-10
#define erg2eV 1.6022e-12

static double *p_bounds;
static inline double bp_energy_integral( double, double, double, float, float );
static inline double bp_density_integral(double, double, double, float, float );
#endif

static double (*spectrum_ptr) (double, int);

static double X_table[TABLE_SIZE] = { 0 },Kernel_Table1[TABLE_SIZE] = { 0 },
			  Kernel_Table2[TABLE_SIZE] = { 0 },
			  nu_c_prefac, j_nu_prefac, E_cntr_prefac;

static gsl_spline *Synchro_Spline1 = NULL, *Synchro_Spline2 = NULL;
static gsl_interp_accel *Acc[2] = { NULL };
#pragma omp threadprivate(Synchro_Spline1, Synchro_Spline2, Acc)


/*
 * Return Synchrotron brightness Integrate over energy Ee
 * (and pitch angle theta) using Trapezoidal rule.
 * see Longair, High Energy Astophysics 1994
 */

void synchrotron(int ipart, double *j_nu)
{
	const double nu = Param.Freq;

	double B = 0;		// calc B strength

	if (Param.SynchroPAngInt)
		B = length3(Gas[ipart].Bfld);
	else
		B = length2(Gas[ipart].Bfld);

	if (B == 0)		// catch boring ones
		return;

	for (int i = 0; i < MAXIMAGES; i++)
		j_nu[i] = 0;

	double E[NSTEP] = { 0 }, dE[NSTEP] = { 0 }, x[NSTEP] = { 0 },
		   F[NSTEP] =  { 0 }, F_mid[NSTEP] = { 0 };

#ifdef POLARISATION
	double F_para[NSTEP] = { 0 }, F_mid_para[NSTEP] = { 0 },
		   F_orth[NSTEP] = { 0 }, F_mid_orth[NSTEP] = { 0 };
#endif

	double E_min = sqrt(nu/(nu_c_prefac * B * X_MAX)); // [me*c]
	double E_max = sqrt(nu/(nu_c_prefac * B * X_MIN));

	double di = log(E_max/E_min)/(NSTEP-1);

	x[0] = X_MAX;
	E[0] = E_min;
	dE[0] = E[0] - E_min * exp(-1*di);

	for (int i = 1; i < NSTEP; i++) {

		E[i] = E_min * exp(di * i);
		dE[i] = E[i] - E[i-1];
		x[i] = nu / (nu_c_prefac * p2(E[i]) * B);
	}

	double j_para = 0, j_orth = 0;

	for (int i = 1; i < NSTEP; i++) { // Simpson rule

		double K_orth = 0, K_para = 0;

		double N_E = (*spectrum_ptr)(E[i], ipart);
		double K = synchro_kernel(x[i], &K_orth, &K_para);

		F[i] = N_E * K;

#ifdef POLARISATION
		F_orth[i] = N_E * K_orth;
		F_para[i] = N_E * K_para;
#endif

		double x_mid = 0.5 * (x[i-1] + x[i]);
		double e_mid = 0.5 * (E[i-1] + E[i]);

		double N_E_mid = (*spectrum_ptr)(e_mid, ipart);
		double K_mid = synchro_kernel(x_mid, &K_orth, &K_para);

		F_mid[i] = N_E_mid * K_mid;

#ifdef POLARISATION
		F_mid_orth[i] = N_E_mid * K_orth;
		F_mid_para[i] = N_E_mid * K_para;
#endif

		j_nu[0] += dE[i] / 6.0 * (F[i] + F[i-1] + 4*F_mid[i]);

#ifdef POLARISATION
		j_orth += dE[i] / 6.0 * (F_orth[i] + F_orth[i-1] + 4*F_mid_orth[i]);
		j_para += dE[i] / 6.0 * (F_para[i] + F_para[i-1] + 4*F_mid_para[i]);
#endif
	}

	j_nu[0] *= j_nu_prefac * B;

#ifdef POLARISATION
	j_para *= 0.5 * j_nu_prefac * B;
	j_orth *= 0.5 * j_nu_prefac * B;

	const double bx = Gas[ipart].Bfld[0];
	const double by = Gas[ipart].Bfld[1];

	double sin_2chi = 0, cos_2chi = 0;

	if (bx || by) {

		sin_2chi = -2.0 * bx * by / (bx * bx + by * by);
		cos_2chi = (bx * bx - by * by) / (bx * bx + by * by);
	}

	double j_pol = j_orth - j_para;

	j_nu[1] = j_pol * cos_2chi;	// Q
	j_nu[2] = j_pol * sin_2chi;	// U
#endif

	return;
}


void synchrotron_bp(int ipart, double *j_nu)
{
	const double nu = Param.Freq;

	double B = 0;		// calc B strength
	double GeV_factor = 1.e-9;

	if (Param.SynchroPAngInt)
		B = length3(Gas[ipart].Bfld);
	else
		B = length2(Gas[ipart].Bfld);

	if (B == 0)		// catch boring ones
		return;

	for (int i = 0; i < MAXIMAGES; i++)
		j_nu[i] = 0;

	//double F[BP_REAL_CRs] =  { 0 };
	double F = 0.0;
	double x;
	double hbound_proper;

#ifdef POLARISATION
	// double F_para[BP_REAL_CRs] = { 0 },
	// 	   F_orth[BP_REAL_CRs] = { 0 };
	double F_para = 0.0,
		   F_orth = 0.0;
#endif

	double j_para = 0, j_orth = 0;

	for (int i = 0; i < BP_REAL_CRs; i++) { // Simpson rule

		double K_orth = 0, K_para = 0, K = 0;

		if ( p_bounds[i+1] > Gas[ipart].CReCut )
			hbound_proper = Gas[ipart].CReCut;
		else
			hbound_proper = p_bounds[i+1];

		double norm = Gas[ipart].CReNorm[i] * 1.e20;

		double E_E = bp_energy_integral( p_bounds[i], hbound_proper,
							             norm, Gas[ipart].CReSlope[i],
						  	             P[ipart].Rho )
					 * P[ipart].Mass * Unit.Mass * erg2eV * GeV_factor;

		double N_E = bp_density_integral( p_bounds[i], hbound_proper,
							              norm, Gas[ipart].CReSlope[i],
						  	              P[ipart].Rho )
					 * P[ipart].Mass * Unit.Mass;

		x = nu / (nu_c_prefac * p2(E_E) * B);

		if ( (E_E != 0.0) && (x <= X_MAX) && (x >= X_MIN) )
		{
		//
		// printf("%i\t%g\t%g\t%g\t%g\t%g\n", ipart, x,
		// 								   E_E, N_E,
		// 								   norm,
		// 							   	   Gas[ipart].CReSlope[i]);

		K = synchro_kernel(x, &K_orth, &K_para);
		}

		//F[i] = N_E * K;
		F = N_E * K;

#ifdef POLARISATION
		// F_orth[i] = N_E * K_orth;
		// F_para[i] = N_E * K_para;
		F_orth = N_E * K_orth;
		F_para = N_E * K_para;
#endif


		// j_nu[0] += F[i];
		j_nu[0] += F;

#ifdef POLARISATION
		// j_orth += F_orth[i];
		// j_para += F_para[i];
		j_orth += F_orth;
		j_para += F_para;
#endif
	}

	j_nu[0] *= j_nu_prefac * B;

#ifdef POLARISATION
	j_para *= 0.5 * j_nu_prefac * B;
	j_orth *= 0.5 * j_nu_prefac * B;

	const double bx = Gas[ipart].Bfld[0];
	const double by = Gas[ipart].Bfld[1];

	double sin_2chi = 0, cos_2chi = 0;

	if (bx || by) {

		sin_2chi = -2.0 * bx * by / (bx * bx + by * by);
		cos_2chi = (bx * bx - by * by) / (bx * bx + by * by);
	}

	double j_pol = j_orth - j_para;

	j_nu[1] = j_pol * cos_2chi;	// Q
	j_nu[2] = j_pol * sin_2chi;	// U
#endif

	return;
}


/* choose CRe spectrum according to Spec */
void set_synchro_factors()
{
	/* if you change this, include prepare_synchrotron()  */
	prep_func_ptr = NULL;

	switch (Param.Effect_Flag) {

	case 0:		// High Energy Approx. Secondary (Brunetti 05)
		set_brunetti_factors();
		spectrum_ptr = &cre_secondaries_brunetti;
		break;

	case 1:		// Read spectrum from file N(r)
		prep_func_ptr = &set_tabulated_factors;
		spectrum_ptr = &cre_spectrum_tabulated;
		break;

	case 2:		// analytic power law spectrum
		set_powerlaw_factors(0);
		spectrum_ptr = &cre_spectrum_powerlaw;
		break;

	case 3:
		printf("Calculating spectrum from BP_REAL_CRs\n");
		//set_simulated_factors();
		//spectrum_ptr = &cre_spectrum_simulation;
		break;

	case 4:		// power law with norm from sim
		printf("Calculating spectrum from simulation\n");
		set_powerlaw_factors(1);
		spectrum_ptr = &cre_spectrum_powerlaw;
		break;

	case 10:		// Read compressed Spectrum
		prep_func_ptr = &setup_decompression;
		spectrum_ptr = &cre_spectrum_compressed;
		break;

	default:
		Assert(0, "Selected Effect_Flag %d not handled",
		       Param.Effect_Flag);
		break;

	}

	prepare_kernel();

	E_cntr_prefac = m_e*c*c * sqrt(2.0/3.0 * (2*pi * m_e*c)/e /0.29);
	nu_c_prefac = 3.0 * e / (4.0 * pi * pow(m_e, 3) * pow(c, 5));
	j_nu_prefac = e * e * e * sqrt(3) / (m_e * c * c);

#ifdef BP_REAL_CRs
	// store momentum bin boundaries
	int Nbound;
	size_t nBytes;
	double pstep;

	nBytes = (BP_REAL_CRs + 1) * sizeof(*p_bounds);
	p_bounds = Malloc(nBytes);

	pstep = log10(Param.CR_Emax/Param.CR_Emin) / BP_REAL_CRs;
	for( Nbound = 0; Nbound <= BP_REAL_CRs; Nbound++ )
	  {
		// dimensionless e momenta into cgs
		p_bounds[Nbound] = Param.CR_Emin * CNST_ME * CNST_C * pow(10.0,(pstep*Nbound));
		printf("%i\t%g\n", Nbound, p_bounds[Nbound]);
	  }
	//
	// nu_c_prefac = 3.0 * ELECTRONCHARGE / (4.0 * pi * pow(CNST_ME, 3) * pow(CNST_C, 5));
	// j_nu_prefac = ELECTRONCHARGE * ELECTRONCHARGE * ELECTRONCHARGE * sqrt(3) / (CNST_ME * CNST_C * CNST_C);
#endif

	return;
}

/* Returns Synchrotron Kernel & Polarisations for x using
 * limit formulas @ small and large energies. In between the Integral is
 * tabulated.
 */

static inline double synchro_kernel(double x, double *orthPol, double *paraPol)
{
	double Fx = gsl_spline_eval(Synchro_Spline1, x, Acc[0]);

#ifdef POLARISATION

	double Gx = gsl_sf_synchrotron_2(x); //gsl_spline_eval(Synchro_Spline2, x, Acc[1]);

	*orthPol = (Fx + Gx);	// Rybicki & Lightman p179
	*paraPol = (Fx - Gx);

#endif

	return Fx;
}

void prepare_kernel()
{
	double dx = log(X_MAX*1.1/X_MIN/0.9)/(TABLE_SIZE-1);

	#pragma omp parallel
	{

	#pragma omp for
	for (int i = 0; i < TABLE_SIZE; i++) {

		X_table[i] = X_MIN*0.9 * exp(dx*i);

		Kernel_Table1[i] = gsl_sf_synchrotron_1(X_table[i]);
		Kernel_Table2[i] = gsl_sf_synchrotron_2(X_table[i]);
	}

	Acc[0] = gsl_interp_accel_alloc();
	Acc[1] = gsl_interp_accel_alloc();

	Synchro_Spline1 = gsl_spline_alloc(gsl_interp_cspline, TABLE_SIZE);
	Synchro_Spline2 = gsl_spline_alloc(gsl_interp_cspline, TABLE_SIZE);

	gsl_spline_init(Synchro_Spline1, X_table, Kernel_Table1, TABLE_SIZE);
	gsl_spline_init(Synchro_Spline2, X_table, Kernel_Table2, TABLE_SIZE);

	} // omp parallel

	return;
}

/*
 * Convert stokes parameter to polarisation angle.
 * angle is defined [-pi,pi]
 */

static double stokes2angle(double Q, double U)
{
	return atan2(U, Q) * 0.5;
}

extern void add_Ipol_chi_Pi()
{
	const int npix2 = p2(Param.XYPix);
	size_t i, iQ = npix2, iU = 2 * npix2, iIpol = 3 * npix2,
	    iChi = 4 * npix2, iPi = 5 * npix2;

	rprintf("\nAdding <3> images: Ipol, chi, Pi \n\n");

	Effect.Nimage = 6;

	size_t nBytes = Effect.Nimage * npix2 * sizeof(*Image);

	Image = Realloc(Image, nBytes);

	/* Fill new images */
	for (i = 0; i < npix2; i++) {
		/* pol. Intensity */
		Image[iIpol] = sqrt(Image[iQ] * Image[iQ] + Image[iU] * Image[iU]);
		/* pol. Angle */
		Image[iChi++] = stokes2angle(Image[iQ++], Image[iU++]);
		/* pol. Degree */
		Image[iPi++] = Image[iIpol++] / Image[i];
	}

	return;
}


#ifdef BP_REAL_CRs

static inline double bp_density_integral(double bound_low, double bound_up,
									  double norm, float slope,
									  float density_in)
	/* density integral (eq. 9 M01) */
{
	double nb;
	double density, density2;
	double slope_var;

	nb = 4.0 * M_PI * norm * bound_low * bound_low * bound_low
	     / ( density_in );

	density = nb * (pow((bound_up/bound_low),(3.0-slope)) - 1.0) / (3.0 - slope);
	if( slope > (3.0 - CR_DSlope) && slope < (3.0 + CR_DSlope) )
	  {
	    slope_var = ( slope - 3.0 ) / CR_DSlope;
	    density2 = nb * log(bound_up/bound_low);
	    if( slope_var != 0.0 )
	      density = density * slope_var + density2 * (1.0 - slope_var);
	    else
	      density = density2;
	  }

	return ( density );
}

static inline double bp_energy_integral(double bound_low, double bound_up,
									    double norm, float slope,
								 	    float density_in)
	/* energy integral (eq.21 M01) */
{
	double en;
	double energy, energy2;
	double slope_var;

	en = 4.0 * M_PI * CNST_C * norm * pow(bound_low,4.0)
	     / ( density_in );

	energy = en * (pow((bound_up/bound_low),(4.0-slope)) - 1.0) / (4.0 - slope);

	if( slope > (4.0 - CR_DSlope) && slope < (4.0 + CR_DSlope) )
	  {
	    slope_var = ( slope - 4.0 ) / CR_DSlope;
	    energy2 = en * log(bound_up/bound_low);
	    if( slope_var != 0.0 )
	      energy = energy * slope_var + energy2 * (1.0 - slope_var);
	    else
	      energy = energy2;
	  }

	return ( energy );
}

#endif

#undef X_MAX
#undef THETA_MAX
#undef NSTEP
