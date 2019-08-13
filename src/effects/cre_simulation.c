#include "../proto.h"
#include "../globals.h"
#include "effects.h"

#define CNST_ME  9.10953e-28
#define CNST_MP  1.6726e-24
#define CNST_C   2.9979e10
#define CR_DSlope 1.0e-6
//
// static double Emax, Emin, pstep;
// static double *p_bounds;

double cre_spectrum_simulation( double E, int ipart )
{
	// double hbound_proper;
	// double N;
	// double norm;
	//
	// if (E < Emin || E > Emax) // out of range
	// 	return 0;
	//
	// int ibin = floor(log10(E / Emin) / pstep);
	//
	// if ( p_bounds[ibin+1] > Gas[ipart].CReCut )
	// 	hbound_proper = Gas[ipart].CReCut;
	// else
	// 	hbound_proper = p_bounds[ibin+1];
	//
	// norm = Gas[ipart].CReNorm[ibin] * 1.e20;
	//
	// N = density_integral( p_bounds[ipart], hbound_proper,
	// 					    norm , Gas[ipart].CReSlope[ibin],
	// 				  	     P[ipart].Rho);

	return ;
}

/* spectrum input */
void set_simulated_factors()
{
	// store momentum bin boundaries
	// int Nbound;
	// size_t nBytes;
	//
	// Emin = Param.CR_Emin;
	// Emax = Param.CR_Emax;
	//
	// nBytes = (BP_REAL_CRs + 1) * sizeof(*p_bounds);
	// p_bounds = Malloc(nBytes);
	//
	// pstep = log10(Emax/Emin) / BP_REAL_CRs;
	// for( Nbound = 0; Nbound <= BP_REAL_CRs; Nbound++ )
	//   {
	// 	// dimensionless e momenta into cgs
	// 	p_bounds[Nbound] = Param.CR_Emin * CNST_ME * CNST_C * pow(10.0,(pstep*Nbound));
	// 	printf("%i\t%g\n", Nbound, p_bounds[Nbound]);
	//   }

}
