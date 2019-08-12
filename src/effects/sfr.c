#include "../proto.h"
#include "../globals.h"
#include "effects.h"

void Sfr(int ipart, double *result)
{
#ifdef SFR
	result[0] = Gas[ipart].Sfr;
#endif
	return;
}
