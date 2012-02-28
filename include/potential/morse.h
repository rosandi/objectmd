/* omdlib
 ********************************************************
 *
 * ObjectMD
 * Molecular Dynamics Class Library
 * Version 2.0 (2009)
 *
 * Released under:
 *       GNU GENERAL PUBLIC LICENSE Version 3
 *       Copyright (C) 2005  Yudi Rosandi
 *
 * Read license.txt file in the source's root directory.
 *
 ********************************************************
 * ObjectMD header file
 *
 * Morse potential implementation
 *
 *
*/

// FIXME! update needed

#include <cstdio>
#include <cmath>
#include <omd/forcekernel.h>

//!Morse force kernel

/**
 *
*/

#define V(R)  D*(exp(-2.0*alp*(R-r0))-2.0*exp(-alp*(R-r0)))
#define dV(R) -2.0*D*alp*(exp(-2.0*alp*(R-r0))-exp(-alp*(R-r0)))
	
class Morse: public ForceKernel {
	
public:	
	double alp, r0, D;

	Morse(double Alpha, double R_0, double _D, double R_CUT=5.1) {
		strcpy(Name, "Morse-Potential");
		alp=Alpha; r0=R_0; D=_D;
		CutRadius=R_CUT;
	}
	
	void ComputeHalf(int at, int to) {
		double dx, dy, dz;
		double r = sqrt(CalcSqrDistance(at, to, dx, dy, dz));
		double ff;
        	if (r < CutRadius) {
			potential = V(r); 
			ff=(-(drpot=dV(r)))/r;
			RETURN_FORCE(dx*ff, dy*ff, dz*ff);
		}
	}
};
