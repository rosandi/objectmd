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
#include <omd/forcekernel.hpp>

//!Morse force kernel

/**
 *
*/

#define V(R)  D*(exp(-2.0*alp*(R-r0))-2.0*exp(-alp*(R-r0)))
#define dV(R) -2.0*D*alp*(exp(-2.0*alp*(R-r0))-exp(-alp*(R-r0)))
	
class Morse: public ForceKernel {
	
public:	
	OMD_FLOAT alp, r0, D;

	Morse(OMD_FLOAT Alpha, OMD_FLOAT R_0, OMD_FLOAT _D, OMD_FLOAT R_CUT=5.1) {
		strcpy(Name, "Morse-Potential");
		alp=Alpha; r0=R_0; D=_D;
		CutRadius=R_CUT;
	}
	
	void Compute(int at, int to) {
		OMD_FLOAT dx, dy, dz;
		OMD_FLOAT r = sqrt(CalcSqrDistance(at, to, dx, dy, dz));
		OMD_FLOAT ff;
        	if (r < CutRadius) {
			potential = V(r); 
			ff=(-(drpot=dV(r)))/r;
			RETURN_FORCE(dx*ff, dy*ff, dz*ff);
		}
	}
};
