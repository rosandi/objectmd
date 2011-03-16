/* omdlib
 ********************************************************
 *
 * ObjectMD
 * Molecular Dynamics Class Library
 * Version 3.0 (2009,2011) 
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
 *       Quadratic repulsive potential
 *
*/

#include <cstdio>
#include <cmath>
#include <omd/gadgets.h>

// FIXME! Update needed
// FIXME! give reference

#define QUAD_R_CUT 5.1

#define Q(x)  exp(qa+x*(qb+qc*x))
#define dQ(x)  2.0*qc*qb*x*Q(x)

class Quadratic: public ForceKernel {
	OMD_FLOAT qa, qb, qc;
	public:	

	Quadratic(OMD_FLOAT a, OMD_FLOAT b, OMD_FLOAT c) {
		qa=a; qb=b; qc=c;
		strcpy(Name, "Quadratic O-Pt Potential");
		CutRadius=QUAD_R_CUT;
	}
	
	void Compute(int at, int to) {
		OMD_FLOAT dx, dy, dz;
		OMD_FLOAT r = sqrt(CalcSqrDistance(at, to, dx, dy, dz));
		OMD_FLOAT ff;
        if (r < CutRadius) {			
			potential=Q(r);	drpot=dQ(r); ff=drpot/r;
	        RETURN_FORCE(dx*ff, dy*ff, dz*ff);
		}
	}
};
