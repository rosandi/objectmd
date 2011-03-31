/* omdlib
 ********************************************************
 *
 * ObjectMD
 * Molecular Dynamics Class Library
 * Version 3.1 (2009,2011)
 *
 * Released under:
 *       GNU GENERAL PUBLIC LICENSE Version 3
 *       Copyright (C) 2005  Yudi Rosandi
 *
 * Read license.txt file in the source's root directory.
 *
 ********************************************************
 * Stillinger-Weber Potential
 *
*/
 
#ifndef _T_SW_H_
#define _T_SW_H_

#include <omd/forcekernel.hpp>

// VerletList required...

class StillingerWeber: public ForceKernel {

	// --- parameters read from file ---
	string paramfile;
	OMD_FLOAT eps;
	OMD_FLOAT bigA;
	OMD_FLOAT bigB;
	OMD_FLOAT powerp;
	OMD_FLOAT powerq;
	OMD_FLOAT sigma;
	OMD_FLOAT alpha;
	OMD_FLOAT lambda;
	OMD_FLOAT gamma;

	// --- local registers ---

	double sigma_gamma, lambda_eps, lambda_eps2;
	double c1,c2,c3,c4,c5,c6;
	class VerletListFull* Verlet;
	
public:

	StillingerWeber(string material);
	virtual ~StillingerWeber();
	void Init(class MDSystem* WorkSys);
	void ReadParameter();
	void TwoBodyTerm(Atom& at, Atom& to,
					 OMD_FLOAT r,
					 OMD_FLOAT dx,
					 OMD_FLOAT dy,
					 OMD_FLOAT dz);
	
	void ThreeBodyTerm(Atom& at0, Atom& at1, Atom& at2,
					   OMD_FLOAT r1, // at <- first_nb
					   OMD_FLOAT dx1,
					   OMD_FLOAT dy1,
					   OMD_FLOAT dz1,
					   OMD_FLOAT r2, // at <- second_nb (iterated)
					   OMD_FLOAT dx2,
					   OMD_FLOAT dy2,
					   OMD_FLOAT dz2);
	
	void ComputeHalf(Atom& at, Atom& to);
	void ComputeFull(Atom& at, Atom& to);
	void PrintInfo(ostream& ost);

};

#endif
