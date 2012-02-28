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

#include <omd/forcekernel.h>

// VerletList required...

class StillingerWeber: public ForceKernel {

	// --- parameters read from file ---
	string paramfile;
	double eps;
	double bigA;
	double bigB;
	double powerp;
	double powerq;
	double sigma;
	double alpha;
	double lambda;
	double gamma;
	double icos0;

	// --- local registers ---
	double c1,c2,c3,c4,c5,c6,c7,c8,c9,c10;
	double zc,ze0,ze1,ze2,ze3;
	double sx0,sx1,sa0,sa1,sa2,sa3;
	
	bool  usezbl;
	class VerletListFull* Verlet;
	
public:

	StillingerWeber(string material);
	void Init(class MDSystem* WorkSys);
	void ReadParameter();

	void ZBL(Atom& at, Atom& to,
					 double r,
					 double dx,
					 double dy,
					 double dz);

	void Spline(Atom& at, Atom& to,
					 double r,
					 double dx,
					 double dy,
					 double dz);
	
	void TwoBodyTerm(Atom& at, Atom& to,
					 double r,
					 double dx,
					 double dy,
					 double dz);
	
	void ThreeBodyTerm(Atom& at0, Atom& at1, Atom& at2,
					   double r1, // at <- first_nb
					   double dx1,
					   double dy1,
					   double dz1,
					   double r2, // at <- second_nb (iterated)
					   double dx2,
					   double dy2,
					   double dz2);
	
	void ComputeHalf(Atom& at, Atom& to);
	void ComputeFull(Atom& at, Atom& to);
	void PrintInfo(ostream& ost);

};

#endif
