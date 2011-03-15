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
 */

#ifndef _VERLET_LIST_HPP_
#define _VERLET_LIST_HPP_

#include <fstream>
#include <cmath>
#include <OMD/gadgets.h>

class VerletList: public MDIterator {
	int UpdatePeriod;
	int Step;
	int U, V, W;
	int NeighSize;
	int *Link;
	int *NeighborList;
	int *NeighborIndex;
	double VerletRadius;
	
public:
	VerletList();
	virtual ~VerletList();

	void CellNumber(int, int&, int&, int&);
	void Reallocate(int NewSize);
	void RefreshTable();
	void Init(SimSystem* WorkSys);
	void PrintInfo() {
		cerr << "ID." << id << " " << Name << " (Radius=" << VerletRadius 
		     << "; Update=" << UpdatePeriod << " steps)\n";
	}
	void PreCalculation() {
		if (Step>=UpdatePeriod)	{RefreshTable(); Step=0;}
		Step++;
	}
	void DumpNeighborIndex();

};

#endif
