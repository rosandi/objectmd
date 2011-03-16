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

#include <omd/iterator.hpp>

class VerletList: public MDIterator {
	int U, V, W;
	int nmean;
	int NeighSize;
	int AllocSize;
	OMD_FLOAT VerletRadius;
	int   *Link;
	int *NeighborList;
	int *NeighborIndex;
	SysBox    Box;

public:
	VerletList();
	virtual ~VerletList();
	
	void ReadParameter();
	bool CheckParameter();
	void Update();
	void Refresh();
	void CellNumber(int, int&, int&, int&);
	void Init(MDSystem* WorkSys);
	void CalculateBox();
	void PreIntegration();
	void Iterate(MDGadget* IteratedClass, bool force_update=false);
	void GetNeighborIndex(int ni, int& start, int& end);
	int  GetNeighbor(int ni);
	void Dump(string fname);
	void PrintInfo(ostream& ost);

};

#endif
