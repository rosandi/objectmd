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
	OMD_INT U, V, W;	
	OMD_INT   nmean;
	OMD_SIZET NeighSize;
	OMD_SIZET AllocSize;
	OMD_FLOAT VerletRadius;
	OMD_INT   *Link;
	OMD_SIZET *NeighborList;
	OMD_SIZET *NeighborIndex;
	SysBox    Box;

public:
	VerletList();
	virtual ~VerletList();
	
	void ReadParameter();
	bool CheckParameter();
	void Update();
	void Refresh();
	void CellNumber(OMD_SIZET, OMD_SIZET&, OMD_SIZET&, OMD_SIZET&);
	void Init(MDSystem* WorkSys);
	void CalculateBox();
	void PreIntegration();
	void Iterate(MDGadget* IteratedClass, bool force_update=false);
	void GetNeighborIndex(OMD_SIZET ni, OMD_SIZET& start, OMD_SIZET& end);
	int  GetNeighbor(OMD_SIZET ni);
	void Dump(string fname);
	void PrintInfo();	

};

#endif
