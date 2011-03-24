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

/**
 @ingroup iterator
 @brief Verlet neighbor list (Iterator)
 
 This class is the verlet neighbor list iterator algorithm.
 The calculation is about twice as fast compares to the NeighborCell iterator.
 No periodic boundary condition is considered. Use this class with MDSystemGrid
 for periodic boundary simulation.
 
 Parameters:
 - @e verlet.update : the update period of the neighbor list. This parameter is ignored
                      in parallel OMD, and the update period is controlled by MDSystemGrid.
 - @e verlet.rebuild : the period to rebuild the neighbor grid. This is zero by default, which
					means no rebuilding.
 - @e verlet.radtole : the radius tolerance of the neighbor grid. In parallel MD this is ignored
                       and the radius is defined by MDSystemGrid.
 - @e verlet.nmean : the estimate value of number of neighbors of an atom. The default is 120. 

**/

class VerletList: public MDIterator {
	int U, V, W;
	int nmean; 
	int NeighSize;
	int AllocSize;
	OMD_FLOAT VerletRadius;
	OMD_FLOAT CutRadius;
	int   *Link;
	int *NeighborList;
	int *NeighborIndex;
	SysBox    Box;
	
	// iteration uses class these variables
	int at_idx;
	int to_idx;
	int nl_idx;
	int ls_start;
	int ls_end;

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
	void GetIterationVariables(int &at, int &to, 
							   int& neiglist_index, int& list_start, int& list_end)
	{at=at_idx;to=to_idx;neiglist_index=nl_idx;list_start=ls_start;list_end=ls_end;}
	void Dump(string fname);
	void PrintInfo(ostream& ost);

};

#endif
