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

#ifndef _VERLET_LIST_FULL_HPP_
#define _VERLET_LIST_FULL_HPP_

#include <omd/iterator.h>

namespace omd {

struct NeighborList {
	// index=0 for accessing full neighbor list
	int  start;
	int  end;
	int  size;
	int* list;
	int* flag;
};

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

class VerletListFull: public MDIterator {
	int U, V, W;
	int nmean; 
	int AllocSize;
	double VerletRadius;
	double CutRadius;
	int *Link;
	
	NeighborList* Neighbors;
	
	SysBox    Box;
	
	// iteration uses class these variables
	int at_idx;
	int to_idx;
	int nl_idx;

public:
	VerletListFull();
	virtual ~VerletListFull();
	
	void ReadParameter();
	bool CheckParameter();
	void GrowNeighborList(int idx);
	void Update();
	void Refresh();
	void CellNumber(int, int&, int&, int&);
	void Init(MDSystem* WorkSys);
	void CalculateBox();
	void PreIntegration();
	void IterateHalf(MDGadget* IteratedClass);
	void IterateFull(MDGadget* IteratedClass);
	void GetIterationVariables(int& at, int& to, // the index of iterating atom
							   int& at_nbidx, // the current neighborlist index
							   NeighborList*& at_nblist); // complete neighbor list
	NeighborList* GetNeighborList(){return Neighbors;}
	void Dump(string fname);
	void PrintInfo(ostream& ost);

};

}

#endif
