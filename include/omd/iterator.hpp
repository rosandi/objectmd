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
 * Iterator class -> keep the iteration method through
 * all particles in the simulation system.
 *
*/

#ifndef _ITERATOR_H_
#define _ITERATOR_H_

#include <omd/conditioner.hpp>

//-------------------Iterator---------------------//

/**
 * @ingroup iterator
 * @brief Do iteration loops
 * 
 * The iterator calls the IterationNode() function of the
 * MDGadget class. 
 * 
 */

class MDIterator:public Conditioner {
protected:
	OMD_FLOAT RadiusTolerance;
	int   RebuildPeriod;
	int   UpdatePeriod;
	bool  dirty;
	
public:
	bool  full_loop;
	bool  half_loop;
	bool  looping_full;
	
	MDIterator(){
		set_name("iterator");
		register_class(get_name());
		Active=0;
		ActiveCode=-1;
		RadiusTolerance=0.0;
		RebuildPeriod=0;
		UpdatePeriod=0;
		dirty=true;
		full_loop=false;
		half_loop=false;
	}
	
	virtual ~MDIterator(){}
	
	/** 
	 Called before the iteration if the iterator is dirty, set by SetDirty(), 
	 or forced implicitly by the argument of Iterate() function.
	 */
	
	virtual void Update(){}
	
	/** 
	 Preserved for recreating/recalculating the neighbor cell, for instance
	 upon box size change, etc. This function is activated by activating the
	 COND_PRE_CALCULATION conditioner type, SetConditionerType().
	 */
	
	virtual void Refresh(){}
	
	virtual void SetRadiusTolerance(OMD_FLOAT tole){
		RadiusTolerance=tole;
	}
	
	void EnableFullLoop(){full_loop=true;}
	void EnableHalfLoop(){half_loop=true;}
	void DisableFullLoop(){full_loop=false;}
	void DisableHalfLoop(){half_loop=false;}

	void SetDirty(){dirty=true;}
	void SetUpdatePeriod(int up){UpdatePeriod=up;}
	void SetRebuildPeriod(int bp){RebuildPeriod=bp;}
	
	void PreCalculation() {
		if(RebuildPeriod) {
			if(NCalls)
				if(!(System->Step%RebuildPeriod)) Refresh();
		}
		if(UpdatePeriod) {
			if(NCalls)
				if(!(System->Step%UpdatePeriod)) SetDirty();
		}
	}

	/**
	 * Perform half iteration on all atoms in the system. When implementing an
	 * Iterrator class, the descendant must ensure to handle the inactive and
	 * ghost atoms properly. The corresponding status can be checked by
	 * CheckActive and CheckGhost() class member of MGGadget.
	 * 
	 */

	virtual void IterateHalf(MDGadget* IteratedClass) {
		if(dirty) Update();
		looping_full=false; // FALSE -> half loop
		
		int na=GetNAtom();
		for(int at=0;at<na-1;at++) {
			if(!IteratedClass->PreIterationNode(at))continue;
			for(int to=(at+1);to<na;to++) {
				IteratedClass->IterationNode(at,to);
			}
		}
	}
	
	/**
	 * Perform full iteration on all atoms in the system.
	 * 
	 */
	
	virtual void IterateFull(MDGadget* IteratedClass) {
		if(dirty) Update();
		looping_full=true; // doing full loop
		
		int na=GetNAtom();
		for(int at=0;at<na;at++) {
			if(!IteratedClass->PreIterationNode(at))continue;
			for(int to=0;to<na;to++) {
				if(at==to) continue;
				IteratedClass->IterationNode(at,to);
			}
		}
	}
	
	virtual void Iterate(MDGadget* IteratedClass) {
		if(half_loop) IterateHalf(IteratedClass);
		if(full_loop) IterateFull(IteratedClass);
	}
	
};

#endif
