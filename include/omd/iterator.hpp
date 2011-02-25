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
	bool dirty;
	
public:
	
	MDIterator(){
		set_name("iterator");
		register_class(get_name());
		Active=0;
		ActiveCode=-1;
		RadiusTolerance=0.0;
		dirty=true;
	}
	
	virtual ~MDIterator(){}
	
	/** Update() is called every time before iteration loop */
	virtual void Update(){}
	
	/** Refresh() is preserved for recreating/recalculating the neighbor cell */
	virtual void Refresh(){}
	
	virtual void SetRadiusTolerance(OMD_FLOAT tole){RadiusTolerance=tole;}
	virtual void SetDirty(){dirty=true;}

	/**
	 * Perform iteration on all atoms in the system. When implementing an
	 * Iterrator class, the descendant must ensure to handle the inactive and
	 * ghost atoms properly. The corresponding status can be checked by
	 * CheckActive and CheckGhost() class member of MGGadget.
	 * 
	 */

	virtual void Iterate(MDGadget* IteratedClass, bool force_update=false) {
		if(dirty||force_update) Update();
		OMD_SIZET na=GetNAtom();
		for(OMD_SIZET at=0;at<na-1;at++) {
			if(!IteratedClass->PreIterationNode(at))continue;
			for(OMD_SIZET to=(at+1);to<na;to++) {
				IteratedClass->IterationNode(at,to);
			}
		}
	}
};

#endif
