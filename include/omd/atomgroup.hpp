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
 *
*/

#ifndef _ATOM_GROUP_H_
#define _ATOM_GROUP_H_

#include <omd/container.hpp>

class AtomGroup: public AtomContainer {
	friend class MDSystem;
	class MDSystem* System;
	int group_flagmask;
	AtomContainer* source;
	AtomKeeper scratch_ak;
	
	// called only by MDSystem...
	void SetGroupMask(int bitpos){group_flagmask=1<<bitpos;}
	void SyncAtomGroupMask();
	
public:

	// key=tag -> a group member is tagged.
	// key=geometry -> membership is checked every step according
	//     to the defined range.
	enum GroupKey {tag, geometry} key;

	AtomGroup(string group_name, AtomContainer* ac, MDSystem* WorkSys);
	AtomGroup* Create();
	AtomGroup* Source(AtomContainer* a);
	AtomGroup* Union(AtomContainer& a);
	AtomGroup* Insert(Atom& a);
	AtomGroup* Select(int* idx_array);
	AtomGroup* SelectType(int type_id);
	AtomGroup* SelectGID(int gid);
	AtomGroup* SelectGT(OMD_FLOAT x, OMD_FLOAT y, OMD_FLOAT z); // greater than
	AtomGroup* SelectGE(OMD_FLOAT x, OMD_FLOAT y, OMD_FLOAT z); // greater or equal
	AtomGroup* SelectLT(OMD_FLOAT x, OMD_FLOAT y, OMD_FLOAT z); // less than
	AtomGroup* SelectLE(OMD_FLOAT x, OMD_FLOAT y, OMD_FLOAT z); // less or equal
	
	AtomGroup* SelectBox(OMD_FLOAT x0, OMD_FLOAT x1, 
						 OMD_FLOAT y0, OMD_FLOAT y1, 
						 OMD_FLOAT z0, OMD_FLOAT z1);
	
	AtomGroup* SelectInverseBox(OMD_FLOAT x0, OMD_FLOAT x1, 
								 OMD_FLOAT y0, OMD_FLOAT y1, 
								 OMD_FLOAT z0, OMD_FLOAT z1);
	
	AtomGroup* Commit(){return Create();}

	int  GetGroupMask() {return group_flagmask;}

	OMD_FLOAT GetMass(int idx);
};

#endif
