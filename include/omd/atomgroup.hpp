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
	class MDSystem* System;
	OMD_SIZET group_flagmask;
	AtomContainer* source;
	AtomKeeper scratch_ak;
public:
	AtomGroup(string group_name, AtomContainer* ac, MDSystem* WorkSys);
	AtomContainer* Create();
	AtomGroup* Source(AtomContainer* a);
	AtomGroup* Union(AtomContainer& a);
	AtomGroup* Insert(Atom& a);
	AtomGroup* Select(int* idx_array);
	AtomGroup* SelectType(int type_id);
	AtomGroup* SelectXID(int xid);
	AtomGroup* SelectGT(OMD_FLOAT x, OMD_FLOAT y, OMD_FLOAT z); // greater than
	AtomGroup* SelectGE(OMD_FLOAT x, OMD_FLOAT y, OMD_FLOAT z); // greater or equal
	AtomGroup* SelectLT(OMD_FLOAT x, OMD_FLOAT y, OMD_FLOAT z); // less than
	AtomGroup* SelectLE(OMD_FLOAT x, OMD_FLOAT y, OMD_FLOAT z); // less or equal
	AtomGroup* Commit(){return dynamic_cast<AtomGroup*>(Create());}

	OMD_SIZET GetGroupFlagMask() {return group_flagmask;}
	OMD_FLOAT GetMass(OMD_SIZET idx);
};

#endif
