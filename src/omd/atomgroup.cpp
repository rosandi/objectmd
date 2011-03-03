/* omdlib
 ********************************************************
 *
 * ObjectMD
 * Molecular Dynamics Class Library
 * Version 3.0 (2009, 2011)
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
 * Atom group
 *
*/

#include <omd/atomgroup.hpp>
#include <omd/system.hpp>

AtomGroup::AtomGroup(string group_name, AtomContainer* ac, MDSystem* WorkSys) {
	set_name(group_name);
	register_class("atom_group");
	source=ac;
	System=WorkSys;
	scratch_ak.set_name("ATOM GROUP SCRATCH");
	scratch_ak.ReferArray(ac->GetAtomStorage()); // this is the scratch...
	group_flagmask=0;
}

AtomContainer* AtomGroup::Create() { // commit.. called by MDSystem...
	if(!created) group_flagmask=System->ClaimFlagBit(this, "group");
	AtomStorage.ReferArray(scratch_ak);
	created=true;
	return this;
}

AtomGroup* AtomGroup::Source(AtomContainer* ac) {
	source=ac;
	scratch_ak.ReferArray(ac->GetAtomStorage());
}

AtomGroup* AtomGroup::Union(AtomContainer& a) { // warning! combining may take time!
	AtomKeeper tmpa(AtomKeeper::Referral);
	tmpa.disable_log(LOGMEMORY);
	OMD_SIZET na=a.GetNAtom();
	OMD_SIZET nm=scratch_ak.GetNAtom();

	for(OMD_SIZET i=0;i<na;i++) {
		Atom* A=a.AtomPtr(i);
		bool doit=true;
		for(OMD_SIZET m=0;m<nm;m++) if(A==scratch_ak.AtomPtr(m)) {doit=false;break;}
		if(doit) tmpa.IndexBook.push(i);
	}

	tmpa.AssignByIndex(a.GetAtomStorage());
	scratch_ak.Append(tmpa);

	return this;
}

AtomGroup* AtomGroup::Insert(Atom& a) {
	Atom *pA=&a;
	bool doit=true;
	OMD_SIZET nm=scratch_ak.GetNAtom();

	for(OMD_SIZET m=0;m<nm;m++)
		if(pA==scratch_ak.AtomPtr(m)) {doit=false;break;}

	if(doit)scratch_ak.Attach(a);
	return this;
}

/**
 * idx_array must be terminated by a negative value
 */

AtomGroup* AtomGroup::Select(OMD_INT* idx_array) {
	AtomKeeper tmpa(AtomKeeper::Referral);
	tmpa.disable_log(LOGMEMORY);
	OMD_INT *a=idx_array;

	while(*a>=0){

		assert(*a<scratch_ak.GetNAtom(),
				"index out of range while selecting atoms to group '"+
				get_name()+"' index "+as_string(*a));

		tmpa.IndexBook.push(*a);
		a++;
	}

	tmpa.AssignByIndex(scratch_ak);
	scratch_ak.ReferArray(tmpa);
	return this;
}

AtomGroup* AtomGroup::SelectType(int type_id) {
	AtomKeeper tmpa(AtomKeeper::Referral);
	tmpa.disable_log(LOGMEMORY);
	OMD_SIZET na=scratch_ak.GetNAtom();

	for(OMD_SIZET i=0;i<na;i++) {
		if(scratch_ak[i].id==type_id) tmpa.IndexBook.push(i);
	}

	tmpa.AssignByIndex(scratch_ak);
	scratch_ak.ReferArray(tmpa);
	return this;
}

AtomGroup* AtomGroup::SelectXID(int xid) {
	AtomKeeper tmpa(AtomKeeper::Referral);
	tmpa.disable_log(LOGMEMORY);
	OMD_SIZET na=scratch_ak.GetNAtom();

	for(OMD_SIZET i=0;i<na;i++) {
		if(scratch_ak[i].xid==xid) tmpa.IndexBook.push(i);
	}

	tmpa.AssignByIndex(scratch_ak);
	scratch_ak.ReferArray(tmpa);
	return this;
}

AtomGroup* AtomGroup::SelectLT(OMD_FLOAT x, OMD_FLOAT y, OMD_FLOAT z) {
	AtomKeeper tmpa(AtomKeeper::Referral);
	tmpa.disable_log(LOGMEMORY);
	OMD_SIZET na=scratch_ak.GetNAtom();

	for(OMD_SIZET i=0;i<na;i++) {
		if((scratch_ak[i].x<x)&&(scratch_ak[i].y<y)&&(scratch_ak[i].z<z)) tmpa.IndexBook.push(i);
	}

	tmpa.AssignByIndex(scratch_ak);
	scratch_ak.ReferArray(tmpa);

	return this;
}

AtomGroup* AtomGroup::SelectGT(OMD_FLOAT x, OMD_FLOAT y, OMD_FLOAT z)  {
	AtomKeeper tmpa(AtomKeeper::Referral);
	tmpa.disable_log(LOGMEMORY);
	OMD_SIZET na=scratch_ak.GetNAtom();
	for(OMD_SIZET i=0;i<na;i++) {
		if((scratch_ak[i].x>x)&&(scratch_ak[i].y>y)&&(scratch_ak[i].z>z)) tmpa.IndexBook.push(i);
	}
	tmpa.AssignByIndex(scratch_ak);
	scratch_ak.ReferArray(tmpa);
	return this;
}

AtomGroup* AtomGroup::SelectGE(OMD_FLOAT x, OMD_FLOAT y, OMD_FLOAT z) {
	AtomKeeper tmpa(AtomKeeper::Referral);
	tmpa.disable_log(LOGMEMORY);
	OMD_SIZET na=scratch_ak.GetNAtom();

	for(OMD_SIZET i=0;i<na;i++) {
		if((scratch_ak[i].x>=x)&&(scratch_ak[i].y>=y)&&(scratch_ak[i].z>=z)) tmpa.IndexBook.push(i);
	}

	tmpa.AssignByIndex(scratch_ak);
	scratch_ak.ReferArray(tmpa);
	return this;
}

AtomGroup* AtomGroup::SelectLE(OMD_FLOAT x, OMD_FLOAT y, OMD_FLOAT z) {
	AtomKeeper tmpa(AtomKeeper::Referral);
	tmpa.disable_log(LOGMEMORY);
	OMD_SIZET na=scratch_ak.GetNAtom();

	for(OMD_SIZET i=0;i<na;i++) {
		if((scratch_ak[i].x<=x)&&(scratch_ak[i].y<=y)&&(scratch_ak[i].z<=z)) tmpa.IndexBook.push(i);
	}

	tmpa.AssignByIndex(scratch_ak);
	scratch_ak.ReferArray(tmpa);
	return this;
}

OMD_FLOAT AtomGroup::GetMass(OMD_SIZET idx) {
	return System->GetMass(idx);
}
