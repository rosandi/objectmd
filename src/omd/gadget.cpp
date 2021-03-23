/* omdlib
 ********************************************************
 *
 * ObjectMD
 * Molecular Dynamics Class Library
 * Version 2.0 (2009)
 *
 * Released under:
 *       GNU GENERAL PUBLIC LICENSE Version 3
 *       Copyright (C) 2005  Yudi Rosandi
 *
 * Read license.txt file in the source's root directory.
 *
 ********************************************************
 *
 * Implementation of Crystal class
 *
*/

#include <cmath>
#include <omd/omdtool.h>
#include <omd/base.h>
#include <omd/gadget.h>
#include <omd/unit.h>
#include <omd/detector.h>
#include <omd/container.h>
#include <omd/integrator.h>
#include <omd/system.h>

using namespace omd;

MDGadget::MDGadget() {
	System= NULL;
	Target=NULL;
	Iterator=NULL;
	Integrator=NULL;
	IsReady = false;
	Active=ActiveCode=1;
	AuxIdx=-1;
	NCalls=0;
	TargetName="";
	Unit=NULL;
	set_name("GADGET");
	register_class("GADGET");
}

int MDGadget::GetContainerID(string name){
	return System->GetContainerID(name);
}

AtomContainer* MDGadget::SearchTarget(string name, bool strict){
	if(name==System->get_name()) return System;

	// Atom group has priority...
	for(int i=0;i<(int)System->SystemAtomGroups.size();i++)
		if(name==System->SystemAtomGroups[i]->get_name())
			return System->SystemAtomGroups[i];

	for(int i=0;i<(int)System->SystemAtoms.size();i++)
		if(name==System->SystemAtoms[i]->get_name())
			return System->SystemAtoms[i];

	if(strict) die("can not find target atom '"+name+"'");
	return NULL; // avoids warning..
}

MDGadget* MDGadget::SearchGadget(string name, bool strict) {
	for(int i=0;i<(int)System->Detectors.size();i++) {
		if(System->Detectors[i]->get_name()==name)
			return System->Detectors[i];
	}
	
	for(int i=0;i<(int)System->Modifies.size();i++) {
		if(System->Modifies[i]->get_name()==name)
			return System->Modifies[i];
	}
	
	if(strict) die("can not find gadget '"+name+"'");
	return NULL;
}

// case insensitive
MDGadget* MDGadget::SearchGadgetType(string type, bool strict) {
	type=(replace_char(lower_case(type), ' ', '_'));
	for(int i=0;i<(int)System->Detectors.size();i++)
		if(System->Detectors[i]->type_of(type))
			return System->Detectors[i];

	for(int i=0;i<(int)System->Modifies.size();i++)
		if(System->Modifies[i]->type_of(type))
			return System->Modifies[i];

	if(strict) die("can not find gadget of type '"+type+"'");
	return NULL;
}

void MDGadget::Init(MDSystem* WorkSys){
	if(!System)System=WorkSys;
	if(!Unit)Unit=System->Unit;
	SysParam=&(System->param);
	if (!Iterator) Iterator=WorkSys->GetIterator();
	if (!Integrator) Integrator=WorkSys->GetIntegrator();
	ReadParameter();
	mdassert(CheckParameter(), "check parameter failed");
	
	if(Target==NULL){
		string nst=lower_case(TargetName);
		if(nst=="-"||nst==""||nst=="all"||nst=="system")Target=WorkSys;
		else Target=SearchTarget(TargetName);
	}
	
	IsReady  = true;
}

// redirect atom access to Target...    
Atom&  MDGadget::Atoms(int idx) {return Target->Atoms(idx);}
Atom*  MDGadget::AtomPtr(int idx) {return Target->AtomPtr(idx);}

int   MDGadget::GetNAtom() {return Target->GetNAtom();}
OMD_FLOAT MDGadget::GetTimeStep() {return System->Integrator->TimeStep;}
OMD_FLOAT MDGadget::GetElapsedTime() {return System->ElapsedTime;}
OMD_FLOAT MDGadget::GetMass(Atom &a) {return System->SystemAtoms[a.tid]->M;}
OMD_FLOAT MDGadget::GetMass(Atom *a) {return System->SystemAtoms[a->tid]->M;}
OMD_FLOAT MDGadget::GetMass(int idx) {return System->GetMass(idx);}

OMD_FLOAT MDGadget::GetZ(Atom &a) {return System->SystemAtoms[a.tid]->Z;}
OMD_FLOAT MDGadget::GetZ(int idx) {return System->GetZ(idx);}

int   MDGadget::ClaimFlagBit() {return System->ClaimFlagBit(this);}
bool  MDGadget::OnTime(OMD_FLOAT tm){return System->OnTime(tm);}
bool  MDGadget::OnStep(int step){return System->OnStep(step);}
OMD_FLOAT&MDGadget::AuxVariable(int i){return Atoms(i).aux[AuxIdx];}
int   MDGadget::ClaimAuxVariable(bool printable, const char* tag, const char* sformat) {
	AuxIdx=System->ClaimAuxVariable(this,printable,tag,sformat);
	return (AuxIdx);
}

OMD_FLOAT MDGadget::CalcSqrDistance(Atom &at, Atom &to, OMD_FLOAT &dx, OMD_FLOAT &dy, OMD_FLOAT &dz, bool check) {
	
	dx=at.x-to.x;
	dy=at.y-to.y;
	dz=at.z-to.z;

	System->BoundaryCorrectDistances(dx, dy, dz);
	OMD_FLOAT r2=(dx*dx+dy*dy+dz*dz);
	
	if(!check) return r2;
	
	if(r2<OMD_PROXIMITY){
		die("CalcSqrDistance: overlapping atoms: r="+as_string(sqrt(r2),"%1.1E")+
			" between atom nid="+
			as_string(at.nid)+" and nid="+as_string(to.nid)+" ("+
			as_string(at.x)+","+
			as_string(at.y)+","+
			as_string(at.z)+") -- ("+
			as_string(to.x)+","+
			as_string(to.y)+","+
			as_string(to.z)+")");
	}

	return r2;
}

int MDGadget::IsActive(int Code) {
	if(Code) return Code&Active;
	return ActiveCode&Active;
}

DataSlot* MDGadget::RegisterMessageSlot(DataSlot* slot) {
	mdassert(System, "attempt to register slot before initialization");
	System->MessageSlots.push_back(slot);
	return slot;
}


void MDGadget::RestartVariable(string tag, OMD_FLOAT &val){
	mdassert(System, "attempt to read restart variable before initialization");
	DataSlot *m=NULL;

	for(int i=0;i<(int)System->RestartVars.size();i++) {
		if(System->RestartVars[i]->GetLabel()==tag) m=System->RestartVars[i];
	}
	
	if(!m){
		m=new DataSlot(tag);
		System->RestartVars.push_back(m);
	} else {
		val=m->AsDouble();
	}

	m->SetData(val);

}

void MDGadget::RestartVariable(string tag, int &val){
	mdassert(System, "attempt to read restart variable befor initialization");
	DataSlot *m=NULL;

	for(int i=0;i<(int)System->RestartVars.size();i++) {
		if(System->RestartVars[i]->GetLabel()==tag) m=System->RestartVars[i];
	}
	
	if(!m){
		m=new DataSlot(tag);
		System->RestartVars.push_back(m);
	} else {
		val=m->AsInt();
	}

	m->SetData(val);
}
