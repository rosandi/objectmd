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
#include <omd/base.hpp>
#include <omd/gadget.hpp>
#include <omd/unit.hpp>
#include <omd/detector.hpp>
#include <omd/container.hpp>
#include <omd/integrator.hpp>
#include <omd/system.hpp>

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

OMD_INT MDGadget::GetContainerID(string name){
	return System->GetContainerID(name);
}

AtomContainer* MDGadget::SearchTarget(string name, bool strict){
	if(name==System->get_name()) return System;
	for(OMD_SIZET i=0;i<System->SystemAtoms.size();i++)
		if(name==System->SystemAtoms[i]->get_name())
			return System->SystemAtoms[i];
	if(strict) die("can not find target atom '"+name+"'");
	return NULL; // avoids warning..
}

MDGadget* MDGadget::SearchGadget(string name, bool strict) {
	for(OMD_SIZET i=0;i<System->Detectors.size();i++) {
		if(System->Detectors[i]->get_name()==name)
			return System->Detectors[i];
	}
	
	for(OMD_SIZET i=0;i<System->Conditioners.size();i++) {
		if(System->Conditioners[i]->get_name()==name)
			return System->Conditioners[i];
	}
	
	if(strict) die("can not find gadget '"+name+"'");
	return NULL;
}

// case insensitive
MDGadget* MDGadget::SearchGadgetType(string type, bool strict) {
	type=(replace_char(lower_case(type), ' ', '_'));
	for(OMD_SIZET i=0;i<System->Detectors.size();i++)
		if(System->Detectors[i]->type_of(type))
			return System->Detectors[i];

	for(OMD_SIZET i=0;i<System->Conditioners.size();i++)
		if(System->Conditioners[i]->type_of(type))
			return System->Conditioners[i];

	if(strict) die("can not find gadget of type '"+type+"'");
	return NULL;
}

void MDGadget::Init(MDSystem* WorkSys){
	if(!System)System=WorkSys;
	if(!Unit)Unit=System->Unit;
	SysParam=&(System->param);
	Iterator=WorkSys->GetIterator();
	Integrator=WorkSys->GetIntegrator();
	if(Target==NULL){
		if(TargetName=="-"||TargetName=="")Target=WorkSys;
		else Target=SearchTarget(TargetName);
	}

	IsReady  = true;	
}

// redirect atom access to Target...    
Atom&  MDGadget::Atoms(OMD_INT idx) {return Target->Atoms(idx);}
Atom*  MDGadget::AtomPtr(OMD_INT idx) {return Target->AtomPtr(idx);}
OMD_SIZET   MDGadget::GetNAtom() {return Target->GetNAtom();}
OMD_FLOAT MDGadget::GetTimeStep() {return System->Integrator->TimeStep;}
OMD_FLOAT MDGadget::GetElapsedTime() {return System->ElapsedTime;}
OMD_FLOAT MDGadget::GetMass(Atom &a) {return System->GetMass(a.id);}
OMD_FLOAT MDGadget::GetMass(Atom *a) {return System->GetMass(a->id);}
OMD_FLOAT MDGadget::GetMass(OMD_INT atomid) {return System->GetMass(atomid);}
OMD_FLOAT MDGadget::GetNumber(Atom &a) {return System->GetNumber(a.id);}
OMD_FLOAT MDGadget::GetNumber(OMD_INT atomid) {return System->GetNumber(atomid);}
OMD_SIZET   MDGadget::ClaimFlagBit() {return System->ClaimFlagBit(this);}
bool   MDGadget::OnTime(OMD_FLOAT tm){return System->OnTime(tm);}
OMD_FLOAT&MDGadget::AuxVariable(OMD_SIZET i){return Atoms(i).aux[AuxIdx];}
OMD_SIZET   MDGadget::ClaimAuxVariable(bool printable, const OMD_CHAR* tag, const OMD_CHAR* sformat) {
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
	
	// FIXME! proximity limit....
	
	if(r2<1e-6){
		die("overlapping atoms: r="+as_string(sqrt(r2),"%1.1E")+
			" between atoms ("+
			as_string(at.x)+","+
			as_string(at.y)+","+
			as_string(at.z)+") and ("+
			as_string(to.x)+","+
			as_string(to.y)+","+
			as_string(to.z)+")");
	}

	return r2;
}

OMD_INT MDGadget::IsActive(OMD_INT Code) {
	if(Code) return Code&Active;
	return ActiveCode&Active;
}

DataSlot* MDGadget::RegisterMessageSlot(DataSlot* slot) {
	assert(System, "attempt to register slot before initialization");
	System->MessageSlots.push_back(slot);
	return slot;
}


void MDGadget::RestartVariable(string tag, OMD_FLOAT &val){
	assert(System, "attempt to read restart variable before initialization");
	DataSlot *m=NULL;

	for(OMD_SIZET i=0;i<System->RestartVars.size();i++) {
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

void MDGadget::RestartVariable(string tag, OMD_INT &val){
	assert(System, "attempt to read restart variable befor initialization");
	DataSlot *m=NULL;

	for(OMD_SIZET i=0;i<System->RestartVars.size();i++) {
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
