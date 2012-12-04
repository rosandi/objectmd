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
 * Tabular EAM potential:
 * 
*/

#include <fstream>
#include <cstdio>
#include <cstring>
#include <omd/system.h>
#include <omd/iterator.h>
#include <potential/team.h>

using namespace omd;

//--Electron-Density-------------------------------

// AtomID Matrix format:
// Row/Col Real ID 
// Value = -1 --> not belongs to EAM
// otherwise = table index

#define CheckID(A) (AtomID[A][A]>=0)

TEmbedding::TEmbedding() {
	id_size=0;
	Active=ActiveCode=MODIFY_PRE_CALCULATION|MODIFY_POST_FORCE;
	register_class("embedding_function");
	set_name("TABLE_EMBEDDING_FUNCTION");
	for (int i=0; i<MAX_ALLOWED_SPECIES; i++) 
		for (int j=0; j<MAX_ALLOWED_SPECIES; j++) AtomID[i][j]=-1;
}

// for the same atom type 
void TEmbedding::AddTable(int aid, string tabfile) {
	mdassert(aid<MAX_ALLOWED_SPECIES, "attempt to insert table for atom.tid="+as_string(aid)+
	       " "+as_string(MAX_ALLOWED_SPECIES)+" maximum allowed atom species exeeded!");
	AtomID[aid][aid]=id_size;
	embed[id_size  ].open(tabfile,"EMBED");
	edens[id_size++].open(tabfile,"EDENS");
}

// for different atom type
void TEmbedding::AddTable(int aid_at, int aid_to, string tabfile) {
	mdassert(aid_at<MAX_ALLOWED_SPECIES&&aid_to<MAX_ALLOWED_SPECIES,
	       "attempt to insert table for atom.tid=("+as_string(aid_at)+","+as_string(aid_to)+
	       ") "+as_string(MAX_ALLOWED_SPECIES)+" maximum allowed atom species exeeded!");

	if (aid_at==aid_to) AddTable(aid_at, tabfile);
	else {
		AtomID[aid_at][aid_to] = AtomID[aid_to][aid_at] = id_size;
		edens[id_size++].open(tabfile,"EDENS");
	}
}

void TEmbedding::Init(MDSystem* WorkSys) {
	Modify::Init(WorkSys);	
	
	mdassert(MAX_ALLOWED_SPECIES>=WorkSys->SystemAtoms.size(),
	       "atom type limit "+as_string(MAX_ALLOWED_SPECIES)+" exceeded");

	// Claim one auxilary variable
	ED_idx=ClaimAuxVariable();

	mdassert(id_size>0, "uninitialized electron density table");	
	for (int i=0; i<id_size; i++) CutRadius[i]=edens[i].max_range();
}

void TEmbedding::IterationNode(Atom &at, Atom &to) {	
	OMD_FLOAT R, Dens, Dx, Dy, Dz;
	int   Tab;	
	Tab=AtomID[(int)(at.tid)][(int)(to.tid)]; // Take edens table index
	if (Tab<0) return; // avoids non eam interactions
	R=sqrt(CalcSqrDistance(at, to, Dx, Dy, Dz));
	if (R<=CutRadius[Tab]) {
		Dens=edens[Tab].read(R);
		at.aux[ED_idx] += Dens;
		to.aux[ED_idx] += Dens;
	}
}

/**
 * Refresh and update the atom density.
 * Atom data is re-synchronize after updating atom density
 */

void TEmbedding::PreCalculation() {
	int na=GetNAtom();
	for (int i = 0; i<na; i++) Atoms(i).aux[ED_idx]=0.0;
	Iterator->IterateHalf(this);
	SyncData(SYNC_AUX,ED_idx);
}

void TEmbedding::Dump() {
	ofstream fl("dump.rho");
	for (int i = 0; i < GetNAtom(); i++) 
		fl << Atoms(i).x << " " 
		   << Atoms(i).y << " " 
		   << Atoms(i).z << " "
		   << Atoms(i).aux[ED_idx]<< "\n";
	fl.close();
}

OMD_FLOAT TEmbedding::GetEmb(Atom &at) {
	int tid=at.tid;
	return embed[AtomID[tid][tid]].read(at.aux[ED_idx]);
}

OMD_FLOAT TEmbedding::GetEmbDeriv(Atom &at) {
	int tid=at.tid;
	return embed[AtomID[tid][tid]].dread(at.aux[ED_idx]);
}

OMD_FLOAT TEmbedding::GetRho(Atom &at) {return at.aux[ED_idx];}

OMD_FLOAT TEmbedding::GetRhoDeriv(OMD_FLOAT r, Atom &at) {
	int tid=at.tid;
	if (r>CutRadius[AtomID[tid][tid]]) return 0.0;
	return edens[AtomID[tid][tid]].dread(r);
}

// The correction by electron density is done here.. 
void TEmbedding::PostForce() {
	int na=GetNAtom();
	for(int i=0;i<na;i++) {
		Atom* a=AtomPtr(i);
		if(a->flag&FLAG_GHOST)continue;
		if(a->flag&FLAG_ACTIVE){
			int id=a->tid;
			if(0<=AtomID[id][id]) a->potential+=GetEmb(*a);
		}
	}
}

void TEmbedding::PrintInfo(ostream& ost) {
	ost<<"id."<<id<<" "<<get_name()<< "\n";
	ost<<" electron_density:\n";
	for(int i=0;i<id_size;i++)
		for(int j=i;j<id_size;j++)
		if(AtomID[i][j]>=0) 
			ost <<"  "<<System->SystemAtoms[i]->get_name()<<"."<<i<<"<-->"
			    <<System->SystemAtoms[j]->get_name()<<"."<<j<<" "
		        <<edens[AtomID[i][j]].filename<< "\n";
	ost <<" embedding function:\n";
	for(int i=0;i<id_size;i++)
		if(AtomID[i][i]>=0) 
			ost <<"  "<<System->SystemAtoms[i]->get_name()<<"."<<i
		        <<" "<<embed[AtomID[i][i]].filename << "\n";
	ost << "\n";
}

//--Force-EAM----------------

TForceEAM::TForceEAM(string PhiFile, TEmbedding *EM) {
	emb=EM;
	tablefile.assign(PhiFile);	
	set_name("TABLE_FORCE_EAM");
	register_class(get_name());
}

TForceEAM* TForceEAM::SetTable(const char* PhiFile) {
	tablefile.assign(PhiFile);
	return this;
}

bool TForceEAM::SearchAttachEmbeddingClass() {

	if(emb){
		blog("already inserted: "+emb->get_name()+"."+as_string(emb->get_id())+
		    "embedding function handler for "+as_string(AtomTypeA)+"<-->"+
			 as_string(AtomTypeB), LOGCREATE);
		return true;
	}

	// find first the Embedding energy handler in modify list
	for (int i=0; i<(int)System->Modifies.size(); i++) {
			
		if(System->Modifies[i]->type_of("EMBEDDING_FUNCTION")){
			emb=dynamic_cast<TEmbedding*>(System->Modifies[i]);
			emb->AddTable(AtomTypeA, AtomTypeB, phi.filename);
			blog("using '"+emb->get_name()+"' "+
			    "embedding function for "+as_string(AtomTypeA)+
				 "<-->"+as_string(AtomTypeB), LOGCREATE);
			return true;
		}
	}

	return false;
}

// ForceKernel is initialized before Modifies. 
// Modify insertion is save here.

void TForceEAM::Init(MDSystem* WorkSys) {
	ForceKernel::Init(WorkSys);
	// saver: let Worksys do the work...
	tablefile=WorkSys->search_path("$OMD_TABLE", "eam."+tablefile);

	phi.open(tablefile, "PAIR");
	CutRadius=phi.max_range();
	CutRadiusSqr=CutRadius*CutRadius;
	
	if(!SearchAttachEmbeddingClass()) {
		// If no EmbeddingEnergyHandler is found... create one.
		emb = new TEmbedding();
		emb->set_logger(logger);
		emb->AddTable(AtomTypeA, AtomTypeB, phi.filename);
		WorkSys->AddModify(emb);
	}
}

void TForceEAM::ComputeHalf(Atom &at, Atom &to) {
	OMD_FLOAT EmbAt, EmbTo, DRhoAt, DRhoTo;
	OMD_FLOAT fr, R, dx, dy, dz, p, dp;
	R = sqrt(CalcSqrDistance(at, to, dx, dy, dz));
	if (R<=CutRadius) {
		phi.read(R, p, dp);
		EmbAt     = emb->GetEmbDeriv(at);
		EmbTo     = emb->GetEmbDeriv(to);
		DRhoAt    = emb->GetRhoDeriv(R, at);
		DRhoTo    = emb->GetRhoDeriv(R, to);
		fr=-(EmbAt*DRhoTo+EmbTo*DRhoAt+dp)/R;
		ReturnForce(at,to,dx,dy,dz,fr,p);
	}
}

void TForceEAM::PrintInfo(ostream& ost){
	ost  <<"id."<<id<<" "<<get_name()<<"\n"
	     <<" atoms " << System->SystemAtoms[AtomTypeA]->get_name()
	     <<"<->"<< System->SystemAtoms[AtomTypeB]->get_name()<< "\n"
	     <<" table="<<phi.filename<<" cut_radius="<<CutRadius << "\n"
	     <<" embedding class id="<<emb->get_id()<<" "<< emb->get_name()<<"\n";
}
