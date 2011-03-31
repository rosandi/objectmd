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
 *
*/

#include <conditioner/VerletListFull.hpp>

using std::cerr;
using std::ofstream;

// FIXME! make a note on VerletRadius and Radius tolerance. CommHandler alter these
// variable. What is the effect? ...

VerletListFull::VerletListFull() {
	set_name("VERLET LIST FULL NEIGHBOR");
	register_class(get_name());
	SetConditionerType(COND_PRE_INTEGRATION|COND_PRE_CALCULATION);
	AllocSize=0;
	Link=NULL;
	Neighbors=NULL;
	half_loop=true;
	full_loop=true;
}

VerletListFull::~VerletListFull() {
	MemFree(Link);
	for(int i=0;i<AllocSize;i++) {
		MemFree(Neighbors[i].list);
		MemFree(Neighbors[i].flag);		
	}
	MemFree(Neighbors);
}

void VerletListFull::ReadParameter() {
	SysParam->peek("verlet.update", UpdatePeriod, 5);
	SysParam->peek("verlet.rebuild", RebuildPeriod, 0);
	SysParam->peek("verlet.radius", CutRadius, -1.0);
	SysParam->peek("verlet.radtole", RadiusTolerance, -1.0);
	SysParam->peek("verlet.n_mean", nmean, 120);
}

bool VerletListFull::CheckParameter() {
	if(RadiusTolerance<=0.0) {
		warn("no radius tolerance for verlet list (verlet.radtole). Setting update period to 1");
		RadiusTolerance=0.0;
		UpdatePeriod=1.0;
	}
	return true;
}

void VerletListFull::Init(MDSystem* WorkSys) {
	Conditioner::Init(WorkSys);
	
	if (System->PBoundary) {
		assert(System->type_of("simulation_system_grid"), 
			   "MDSystemGrid is required for periodic boundary simulation");
	}
	
	if(CutRadius<=0.0) CutRadius=System->GetMaxCutRadius();
	VerletRadius=CutRadius+RadiusTolerance;

}

void VerletListFull::CellNumber(int at, int& xid, int& yid, int& zid) {
	xid=(int)((Atoms(at).x-Box.x0)/VerletRadius);
	yid=(int)((Atoms(at).y-Box.y0)/VerletRadius);
	zid=(int)((Atoms(at).z-Box.z0)/VerletRadius);
	if(xid>(U-1))xid=U-1;
	if(xid<0)xid=0;
	if(yid>(V-1))yid=V-1;
	if(yid<0)yid=0;
	if(zid>(W-1))zid=W-1;
	if(zid<0)zid=0;
}

void VerletListFull::CalculateBox() {
	OMD_FLOAT MinX, MaxX, MinY, MaxY, MinZ, MaxZ;        
	int na=GetNAtom();
	MinX=MinY=MinZ= DBL_MAX;
	MaxX=MaxY=MaxZ=-DBL_MAX;

	VerletRadius=CutRadius+RadiusTolerance;

	for(int i=0; i<na; i++) {
		if (!CheckActive(i)) continue; // consider only active atoms...
		if (MinX>Atoms(i).x)MinX=Atoms(i).x;
		if (MinY>Atoms(i).y)MinY=Atoms(i).y;
		if (MinZ>Atoms(i).z)MinZ=Atoms(i).z;
		if (MaxX<Atoms(i).x)MaxX=Atoms(i).x;
		if (MaxY<Atoms(i).y)MaxY=Atoms(i).y;
		if (MaxZ<Atoms(i).z)MaxZ=Atoms(i).z;
	}

	Box.x0=MinX;Box.x1=MaxX;
	Box.y0=MinY;Box.y1=MaxY;
	Box.z0=MinZ;Box.z1=MaxZ;	
	Box.lx=fabs(Box.x1-Box.x0);
	Box.ly=fabs(Box.y1-Box.y0);
	Box.lz=fabs(Box.z1-Box.z0);
	Box.hlx=Box.lx/2.0;
	Box.hly=Box.ly/2.0;
	Box.hlz=Box.lz/2.0;
	U = (int)(Box.lx/VerletRadius); if(U==0) U=1;
	V = (int)(Box.ly/VerletRadius); if(V==0) V=1;
	W = (int)(Box.lz/VerletRadius); if(W==0) W=1;
	blog("linked cell dimension ("+as_string(U)+","+as_string(V)+","+as_string(W)+")");
}

void VerletListFull::Refresh(){
	CalculateBox();
	SetDirty();
}

void VerletListFull::GrowNeighborList(int idx) {
	int newsize=Neighbors[idx].size+nmean;
	MemRealloc(Neighbors[idx].list, newsize*sizeof(int));
	MemRealloc(Neighbors[idx].flag, newsize*sizeof(int));
	for(int i=Neighbors[idx].size;i<newsize;i++) {
		Neighbors[idx].list[i]=-1;
		Neighbors[idx].flag[i]=-1;
	}
	Neighbors[idx].size=newsize;
}

void VerletListFull::Update() {

	int na=GetNAtom();
	
	if(AllocSize<na) {
		MemRealloc(Link, na*sizeof(int));
		MemRealloc(Neighbors, na*sizeof(NeighborList));
		for(int i=AllocSize;i<na;i++) {
			Neighbors[i].size=0.0;
			MemAlloc(Neighbors[i].list,nmean*sizeof(int));
			MemAlloc(Neighbors[i].flag,nmean*sizeof(int));
			memset(Neighbors[i].list,-1,nmean*sizeof(int));
			memset(Neighbors[i].flag,-1,nmean*sizeof(int));			
		}
		AllocSize=na;
	}

	int CellHeaders[U][V][W];

	// Initiate the Link cells and header tabel	
	for (int at=0;at<na;at++) Link[at]=-1;

	for (int i=0;i<U;i++)
		for(int j=0;j<V;j++)
			for(int k=0;k<W;k++)
				CellHeaders[i][j][k]=-1;
	
	// Find membership of every atom in the system
	// Remember! smallest number in cell must be the header
	// This algorithm ensures that.
	// Link[i] --> next index after i
	
	for (int at=na-1;at>=0;at--) {
		int u, v, w;
		CellNumber(at, u, v, w);
		Link[at]=CellHeaders[u][v][w];
		CellHeaders[u][v][w]=at;
	}

	// Create Verlet neighbor list
	// Every atom only search to atoms in it cell and 
	// other 26 neighboring cells
	OMD_FLOAT  xto, yto, zto;
	OMD_FLOAT  rd, Rsq=VerletRadius*VerletRadius;
//	int NCounter=0;

	for(int i=0;i<na;i++) Neighbors[i].end=0;
	
	for (int me=0;me<na;me++) {
		int u, v, w;
		CellNumber(me, u, v, w);
		
		if(me!=CellHeaders[u][v][w])
			die("wrong header at index "+as_string(me)+
			    " cell index ("+as_string(u)+","+as_string(v)+","+as_string(w)+")");

		//start=NeighborIndex[ni]; end=NeighborIndex[ni+1];
		CellHeaders[u][v][w]=Link[me];
		Neighbors[me].start=Neighbors[me].end;

		for (int k=w-1; k<=w+1; k++) {
			if(k==-1||k==W) continue;
			for (int j=v-1; j<=v+1; j++) {
				if(j==-1||j==V) continue;
				for (int i=u-1; i<=u+1; i++) {
					if(i==-1||i==U) continue;
					int to=CellHeaders[i][j][k];
					while (to>=0) {
						if(me==to) continue;
						rd=CalcSqrDistance(me, to, xto, yto, zto);
						if(rd<=Rsq) {
							
							if(Neighbors[me].end>=Neighbors[me].size)
								GrowNeighborList(me);
							if(Neighbors[to].end>=Neighbors[to].size)
								GrowNeighborList(to);
							
							Neighbors[me].list[Neighbors[me].end]=to;
							Neighbors[me].flag[Neighbors[me].end++]=0;
							
							Neighbors[to].list[Neighbors[to].end]=me;
							Neighbors[to].flag[Neighbors[to].end++]=1;
							
						}
						to=Link[to]; // next atom in this cell
					}
				}
			} 
		}
	}
	
	dirty=false;
}

void VerletListFull::Dump(string fname) {
	ofstream fl(fname.c_str());
	for (int i=0;i<GetNAtom();i++) {
		fl << i <<": ";
		for(int j=0;j<Neighbors[i].end;j++)
			fl << "("<<Neighbors[i].list[j]<<","<<Neighbors[i].flag[j]<< ") ";
		fl << std::endl;
	}
	fl.close();
}

// one-shot...
void VerletListFull::PreIntegration() {
	AllocSize=GetNAtom();
	MemRealloc(Link, AllocSize*sizeof(int));
	MemRealloc(Neighbors, AllocSize*sizeof(NeighborList));
	
	for(int i=0;i<AllocSize;i++) {
		Neighbors[i].start=0.0;
		Neighbors[i].end=0.0;
		Neighbors[i].size=nmean;
		MemAlloc(Neighbors[i].list, nmean*sizeof(int));
		MemAlloc(Neighbors[i].flag, nmean*sizeof(int));
		memset(Neighbors[i].list,-1,nmean*sizeof(int));
		memset(Neighbors[i].flag,-1,nmean*sizeof(int));		
	}
	
	Refresh();
	SetConditionerType(COND_PRE_CALCULATION);
}

void VerletListFull::IterateHalf(MDGadget* IteratedClass) {
	
	if(!half_loop) return;
	
	looping_full=false;
	int na=GetNAtom();
	
	if(dirty) Update();
	
	for(at_idx=0; at_idx<na; at_idx++) { // outer loop...
		if(!CheckActive(at_idx)) continue;
		if(!(IteratedClass->PreIterationNode(at_idx))) continue;
		for(
			nl_idx=Neighbors[at_idx].start;
			nl_idx<Neighbors[at_idx].end;
			nl_idx++
			) 
		{ // inner loop...
			to_idx=Neighbors[at_idx].list[nl_idx];
			if(!CheckActive(to_idx)) continue;
			if((CheckGhost(at_idx)&&CheckGhost(to_idx))) continue;
			IteratedClass->IterationNode(at_idx,to_idx);
		}
	}
}

void VerletListFull::IterateFull(MDGadget* IteratedClass) {
	
	if(!full_loop) return;

	looping_full=true;
	int na=GetNAtom();
	
	if(dirty) Update();
	
	for(at_idx=0; at_idx<na; at_idx++) { // outer loop...
		if(!CheckActive(at_idx)) continue;
		if(!(IteratedClass->PreIterationNode(at_idx))) continue;
		for(
			nl_idx=0; // start from the begining of the list...
			nl_idx<Neighbors[at_idx].end;
			nl_idx++
			) 
		{ // inner loop...
			to_idx=Neighbors[at_idx].list[nl_idx];
			if(!CheckActive(to_idx)) continue;
			if((CheckGhost(at_idx)&&CheckGhost(to_idx))) continue;
			IteratedClass->IterationNode(at_idx,to_idx);
		}
	}
	
	// DEBUG
	/*
	for(int i=0;i<na;i++) {
		std::cerr << i<<": ("<<Atoms(i).x<<","<<Atoms(i).y<<","<<Atoms(i).z<<") ("
		<< Atoms(i).fx<<","<<Atoms(i).fy<<","<<Atoms(i).fz<<")\n";
	}
	
	if(System->Step==5) die("hope..hope..");
	*/
}

void VerletListFull::GetIterationVariables(int& at, int& to, 
										   int& at_nbidx, 
										   NeighborList*& at_nblist)
{
	at=at_idx;
	to=to_idx;
	at_nbidx=nl_idx;
	at_nblist=&(Neighbors[at]);
}

void VerletListFull::PrintInfo(ostream& ost) {
	ost<<"id."<<id<<" "<<get_name()
	<<": rebuild_period="<<RebuildPeriod<<" steps "
	<<"cut_radius="<<CutRadius<<" tolerance="<<RadiusTolerance<<std::endl;
}
