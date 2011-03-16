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

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <conditioner/VerletList.hpp>

using std::cerr;
using std::ofstream;

VerletList::VerletList() {
	set_name("VERLET LIST");
	register_class(get_name());
	SetConditionerType(COND_PRE_INTEGRATION|COND_PRE_CALCULATION);
	NeighSize=0;
	AllocSize=0;
	nmean=120;
	Link=NULL;
	NeighborList=NULL;
	NeighborIndex=NULL;
}

VerletList::~VerletList() {
	MemFree(Link);
	MemFree(NeighborIndex);
	MemFree(NeighborList);
}

void VerletList::ReadParameter() {
	SysParam->peek("verlet.update", UpdatePeriod, 5);
	SysParam->peek("verlet.refresh", RebuildPeriod, 0);
	SysParam->peek("verlet.radtole", RadiusTolerance, -1.0);
	SysParam->peek("verlet.n_mean", nmean, 120);
}

bool VerletList::CheckParameter() {
	if(RadiusTolerance<=0.0) {
		warn("no radius tolerance for verlet list. (verlet.radtole)");
		RadiusTolerance=0.0;
	}
	return true;
}

void VerletList::Init(MDSystem* WorkSys) {
	Conditioner::Init(WorkSys);
	
	if (System->PBoundary) {
		assert(System->type_of("simulation_system_grid"), 
			   "MDSystemGrid is required for periodic boundary simulation");
	}

	CalculateBox();	
}

void VerletList::CellNumber(OMD_SIZET at, OMD_SIZET& xid, OMD_SIZET& yid, OMD_SIZET& zid) {
	OMD_INT vxid=(OMD_INT)((Atoms(at).x-Box.x0)/VerletRadius);
	OMD_INT vyid=(OMD_INT)((Atoms(at).y-Box.y0)/VerletRadius);
	OMD_INT vzid=(OMD_INT)((Atoms(at).z-Box.z0)/VerletRadius);
	if(vxid>(U-1))vxid=U-1;if(vxid<0)vxid=0;
	if(vyid>(V-1))vyid=V-1;if(vyid<0)vyid=0;
	if(vzid>(W-1))vzid=W-1;if(vzid<0)vzid=0;
	xid=vxid;yid=vyid;zid=vzid;
}

void VerletList::CalculateBox() {
	OMD_FLOAT MinX, MaxX, MinY, MaxY, MinZ, MaxZ;        
	OMD_INT na=GetNAtom();
	MinX=MinY=MinZ= DBL_MAX;
	MaxX=MaxY=MaxZ=-DBL_MAX;
	
	VerletRadius=System->GetMaxCutRadius()+RadiusTolerance;

	for(OMD_INT i=0; i<na; i++) {
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
}

void VerletList::Refresh(){
	CalculateBox();
	SetDirty();
}

void VerletList::Update() {

	OMD_SIZET i, j, k, ni, nj, nk;
	OMD_SIZET at, to, NCounter;
	double  rd, Rsq=VerletRadius*VerletRadius;
	double  xto, yto, zto;
	
	OMD_SIZET na=GetNAtom();
	
	if(AllocSize<na) {
		MemRealloc(Link, na*sizeof(int));
		MemRealloc(NeighborIndex, (na+1)*sizeof(int));
		AllocSize=na;
	}

	OMD_INT CellHeaders[U][V][W];

	// Initiate the Link cells and header tabel	
	for (at=0;at<na;at++) Link[at]=-1;
	for (int i=0;i<U;i++)
		for(int j=0;j<V;j++)
			for(int k=0;k<W;k++) CellHeaders[i][j][k]=-1;
	
	// Find membership of every atom in the system
	// Remember! smallest number in cell must be the header
	// This algorithm ensures that.
	
	for (at=na-1;at>=0;at--) {
		OMD_SIZET u, v, w;
		CellNumber(at, u, v, w);
		Link[at]=CellHeaders[u][v][w];
		CellHeaders[u][v][w]=at;
	}

	// Create Verlet neighbour list
	// Every atom only search to atoms in it cell and 
	// other 26 neighboring cells
	NCounter=0;
	for (OMD_SIZET me=0;me<na;me++) {
		OMD_SIZET u, v, w;
		CellNumber(me, u, v, w);
		
		if(me!=(OMD_SIZET)CellHeaders[u][v][w]) 
			die("wrong header at index "+as_string(me)+
			    " cell index ("+as_string(u)+","+as_string(v)+","+as_string(w)+")");

		//start=NeighborIndex[ni]; end=NeighborIndex[ni+1];
		CellHeaders[u][v][w]=Link[me];
		NeighborIndex[me]=NCounter;
		for (k=w-1; k<=w+1; k++) {			
			for (j=v-1; j<=v+1; j++) {
				for (i=u-1; i<=u+1; i++) {
					to=CellHeaders[ni][nj][nk];
					while (to>=0) {
						rd=CalcSqrDistance(me, to, xto, yto, zto);
						if(rd<=Rsq) {
							if(NCounter>=NeighSize) {
								NeighSize+=nmean;
								MemRealloc(NeighborList, NeighSize*sizeof(int));
							}
							NeighborList[NCounter++]=to;
						}
						to=Link[to]; // next atom in this cell
					}
				}
			} 
		}
	}
	NeighborIndex[na]=0;
}

void VerletList::Dump(string fname) {
	ofstream fl(fname.c_str());
	for (OMD_SIZET i = 0; i < GetNAtom(); i++) fl << NeighborIndex[i] << '\n';
	fl.close();
}

// one-shot...
void VerletList::PreIntegration() {
	CalculateBox();
	AllocSize=GetNAtom();
	NeighSize=AllocSize*nmean;
	MemRealloc(Link, AllocSize*sizeof(int));
	MemRealloc(NeighborIndex, (AllocSize+1)*sizeof(int));
	MemRealloc(NeighborList,NeighSize*sizeof(int));
	SetConditionerType(COND_PRE_CALCULATION);
	Refresh();
}

void VerletList::GetNeighborIndex(OMD_SIZET ni, OMD_SIZET &start, OMD_SIZET &end) {
	start=NeighborIndex[ni];
	end=NeighborIndex[ni+1];
}

int VerletList::GetNeighbor(OMD_SIZET ni) {
	return NeighborList[ni];
}

void VerletList::Iterate(MDGadget* IteratedClass, bool force_update) {
	if(dirty||force_update) Update();
	OMD_SIZET na=GetNAtom();
	OMD_SIZET start, end;
	
	for(OMD_SIZET at=0; at<na; at++) {
		if(!CheckActive(at)) continue;
		if(!(IteratedClass->PreIterationNode(at))) continue;
		GetNeighborIndex(at,start,end);
		for(OMD_SIZET i=start;i<end;i++) {
			OMD_SIZET to=GetNeighbor(i);
			if(!CheckActive(to)) continue;
			if((CheckGhost(at)&&CheckGhost(to))) continue;
			IteratedClass->IterationNode(at,to);
		}
	}
}
