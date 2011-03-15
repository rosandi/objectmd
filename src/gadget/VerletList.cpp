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

#define MEANNB 120

using std::cerr;
using std::ofstream;

VerletList::VerletList() {
	set_name("VERLET LIST");
	register_class(get_name());
	SetConditionerType(COND_PRE_INTEGRATION|COND_PRE_CALCULATION);
	UpdatePeriod=5; 
	Step=0;
	NeighSize=0;
	NeighborList=NULL;
	NeighborIndex=NULL;
	VerletRadius=0.0;
}

void VerletList::CellNumber(int at, int& xid, int& yid, int& zid)
{
	xid=(int)((Atoms[at].x-MDSystem->Box.x0)/VerletRadius);
	yid=(int)((Atoms[at].y-MDSystem->Box.y0)/VerletRadius);
	zid=(int)((Atoms[at].z-MDSystem->Box.z0)/VerletRadius);
	if(xid>(U-1))xid=U-1;if(xid<0)xid=0;
	if(yid>(V-1))yid=V-1;if(yid<0)yid=0;
	if(zid>(W-1))zid=W-1;if(zid<0)zid=0;
}

void VerletList::Init(SimSystem* WorkSys)
{
	Conditioner::Init(WorkSys);
	// Accept double parameter[0] as Verlet radius
	VerletRadius=(par.size()==1)?par[0]:(2.*UpdatePeriod+1)*MDSystem->MaxPath+MDSystem->Forces->MaxCutRadius;
	
	// Take one less cells. Last cell may have more atoms then others!!!
	// note: this imbalance cells content would decrease performance for parallel program!
	// 0 dimension is not allowed!!!
	
	U = (int)(MDSystem->Box.lx/VerletRadius); if(U==0) U=1;
	V = (int)(MDSystem->Box.ly/VerletRadius); if(V==0) V=1;
	W = (int)(MDSystem->Box.lz/VerletRadius); if(W==0) W=1;
	
	if (USE_BOUNDARY(PERIODIC_X) && U<3) THROW("X layer is too small to use periodic boundary condition");
	if (USE_BOUNDARY(PERIODIC_Y) && V<3) THROW("Y layer is too small to use periodic boundary condition");
	if (USE_BOUNDARY(PERIODIC_Z) && W<3) THROW("Z layer is too small to use periodic boundary condition");

	Link          = new int[NAtom];
	NeighborIndex = new int[NAtom+1];
	
	// I dont believe in 'new' for this one..
	NeighSize=NAtom*MEANNB;
	NeighborList  = (int*) malloc(NeighSize*sizeof(int));
	
	// Link it to Forces
	MDSystem->Forces->NeighborIndex=NeighborIndex;
	MDSystem->Forces->NeighborList=NeighborList;
	// refresh it for the first time
	RefreshTable();
}

void VerletList::Reallocate(int NewSize) {
	NeighborList=(int*) realloc(NeighborList, NewSize*sizeof(int));
	if (NeighborList==NULL) THROW("Fail to reallocate neighbor list");
	MDSystem->Forces->NeighborList=NeighborList;
	NeighSize=NewSize;
}

//---------------------Local-macros-for-RefreshTable()----------------------//
#define CHECK_BOUNDARY_Z {if(USE_BOUNDARY(PERIODIC_Z))\
{nk=(k==-1)?(W-1):(k==W)?0:k;}else if(k==-1||k==W)continue;else nk=k;}
#define CHECK_BOUNDARY_Y {if(USE_BOUNDARY(PERIODIC_Y))\
{nj=(j==-1)?(V-1):(j==V)?0:j;}else if(j==-1||j==V)continue;else nj=j;}
#define CHECK_BOUNDARY_X {if(USE_BOUNDARY(PERIODIC_X))\
{ni=(i==-1)?(U-1):(i==U)?0:i;}else if(i==-1||i==U)continue;else ni=i;}							  
//--------------------------------------------------------------------------//

void VerletList::RefreshTable()
{

	int     i, j, k, ni, nj, nk;
	int     at, to, NCounter;
	double  rd, Rsq=VerletRadius*VerletRadius;
	double  xto, yto, zto;
	int     CellHeaders[U][V][W];

	// Initiate the Link cells and header tabel	
	for (at=0;at<NAtom;at++) Link[at]=-1;
	for (int i=0;i<U;i++)for(int j=0;j<V;j++)for(int k=0;k<W;k++)CellHeaders[i][j][k]=-1;
	
	// Find membership of every atom in the system
	// Remember! smallest number in cell must be the header
	// This algorithm ensures that.
	
	for (at=NAtom-1;at>=0;at--) {
		int u, v, w;
		CellNumber(at, u, v, w);
		Link[at]=CellHeaders[u][v][w];
		CellHeaders[u][v][w]=at;
	}
//	DEBUG3INT("Dimension", U, V, W);
//	DEBUGEXIT
	// Create Verlet neighbour list
	// Every atom only search to atoms in it cell and 
	// other 26 neighboring cells
	NCounter=0;

	for (int me=0;me<NAtom;me++) {
		int u, v, w;
		CellNumber(me, u, v, w);
		
		if(me!=CellHeaders[u][v][w]) {
			char st[128];
			sprintf(st, "Wrong header!!\n(%d)(%d, %d, %d)\n", me, u, v, w);
			THROW(st);
		}

		//start=NeighborIndex[ni]; end=NeighborIndex[ni+1];
		CellHeaders[u][v][w]=Link[me];
		NeighborIndex[me]=NCounter;
		for (k=w-1; k<=w+1; k++) {			
			CHECK_BOUNDARY_Z;					
			for (j=v-1; j<=v+1; j++) {
				CHECK_BOUNDARY_Y;
				for (i=u-1; i<=u+1; i++) {
					CHECK_BOUNDARY_X;	
					to=CellHeaders[ni][nj][nk];
					while (to>=0) {
						rd=CalcSqrDistance(me, to, xto, yto, zto);
						if(rd<=Rsq) {
							if (NCounter>=NeighSize) Reallocate(NeighSize+MEANNB);
							NeighborList[NCounter++]=to;
						}
						to=Link[to]; // next atom in this cell
					}
				}
			} 
		}
	}
	NeighborIndex[NAtom]=0;
}

void VerletList::DumpNeighborIndex() {
	ofstream fl("neighborindex.dump");
	for (int i = 0; i < NAtom; i++)
	{
		fl << NeighborIndex[i] << '\n';
	}
	fl.close();
}

VerletList::~VerletList() {
	delete[] Link;
	delete[] NeighborIndex;
	free(NeighborList);
}
