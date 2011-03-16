/* lib
 ************************************************
 * Object-MD
 * Object Oriented Molecular Dynamics Program
 * 
 * (c)Yudi Rosandi, 2005
 *
 * Version 2.0 (06.10.07)
 *
 * Project started on July 2005.
 *
 ************************************************
 * Object-MD header file: NeighborCell
 * has basically NON-PERIODIC boundary.
 * The caller must take care of the periodicity.
 *
*/

#ifndef _NEIGHBOR_CELL_HPP_
#define _NEIGHBOR_CELL_HPP_

#include <cfloat>
#include <omd/iterator.hpp>

/**
 @ingroup iterator
 @brief Neighbor Cell Iterator
 
 This class is an iterator using the neighbor cell algorithm, to
 make the force calculation loop fast. The interator is applicable
 for short ranged inteaction potentials like in metals.
 
**/
 

class NeighborCell: public MDIterator {
	OMD_INT GridNX, GridNY, GridNZ, TotalGrid;
	OMD_FLOAT RCell;
	
	struct CellStruct {
		OMD_INT  NbCell[13];
		OMD_INT  Size;
	} *CellGrid;

	struct IndexList {
		OMD_INT I;
		IndexList* Next;
	} **Cell, **CellTail;
	
	OMD_INT *NMember;
	
	SysBox Box;

public:
	
	// refresh=0 means: never refresh
	NeighborCell(OMD_INT refresh=0, OMD_FLOAT rtol=0.0){
	 	set_name("NEIGHBORCELL");
		register_class(get_name());
		RefreshPeriod=refresh;
	 	SetConditionerType(COND_PRE_INTEGRATION|COND_PRE_CALCULATION);
	 	TotalGrid=0;
	 	CellGrid=NULL;
	 	Cell=CellTail=NULL;
	 	NMember=NULL;
	 	RadiusTolerance=rtol;
	}

	void Release(){
		for(OMD_INT i=0;i<TotalGrid;i++){
			DeleteTail(Cell[i]);
			delete Cell[i];
		}
		MemFree(CellGrid);
		MemFree(Cell);
		MemFree(CellTail);
		MemFree(NMember);
	}

	virtual ~NeighborCell() {Release();}
	
	void Init(MDSystem *WorkSys) {
		MDIterator::Init(WorkSys);
		RCell=WorkSys->GetMaxCutRadius()+RadiusTolerance;
	}

	void CalculateSystemBox(){
		OMD_FLOAT MinX, MaxX, MinY, MaxY, MinZ, MaxZ;        
		OMD_INT natom=GetNAtom();
		MinX=MinY=MinZ= DBL_MAX;
		MaxX=MaxY=MaxZ=-DBL_MAX;
		for(OMD_INT i=0; i<natom; i++) {
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
	}
	
	void CreateGrid() {
		RCell=WorkSys->GetMaxCutRadius()+RadiusTolerance;
		CalculateSystemBox();
		GridNX=OMD_INT(Box.lx/RCell)+1;
		GridNY=OMD_INT(Box.ly/RCell)+1;
		GridNZ=OMD_INT(Box.lz/RCell)+1;
		TotalGrid=GridNX*GridNY*GridNZ;
		
		OMD_SIZET s_cellgrid=TotalGrid*sizeof(CellStruct);
		OMD_SIZET s_cell=TotalGrid*sizeof(IndexList*);
		OMD_SIZET s_celltail=TotalGrid*sizeof(IndexList*);
		OMD_SIZET s_nmember=TotalGrid*sizeof(OMD_INT);

		MemAlloc(CellGrid,s_cellgrid);
		MemAlloc(Cell,s_cell);
		MemAlloc(CellTail,s_celltail);
		MemAlloc(NMember,s_nmember);
		
		for(OMD_INT i=0;i<TotalGrid;i++){Cell[i]=NULL;CellTail[i]=NULL;NMember[i]=0;}
		
		for(OMD_INT i=0;i<TotalGrid;i++)
			for(OMD_INT j=0;j<13;j++) {
				CellGrid[i].NbCell[j]=-1;
				CellGrid[i].Size=0;
			}
		
		// -1 means no neighboring cell	
		for(OMD_INT p=0;p<TotalGrid;p++) {
			OMD_INT n=0;			
			OMD_INT px=p%GridNX;
			OMD_INT py=(p/GridNX)%GridNY;
			OMD_INT pz=(p/(GridNX*GridNY))%GridNZ;
			for(OMD_INT k=-1;k<=1;k++) 
				for(OMD_INT j=-1;j<=1;j++)
					for(OMD_INT i=-1;i<=1;i++) {
						OMD_INT nx=px+i;
						OMD_INT ny=py+j;
						OMD_INT nz=pz+k;
						if((nx<0)||(nx>=GridNX))continue;
						if((ny<0)||(ny>=GridNY))continue;
						if((nz<0)||(nz>=GridNZ))continue;
						
						if(nz<pz)continue;
						if(nz==pz){
							if(ny<py)continue;
							if(ny==py&&nx<=px)continue;
						}
						assert(n<=13, "maximum 13 cell neighbors exceeded");
						CellGrid[p].NbCell[n++]=nz*(GridNX*GridNY)+ny*GridNX+nx;
					}
			CellGrid[p].Size=n;
		}
	}
	
	void Dump(ostream &ost) {
		ost << "box (" << Box.x0<<","<<Box.y0<<","<<Box.z0<<");("
		               << Box.x1<<","<<Box.y1<<","<<Box.z1<<")\n";
		ost << "radius " << RCell << " tolerance " << RadiusTolerance <<"\n";
		ost << "dimension " << GridNX << "," << GridNY << "," << GridNZ << "\n"
		    << "total_grid " << TotalGrid << "\n";

		for(OMD_INT p=0;p<TotalGrid;p++){
			ost << "Cell " <<p << ":\n";
			for(OMD_INT n=0;n<CellGrid[p].Size;n++) {
				ost << "  " << CellGrid[p].NbCell[n] << "\n";
			}
		}
	}
	
	OMD_INT GetCellNumber(OMD_INT at)
	{
		Atom* a=AtomPtr(at);
		OMD_INT xid=(OMD_INT)((a->x-Box.x0)/RCell);
		OMD_INT yid=(OMD_INT)((a->y-Box.y0)/RCell);
		OMD_INT zid=(OMD_INT)((a->z-Box.z0)/RCell);
		if(xid>=GridNX)xid=GridNX-1;if(xid<0)xid=0;
		if(yid>=GridNY)yid=GridNY-1;if(yid<0)yid=0;
		if(zid>=GridNZ)zid=GridNZ-1;if(zid<0)zid=0;
		OMD_INT cn=zid*(GridNX*GridNY)+yid*GridNX+xid;
		assert(cn<TotalGrid, "total number of grid exceeded: "+as_string(cn)+
			   " of "+as_string(TotalGrid));
		return cn;
	}

	void Attach(OMD_INT cn, OMD_INT idx) {		
		if(NMember[cn]==0) {
			if(!CellTail[cn]) {
				Cell[cn]=new IndexList;
				Cell[cn]->Next=NULL;
				CellTail[cn]=Cell[cn];
			}
		} else {
			if(!CellTail[cn]->Next) {
				CellTail[cn]->Next=new IndexList;
				CellTail[cn]->Next->Next=NULL;
			}
			CellTail[cn]=CellTail[cn]->Next;
		}
		CellTail[cn]->I=idx;
		NMember[cn]++;
	}

	void ListReset(OMD_INT cn) {
		CellTail[cn]=Cell[cn];
		NMember[cn]=0;
	}

	void DeleteTail(IndexList* il) {
		if(!il) return;
		IndexList *til=il->Next;
		il->Next=NULL;
		while(til) {
			il=til;
			til=til->Next;
			delete il;
		}
	}

	void Update() { // update cell membership...
		OMD_INT na=GetNAtom();
		for(OMD_INT i=0;i<TotalGrid;i++) ListReset(i);
		for(OMD_INT i=0;i<na;i++) {
			Attach(GetCellNumber(i),i);
		}
		for(OMD_INT i=0;i<TotalGrid;i++){
			DeleteTail(CellTail[i]);
			if((NMember[i]==0)&&(CellTail[i]!=NULL)) {
				delete CellTail[i];
				Cell[i]=CellTail[i]=NULL;
			}
		}
		dirty=false;
	}

	void Refresh(){ // refresh and recreate grid...
		Release();
		CreateGrid();
		dirty=true;
	}

	// One-shot --> only to create grid after data synchronization
	void PreIntegration(){
		CreateGrid();
		SetConditionerType(COND_PRE_CALCULATION);
	}
	//-----------------------------------------------------------//

/**
 * Execute iteration loop on the IteratedClass, which is of type MDGadget.
 * The iteration runs through all the active atoms.
 * Interaction between ghosts is ignored.
 * The neighbor list is updated by default in every call. This behavior can be
 * altered by explicitly unset force_update parameter.
 */

	void Iterate(MDGadget* IteratedClass, bool force_update=false) {

		if(dirty||force_update) Update();
		
		IndexList *atptr,*toptr;
		for(OMD_INT c=0;c<TotalGrid;c++) {
			if(!Cell[c])continue;

			// This cell
			atptr=Cell[c];			
			while(atptr->Next){
				if(CheckActive(atptr->I)){
					if(IteratedClass->PreIterationNode(atptr->I)) {
						toptr=atptr->Next;
						while(toptr) {

							// only active atoms and no interaction between ghosts
							if(CheckActive(toptr->I)&&
							 !(CheckGhost(atptr->I)&&CheckGhost(toptr->I))) 
							 
								{IteratedClass->IterationNode(atptr->I,toptr->I);}

							toptr=toptr->Next;
						}
					}
				}
				atptr=atptr->Next;
			}
			
			// Neighbor cell
			atptr=Cell[c];
			while(atptr) {
				if(CheckActive(atptr->I)){
					if(IteratedClass->PreIterationNode(atptr->I)) {
						for(OMD_INT nc=0;nc<CellGrid[c].Size;nc++) {
							OMD_INT neigcell=CellGrid[c].NbCell[nc];
							toptr=Cell[neigcell];
							while(toptr) {
								
								if(CheckActive(toptr->I)&&
								 !(CheckGhost(atptr->I)&&CheckGhost(toptr->I))) 
								 
									{IteratedClass->IterationNode(atptr->I,toptr->I);}
									
								toptr=toptr->Next;
								
							}
						}
					}
				}
				atptr=atptr->Next;
			}
		}
	}

	void PrintInfo(ostream& ost) {
		ost<<"id."<<id<<" "<<get_name()
		<<": refresh_period="<<RefreshPeriod<<" steps, "
		<<"radius="<<RCell<<" tolerance="<<RadiusTolerance<<std::endl;
	}

};

#endif
