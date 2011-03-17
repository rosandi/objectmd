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
 *        MDSystemGrid
 *
 *        This header file contains the definition of
 *        The main class in the library (MDSystemGrid).
 *
*/

#ifndef _SIMGRID_H_
#define _SIMGRID_H_

#define OMD_GRID

#include <iostream>
#include <omd/system.hpp>
#include <omd/comhandler.hpp>

#define MYSELF 13
#define CTRDOWN 4
#define CTRUP   21

#define GRID_AUTOX    -1
#define GRID_AUTOY    -2
#define GRID_AUTOZ    -3

/**
 * \brief Parallel version of Object-MD main class.
 * 
 * 
 */

class MDSystemGrid: public MDSystem {	
	friend class CommunicationHandler;


public:
	// negative values needed for GRID_AUTO?
	int ClusterNX, ClusterNY, ClusterNZ;

protected:
	
	string BinDirectory;

	/**
	 * \brief Information about grid configuration
	 * 
	 * NeighborList: the list of rank of the neighboring processors. 
	 * Cell 13 is the rank of current processor (MYSELF)
	 * Border: the border coordinates of the cell (bottom-south-west, top-north-east)
	 * CellX,Y,Z: the absolute cell coordinate in the grid
	 * AtomNames: The names of atom containers in the system
	 * NumberOfContainers: number of containers
	 * TotalAtom: total number of atoms
	 * LCellX,Y,Z: side lenghts of the grid cells.
	 * Box: The simulation box.
	 *
	 */

	struct StructInfo {
		int   NeighborList[MAXPROC][27];
		SysBox    Border[MAXPROC];
		int CellX[MAXPROC],CellY[MAXPROC],CellZ[MAXPROC];
		char  AtomNames[MAXATOMTYPE][32];
		int NumberOfContainers;
		int TotalAtom;
		OMD_FLOAT LCellX,LCellY,LCellZ;
		SysBox    Box;
	} ProcInfo;

	CommunicationHandler *Communicator;
	AtomKeeper GhostBuffer;
	AtomKeeper LocalBuffer;

	int  LocalAtomNumber;
	int CommRefreshPeriod;
	bool FirstSync;

public:

	//----Constructions and destructions----//
	MDSystemGrid(int &argc, char** &argv, int nx=1, int ny=1, int nz=1);
	MDSystemGrid();
	
	virtual ~MDSystemGrid();

	void CommInit();

	//----utilities----//
	int GetGridSize(){return Communicator->GetSize();}	
	int GetRank(){return Communicator->GetRank();}

	void SetClusterArch(int nx, int ny=1, int nz=1){
		ClusterNX=nx;ClusterNY=ny;ClusterNZ=nz;
	}
	
	int GetCellIndex(OMD_FLOAT x, OMD_FLOAT y, OMD_FLOAT z);
	bool CheckOwnership(OMD_FLOAT x, OMD_FLOAT y, OMD_FLOAT z);

	/**
	 * \brief Inquiry the rank of neighboring processor
	 * 
	 * Parameter the relative index of the neighbor:
	 * [0..8]=lower
	 * [9..17]=middle
	 * [18..26]=upper
	 * [13] = the current processor (in the middle slab).
	 * 
	 * negative value means no neighbor, othewise the rank.
	 * 
	 */

	int GetNeighborRank(int relcoord){
		return ProcInfo.NeighborList[Communicator->GetRank()][relcoord];
	}
	
	SysBox& GetCellBorder(int rank=-1){
		if(rank<0)rank=Communicator->GetRank();
		return ProcInfo.Border[rank];
	}
	
	CommunicationHandler* GetCommunicator(){return Communicator;}
	virtual void PrintHeader(ostream& ost){if(GetRank()==0)MDSystem::PrintHeader(ost);}	
	virtual void CheckBoundary(){}; // boundary checking is not needed...
	
	// cancels boundary's minimum distance convention: image atoms copy is used
	virtual void BoundaryCorrectDistances(OMD_FLOAT& dx, OMD_FLOAT& dy, OMD_FLOAT& dz){}
	
	//----Master(Root) routines----//
	virtual void Root_ArrangeNeighbor();	
	virtual void Root_DistributeAtoms();
	virtual void Root_Prepare();

	//----Communication routines----//
	virtual bool CheckArch();
	virtual bool CheckComm();

	virtual void FlattenAtomBox();
//	virtual void Unpack(int package);
	
	//----Creation----//
	virtual void LoadAtoms();
	virtual void PrintInfo(ostream& ost);
	
	virtual void PrintInfo(string fname);
	
	virtual void PrintTime(ostream& ost){if(GetRank()==ROOT) MDSystem::PrintTime(ost);}
	
	virtual AtomContainer* DumpAtoms(AtomKeeper& ak, 
	                       string fname,
	                       int mode=0,
	                       bool* AuxPrintable=NULL,
	                       char* AuxFormat[]=NULL,
	                       string AuxNames=""); 

	virtual AtomContainer* DumpAtoms(
	                       string fname,
	                       int mode=0, 
	                       bool* AuxPrintable=NULL,
	                       char* AuxFormat[]=NULL,
	                       string AuxNames="");
	
	
	        void SetBinDirectory(string bindir);

	virtual AtomContainer* Save(string binname, string mode="w");
	virtual void SaveSimulation(string binfile="");
	virtual void LoadVariables(FILE* fl);
	virtual void CreationFunction();
	virtual void InitGadgets();
	virtual void SyncData(int syncmode);
	virtual void SyncEnvironment();
	virtual void SyncVariable();
	virtual bool CheckRun();
	virtual void AdjustSystem();
	virtual void FirstRun();
	virtual void Initiate();
	virtual void ErrorHandler(const char* errst);
	virtual void MeasurePotential();
	virtual void MeasureKinetic();
	virtual void PrintMessages(ostream& ost);
	virtual void SetCommRefreshPeriod(int peri) {CommRefreshPeriod=peri;}
	string  GetGridConfiguration();
	virtual void ReadParameters();
	virtual int GetLocalAtomNumber(){return LocalAtomNumber;}
	virtual int GetTotalAtomNumber(){return ProcInfo.TotalAtom;}
	virtual void DistributeContainers();
	virtual void UpdateRadiusTolerance();
};

#endif
