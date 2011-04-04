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
 * Object-MD header file
 *
 *     CommunicationHandler for 
 *     parallel ObjectMD implementation
 *
*/

/**
 * @class CommunicationHandler
 * @ingroup essential
 * @brief A conditioner to handle communication between processors
 * 
 * This class handles message passing between processors in a parallel MD
 * application. It implements basic data transfer functions used by 
 * MDSystemGrid class. The MPI library is used, though the class can be
 * a base for creating another inter-processor communication routines.
 * The class has two sets of buffers, SendBuffer for sending data, and RecvBuffer
 * for receiving. These buffers are manipulated in the main class. 
 * 
 * Albeit its a minimalist encapsulation of message passing, almost all basic 
 * functions needed by Object-MD are implemented: i.e. send and receive, data
 * gathering, broadcast, reduce, and process synchronizing. The main function is 
 * SendReceive function, that sends and receives all data to/from the neigboring
 * cells (27 cell). Exchanging data is performed using conditioner functions 
 * PreIntegration() and ForceModifier().
 * 
 * The class must be instantiated prior to the creation of simulation system, in
 * MDSystemGrid::CreateSystem(). If any of parallel environment information is 
 * needed before the creation stage, then the public function 
 * MDSystemGrid::CommInit() must be explicitely invoked.
 * 
*/

#ifndef _COMHANDLER_H_
#define _COMHANDLER_H_

// funny: these directives must be comment out to let the file scanned
// by doxigen

#include <omd/conditioner.hpp>

#define SIZETAG   0x0100
#define DATATAG   0x0200
#define CMDTAG    0x0300
#define STREAMTAG 0x0400

#define OMD_COMM_CHAR  0
#define OMD_COMM_INT   1
#define OMD_COMM_FLOAT 2

#define OMD_SUM   0
#define OMD_MAX   1
#define OMD_MIN   2

#define DEFAULT_TRANSFER_LENGTH 4096

struct PackageVector {
	int nid;
	OMD_FLOAT x,y,z;
};

struct PackageScalar {
	int nid;
	OMD_FLOAT value;
};

struct PackageSpace {
    int tid;              ///< atom's type id: interaction and container membership
    int gid;             ///< extended id, tagging, etc
    int nid;			   ///< enumerated id
    int flag;            ///< the status and multi-purpose flag
    OMD_FLOAT x,   y,  z;
    OMD_FLOAT vx, vy, vz;
};

class CommunicationHandler:public MDClass {
	int NProc;
	int Rank;
	class MDSystemGrid *System;

	PackageSpace* SpaceSendBuffer[27];
	PackageSpace* SpaceRecvBuffer[27];
	PackageVector* VectorSendBuffer[27];
	PackageVector* VectorRecvBuffer[27];
	PackageScalar* ScalarSendBuffer[27];
	PackageScalar* ScalarRecvBuffer[27];

	int      SendNumber[27];
	int      RecvNumber[27];
	IndexList*   NeigRubix[27];

	string Hostname;
	bool opened;
	int MirrorMask;
	OMD_FLOAT TotalComtime;
	OMD_FLOAT CellRadiusTolerance;

	virtual void DistributeAtomIndexSlab(int slab);

public:
	CommunicationHandler();
	virtual ~CommunicationHandler();

	Atom&    Atoms(int idx);
	Atom*    AtomPtr(int idx);
	int     GetNAtom();

	int GetRank(){return Rank;}
	int GetSize(){return NProc;}

	IndexList** GetRubix(){return NeigRubix;}

	void   SetRadiusTolerance(OMD_FLOAT tole){CellRadiusTolerance=tole;}
	OMD_FLOAT GetRadiusTolerance(){return CellRadiusTolerance;}

	string GetHostname(){return Hostname;}

	virtual void   Link(MDSystemGrid* s);
	virtual bool   CheckLink(){return opened;}
	virtual void   Close();

	virtual void   CollectSendRecvNumber();
	virtual void   DistributeAtomIndex();
	virtual void   SendReceiveData(void** sendpack, void** recvpack, int unitlength);
	virtual void   SendReceive(int mode);
	virtual void   PrepareSpaceBuffers();
	virtual void   PreparePositionBuffers();
	virtual void   PrepareVelocityBuffers();
	virtual void   PrepareForceBuffers();
	virtual void   PrepareAuxBuffers(int aidx);
	virtual void   UnpackSpace();
	virtual void   UnpackPosition();
	virtual void   UnpackVelocity();	
	virtual void   UnpackForce();
	virtual void   UnpackAux(int aidx);
	virtual bool   CheckCellShift(int b, OMD_FLOAT& xshift, OMD_FLOAT& yshift, OMD_FLOAT& zshift);
	virtual void   SyncProcesses();
	virtual void   Abort();

	virtual int    RawSend(int toproc, void *data, int length);
	virtual int    RawReceive(int fromproc, void* data, int length=DEFAULT_TRANSFER_LENGTH);
	virtual void   RootReduceSUM(OMD_FLOAT*source,OMD_FLOAT*dest,int length);
	virtual int    TakeSUM(int a);
	virtual int    TakeMAX(int a);
	virtual int    TakeMIN(int a);
	virtual OMD_FLOAT TakeSUM(OMD_FLOAT a);
	virtual OMD_FLOAT TakeMAX(OMD_FLOAT a);
	virtual OMD_FLOAT TakeMIN(OMD_FLOAT a);
	virtual void   Gather(void *idata, void *buffer, int length);
	virtual void   AllGather(void *send, void *receive, int senlen, int reclen);
	virtual void   Broadcast(void* a, int size);

};

#endif
