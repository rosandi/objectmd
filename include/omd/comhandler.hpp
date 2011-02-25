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
	OMD_INT nid;
	OMD_FLOAT x,y,z;
};

struct PackageScalar {
	OMD_INT nid;
	OMD_FLOAT value;
};

struct PackageSpace {
    OMD_CHAR  id;              ///< atom's type id: interaction and container membership
    OMD_CHAR  xid;             ///< extended id, tagging, etc
    OMD_INT   nid;			   ///< enumerated id
    OMD_INT   flag;            ///< the status and multi-purpose flag
    OMD_FLOAT x,   y,  z;
    OMD_FLOAT vx, vy, vz;
};

class CommunicationHandler:public MDClass {
	OMD_INT NProc;
	OMD_INT Rank;
	class MDSystemGrid *System;

	PackageSpace* SpaceSendBuffer[27];
	PackageSpace* SpaceRecvBuffer[27];
	PackageVector* ForceSendBuffer[27];
	PackageVector* ForceRecvBuffer[27];
	PackageScalar* AuxSendBuffer[27];
	PackageScalar* AuxRecvBuffer[27];

	OMD_INT      SendNumber[27];
	OMD_INT      RecvNumber[27];
	IndexList*   NeigRubix[27];

	string Hostname;
	bool opened;
	OMD_INT MirrorMask;
	OMD_FLOAT TotalComtime;
	OMD_FLOAT CellRadiusTolerance;

	virtual void DistributeAtomIndexSlab(int slab);

public:
	CommunicationHandler();
	virtual ~CommunicationHandler();

	Atom&    Atoms(OMD_INT idx);
	Atom*    AtomPtr(OMD_INT idx);
	OMD_SIZET     GetNAtom();

	OMD_INT GetRank(){return Rank;}
	OMD_INT GetSize(){return NProc;}

	IndexList** GetRubix(){return NeigRubix;}

	void   SetRadiusTolerance(OMD_FLOAT tole){CellRadiusTolerance=tole;}
	OMD_FLOAT GetRadiusTolerance(){return CellRadiusTolerance;}

	string GetHostname(){return Hostname;}

	virtual void   Link(MDSystemGrid* s);
	virtual bool   CheckLink(){return opened;}
	virtual void   Close();

	virtual void   CollectSendRecvNumber();
	virtual void   DistributeAtomIndex();
	virtual void   PrepareSendBuffers(OMD_INT mode);
	virtual void   SendReceiveData(void** sendpack, void** recvpack, OMD_SIZET unitlength);
	virtual void   SendReceive(OMD_INT mode);
	virtual void   UnpackSpace();
	virtual void   UnpackForce();
	virtual void   UnpackAux(OMD_INT aidx);
	virtual void   Unpack(OMD_INT mode);
	virtual void   SyncProcesses();
	virtual void   Abort();

	virtual OMD_INT    RawSend(OMD_INT toproc, void *data, OMD_INT length);
	virtual OMD_INT    RawReceive(OMD_INT fromproc, void* data, OMD_INT length=DEFAULT_TRANSFER_LENGTH);
	virtual void   RootReduceSUM(OMD_FLOAT*source,OMD_FLOAT*dest,OMD_INT length);
	virtual OMD_INT    TakeSUM(OMD_INT a);
	virtual OMD_INT    TakeMAX(OMD_INT a);
	virtual OMD_INT    TakeMIN(OMD_INT a);
	virtual OMD_FLOAT TakeSUM(OMD_FLOAT a);
	virtual OMD_FLOAT TakeMAX(OMD_FLOAT a);
	virtual OMD_FLOAT TakeMIN(OMD_FLOAT a);
	virtual void   Gather(void *idata, void *buffer, OMD_INT length);
	virtual void   AllGather(void *send, void *receive, OMD_INT senlen, OMD_INT reclen);
	virtual void   Broadcast(void* a, OMD_INT size);

};

#endif
