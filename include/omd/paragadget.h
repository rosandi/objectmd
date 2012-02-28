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
 *    ParallelGadget --> Base class of parallel gadgets (multiple inherited)
 */

#ifndef _PARALLEL_GADGET_
#define _PARALLEL_GADGET_

#include <omd/omdtool.h>
#include <omd/config.h>
#include <omd/systemgrid.h>
#include <omd/comhandler.h>

namespace omd {

/**
 * @ingroup baseclass
 * @brief Base class for parallel gadgets
 *
 * This class provides access to the systems communication handler for gadgets 
 * (Conditioners & Detectors). The class must be multiple inherited by gadget child
 * classes, to make it parallel compliant.
 * 
 *
 */


class ParallelGadget {
	MDSystem *sysptr;
	CommunicationHandler* Communicator;
	
public:
	ParallelGadget() {Communicator=NULL;}

	void Init(MDSystem* WorkSys) {
		if(WorkSys->type_of("simulation_system_grid")) {
			Communicator=dynamic_cast<MDSystemGrid*>(WorkSys)->GetCommunicator();
			mdassert(Communicator, "Communication handler is not initialized");
		}
		sysptr=WorkSys;
	}
		
	int RawSend(int toproc, void *data, int length) {
		if(Communicator) return Communicator->RawSend(toproc, data, length);
		else die("Attempt to use RawSend() in a non parallel OMD application");
		return 0;
	}
	
	int RawReceive(int fromproc, void* data, int length=DEFAULT_TRANSFER_LENGTH) {
		if(Communicator) return Communicator->RawReceive(fromproc, data, length);
		else die("Attempt to use RawReceive() in a non parallel OMD application");
		return 0;
	}
	
	void RootReduceSUM(double*source,double*dest,int length){
		if(Communicator) Communicator->RootReduceSUM(source,dest,length);
		else die("Attempt to use RootReduceSUM() in a non parallel OMD application");

	}
	void Gather(void *idata, void *buffer, int length){
		if(Communicator) Communicator->Gather(idata,buffer,length);
		else die("Attempt to use Gather() in a non parallel OMD application");
	}
	

	void AllGather(void *send, void *receive, int senlen, int reclen){
		if(Communicator) Communicator->AllGather(send,receive,senlen,reclen);
		else die("Attempt to use AllGather() in a non parallel OMD application");
	}
	
	// may be use in serial OMD:

	int GetRank() {
		if(Communicator) return Communicator->GetRank();
		return 0;
	}

	
	int TakeSUM(int a) {
		if(Communicator) return Communicator->TakeSUM(a);
		return a;
	}
	
	int TakeMAX(int a) {
		if(Communicator) return Communicator->TakeMAX(a);
		return a;
	}
	
	int TakeMIN(int a) {
		if(Communicator) return Communicator->TakeMIN(a);
		return a;
	}
	
	double TakeSUM(double a) {
		if(Communicator) return Communicator->TakeSUM(a);
		return a;
	}
	
	double TakeMAX(double a) {
		if(Communicator) return Communicator->TakeMAX(a);
		return a;
	}
	
	double TakeMIN(double a) {
		if(Communicator) return Communicator->TakeMIN(a);
		return a;
	}

	void Broadcast(void* a, int size){
		if(Communicator) Communicator->Broadcast(a,size);
	}

	void Broadcast(int& a){
		if(Communicator) Communicator->Broadcast((void*)&a,sizeof(int));
	}

	void Broadcast(double& a){
		if(Communicator) Communicator->Broadcast((void*)&a,sizeof(double));
	}

	void SyncProcesses() {
		if(Communicator) Communicator->SyncProcesses();
	}

};

#define _PARALLEL_COMPLIANT_ public ParallelGadget

}

#endif
