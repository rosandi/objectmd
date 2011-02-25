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

#include <omd/config.hpp>
#include <omd/systemgrid.hpp>
#include <omd/comhandler.hpp>

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
			WorkSys->assert(Communicator, "Communication handler is not initialized");
		}
		sysptr=WorkSys;
	}
		
	OMD_INT RawSend(OMD_INT toproc, void *data, OMD_INT length) {
		if(Communicator) return Communicator->RawSend(toproc, data, length);
		else sysptr->die("Attempt to use RawSend() in a non parallel OMD application");
		return 0;
	}
	
	OMD_INT RawReceive(OMD_INT fromproc, void* data, OMD_INT length=DEFAULT_TRANSFER_LENGTH) {
		if(Communicator) return Communicator->RawReceive(fromproc, data, length);
		else sysptr->die("Attempt to use RawReceive() in a non parallel OMD application");
		return 0;
	}
	
	void RootReduceSUM(OMD_FLOAT*source,OMD_FLOAT*dest,OMD_INT length){
		if(Communicator) Communicator->RootReduceSUM(source,dest,length);
		else sysptr->die("Attempt to use RootReduceSUM() in a non parallel OMD application");

	}
	void Gather(void *idata, void *buffer, OMD_INT length){
		if(Communicator) Communicator->Gather(idata,buffer,length);
		else sysptr->die("Attempt to use Gather() in a non parallel OMD application");
	}
	

	void AllGather(void *send, void *receive, OMD_INT senlen, OMD_INT reclen){
		if(Communicator) Communicator->AllGather(send,receive,senlen,reclen);
		else sysptr->die("Attempt to use AllGather() in a non parallel OMD application");
	}
	
	// may be use in serial OMD:

	OMD_INT GetRank() {
		if(Communicator) return Communicator->GetRank();
		return 0;
	}

	
	OMD_INT TakeSUM(OMD_INT a) {
		if(Communicator) return Communicator->TakeSUM(a);
		return a;
	}
	
	OMD_INT TakeMAX(OMD_INT a) {
		if(Communicator) return Communicator->TakeMAX(a);
		return a;
	}
	
	OMD_INT TakeMIN(OMD_INT a) {
		if(Communicator) return Communicator->TakeMIN(a);
		return a;
	}
	
	OMD_FLOAT TakeSUM(OMD_FLOAT a) {
		if(Communicator) return Communicator->TakeSUM(a);
		return a;
	}
	
	OMD_FLOAT TakeMAX(OMD_FLOAT a) {
		if(Communicator) return Communicator->TakeMAX(a);
		return a;
	}
	
	OMD_FLOAT TakeMIN(OMD_FLOAT a) {
		if(Communicator) return Communicator->TakeMIN(a);
		return a;
	}

	void Broadcast(void* a, OMD_INT size){
		if(Communicator) Communicator->Broadcast(a,size);
	}

	void Broadcast(int& a){
		if(Communicator) Communicator->Broadcast((void*)&a,sizeof(OMD_INT));
	}

	void Broadcast(OMD_FLOAT& a){
		if(Communicator) Communicator->Broadcast((void*)&a,sizeof(OMD_FLOAT));
	}

	void SyncProcesses() {
		if(Communicator) Communicator->SyncProcesses();
	}

};

#endif
