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
 * Implementation of CommunicationHandler class
 *
*/

#include <iostream>
#include <iomanip>
#include <ctime>
#include <sys/time.h>
#include <mpi.h>
#include <omd/systemgrid.hpp>
#include <omd/comhandler.hpp>

// rubix offset:
//   0->down_slab
//   9->middle_slab (here)
//  18->up

#define BottomSlab 0
#define MidSlab 9
#define TopSlab 18

CommunicationHandler::CommunicationHandler(){
	opened=false;
	System=NULL;
	Hostname[0]='\0';
	for(int i=0;i<27;i++) {
		SpaceSendBuffer[i]=NULL;
		SpaceRecvBuffer[i]=NULL;
		VectorSendBuffer[i]=NULL;
		VectorRecvBuffer[i]=NULL;
		ScalarSendBuffer[i]=NULL;
		ScalarRecvBuffer[i]=NULL;
		SendNumber[i]=0;
		RecvNumber[i]=0;
		NeigRubix[i]=new IndexList;
	}

	set_name("CommunicationHandler:MPI");
	register_class("communication_handler");
	TotalComtime=0.0;
	CellRadiusTolerance=0.0;
}


CommunicationHandler::~CommunicationHandler() {
	if(!opened) return;

	for(int i=0;i<27;i++) {
		disable_log(LOGWARNING);
		MemFree(SpaceSendBuffer[i]);
		MemFree(SpaceRecvBuffer[i]);
		MemFree(VectorSendBuffer[i]);
		MemFree(VectorRecvBuffer[i]);
		MemFree(ScalarSendBuffer[i]);
		MemFree(ScalarRecvBuffer[i]);
		enable_log(LOGWARNING);
	}

	char st[DEFAULT_TRANSFER_LENGTH];

	std::cerr.flush();	
	std::cerr << std::fixed << std::setprecision(4);
	
	if(GetRank()==ROOT){
		std::cerr << "\n(Root) Waiting for child processes to finish\n"
				  << "[\n ElapsedTime: " << System->ElapsedTime
				  << "\n Step: " << System->Step << "\n]\n";
	}
	
	SyncProcesses();
	OMD_FLOAT walltime=TakeMAX(System->SimWallTime);
	OMD_FLOAT comtime=TakeMAX(TotalComtime);
	string tmsg("walltime: ");
	tmsg.append(as_string(walltime)+
				" commtime: "+as_string(comtime)+" ("+as_string(comtime/walltime)+") (seconds)");

	if(GetRank()==ROOT) {
		for(int r=1;r<NProc;r++) {
			RawReceive(r,st,DEFAULT_TRANSFER_LENGTH);
			if(strncmp(st,"FINISHED",8)==0)
				blog("Received from proc "+as_string(r)+" Message=\""+st+"\"");
		}
		std::cerr << "Simulation walltime: " << walltime << " seconds\n"
		          << "Communication time: " << comtime << " seconds\n";
		System->PrintTime(std::cerr);
	} else {
		sprintf(st, "FINISHED: PROC %d",Rank);
		RawSend(ROOT,st,DEFAULT_TRANSFER_LENGTH);
	}
	
	blog(tmsg);
	blog("SIMULATION ENDS");
	SyncProcesses();
	Close();
}

void CommunicationHandler::Link(MDSystemGrid* s) {
	int ln;

	System=s;

	char hh[MPI_MAX_PROCESSOR_NAME];
	assert(s->Argc&&s->Argv, "no valid arguments in the caller class");
	MPI_Init(s->Argc,s->Argv);
	MPI_Comm_size(MPI_COMM_WORLD, &NProc);
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
	MPI_Get_processor_name(hh, &ln);
	int nreqp=s->ClusterNX*s->ClusterNY*s->ClusterNZ;
	
	assert(nreqp<=NProc,
		   "failed to initiate sufficient number of processors: proc_required="+as_string(nreqp)+
		   " MPI gives "+as_string(NProc)+" processor(s)");
	
	Hostname.assign(hh);
	opened=true;
}

void CommunicationHandler::Close() {
	MPI_Finalize();
	opened=false;
}

int CommunicationHandler::RawSend(int toproc, void *data, int length){
	MPI_Send(data, length, MPI_CHAR, toproc, STREAMTAG, MPI_COMM_WORLD);
	return length;
}

int CommunicationHandler::RawReceive(int fromproc, void* data, int length){
	MPI_Status sta;
	int received;
	MPI_Recv(data, length, MPI_CHAR, fromproc, STREAMTAG, MPI_COMM_WORLD, &sta);
	MPI_Get_count(&sta, MPI_CHAR, &received);
	return received;
	
}

void CommunicationHandler::RootReduceSUM(OMD_FLOAT*source, OMD_FLOAT*dest, int length){
	MPI_Reduce(source,dest,length,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
}

void CommunicationHandler::PrepareSpaceBuffers(){
	for(int i=0;i<27;i++){
		int na=NeigRubix[i]->length;
		int a;
		if(i==MYSELF||na<=0) continue;
		if(System->GetNeighborRank(i)<0)continue;
		MemRealloc(SpaceSendBuffer[i], na*sizeof(PackageSpace));
		PackageSpace* X=SpaceSendBuffer[i];
		NeigRubix[i]->reset();
		while((a=NeigRubix[i]->fetch())>=0) {
			Atom* A=System->AtomPtr(a);
			X->id=A->id;
			X->xid=A->xid;
			X->nid=A->nid;
			X->flag=A->flag;
			X->x=A->x;
			X->y=A->y;
			X->z=A->z;
			X->vx=A->vx;
			X->vy=A->vy;
			X->vz=A->vz;
			X++;
		}
	}
}

void CommunicationHandler::PreparePositionBuffers(){
	for(int i=0;i<27;i++){
		int na=NeigRubix[i]->length;
		int a;
		if(i==MYSELF||na<=0) continue;
		if(System->GetNeighborRank(i)<0)continue;
		MemRealloc(VectorSendBuffer[i], na*sizeof(PackageVector));
		PackageVector* X=VectorSendBuffer[i];
		NeigRubix[i]->reset();
		while((a=NeigRubix[i]->fetch())>=0) {
			Atom* A=System->AtomPtr(a);
			X->nid=A->nid;
			X->x=A->x;
			X->y=A->y;
			X->z=A->z;
			X++;
		}			
	}
}
		
void CommunicationHandler::PrepareVelocityBuffers(){
	for(int i=0;i<27;i++){
		int na=NeigRubix[i]->length;
		int a;
		if(i==MYSELF||na<=0) continue;
		if(System->GetNeighborRank(i)<0)continue;
		MemRealloc(VectorSendBuffer[i], na*sizeof(PackageVector));
		PackageVector* X=VectorSendBuffer[i];
		NeigRubix[i]->reset();
		while((a=NeigRubix[i]->fetch())>=0) {
			Atom* A=System->AtomPtr(a);
			X->nid=A->nid;
			X->x=A->vx;
			X->y=A->vy;
			X->z=A->vz;
			X++;
		}			
	}		
}
void CommunicationHandler::PrepareForceBuffers(){
	for(int i=0;i<27;i++){
		int na=NeigRubix[i]->length;
		int a;
		if(i==MYSELF||na<=0) continue;
		if(System->GetNeighborRank(i)<0)continue;		
		MemRealloc(VectorSendBuffer[i], na*sizeof(PackageVector));
		PackageVector* X=VectorSendBuffer[i];
		NeigRubix[i]->reset();
		while((a=NeigRubix[i]->fetch())>=0) {
			Atom* A=System->AtomPtr(a);
			X->nid=A->nid;
			X->x=A->fx;
			X->y=A->fy;
			X->z=A->fz;
			X++;
		}
	}
}

void CommunicationHandler::PrepareAuxBuffers(int aidx){
	for(int i=0;i<27;i++){
		int na=NeigRubix[i]->length;
		int a;
		if(i==MYSELF||na<=0) continue;
		if(System->GetNeighborRank(i)<0)continue;
		if(aidx>MAXAUXVAR) die("attempt to send/receive invalid auxvar index="+as_string(aidx));
		MemRealloc(ScalarSendBuffer[i], na*sizeof(PackageScalar));
		PackageScalar* X=ScalarSendBuffer[i];
		NeigRubix[i]->reset();
		while((a=NeigRubix[i]->fetch())>=0) {
			Atom* A=System->AtomPtr(a);
			X->nid=A->nid;
			X->value=A->aux[aidx];
			X++;
		}
	}
}


void CommunicationHandler::CollectSendRecvNumber() {
	MPI_Status  stat[54];
	MPI_Request req[54];
	int nreq=0;

	memset(SendNumber,0,27*sizeof(int));
	memset(RecvNumber,0,27*sizeof(int));

	for(int i=0;i<27;i++){
		int proc=System->GetNeighborRank(i);
		if(i==MYSELF||proc<0)continue;
		SendNumber[i]=NeigRubix[i]->length;
		MPI_Send_init(&SendNumber[i],1,MPI_INT,proc,SIZETAG+(26-i),MPI_COMM_WORLD,&req[nreq++]);
		MPI_Recv_init(&RecvNumber[i],1,MPI_INT,proc,SIZETAG+i,MPI_COMM_WORLD, &req[nreq++]);
	}

	MPI_Startall(nreq, req);
	MPI_Waitall(nreq,req,stat);
	// FIXME! check statuses here
}


bool CommunicationHandler::CheckCellShift(int b, OMD_FLOAT& xshift,OMD_FLOAT& yshift, OMD_FLOAT& zshift) {
	int rank=GetRank();
	// my grid coordinate
	int mx=System->ProcInfo.CellX[rank];
	int my=System->ProcInfo.CellY[rank];
	int mz=System->ProcInfo.CellZ[rank];

	// neigboring cell coordinate relative to the current cell
	int nx=b%3;
	int ny=(b/3)%3;
	int nz=(b/9)%3;

	bool do_shift=false;

	// this occurs only in periodic boundary condition,
	// wraps on a single processor
	if(System->GetNeighborRank(b)==rank) {

		xshift=(nx==0)?(-System->ProcInfo.Box.lx):((nx==2)?System->ProcInfo.Box.lx:0.0);
		yshift=(ny==0)?(-System->ProcInfo.Box.ly):((ny==2)?System->ProcInfo.Box.ly:0.0);
		zshift=(nz==0)?(-System->ProcInfo.Box.lz):((nz==2)?System->ProcInfo.Box.lz:0.0);
		do_shift=true;

	} else {

		// this happens only if PERIODIC_X|Y|Z
		// otherwise neighbor rank b < 0...

		if(mx==0&&nx==0)
		{xshift=-System->ProcInfo.Box.lx; do_shift=true;}
		else if((mx==(int)System->ClusterNX-1)&&(nx==2))
		{xshift=System->ProcInfo.Box.lx; do_shift=true;}
		else xshift=0.0;

		if(my==0&&ny==0)
		{yshift=-System->ProcInfo.Box.ly;do_shift=true;}
		else if((my==(int)System->ClusterNY-1)&&(ny==2))
		{yshift=System->ProcInfo.Box.ly;do_shift=true;}
		else yshift=0.0;

		if(mz==0&&nz==0)
		{zshift=-System->ProcInfo.Box.lz;do_shift=true;}
		else if((mz==(int)System->ClusterNZ-1)&&(nz==2))
		{zshift=System->ProcInfo.Box.lz;do_shift=true;}
		else zshift=0.0;

	}

	return do_shift;
}

void CommunicationHandler::UnpackSpace() {
	int newna=System->LocalAtomNumber;
	OMD_FLOAT xshift, yshift, zshift;
	System->AtomStorage.Cut(System->LocalAtomNumber);

	// appending all received atoms from other procs...
	for(int b=0;b<27;b++) {

		if(System->GetNeighborRank(b)<0||b==MYSELF||RecvNumber[b]==0) continue;
		
		if(CheckCellShift(b,xshift,yshift,zshift)) {
			for(int c=0;c<RecvNumber[b];c++){
				SpaceRecvBuffer[b][c].x+=xshift;
				SpaceRecvBuffer[b][c].y+=yshift;
				SpaceRecvBuffer[b][c].z+=zshift;
			}
		}

		PackageSpace* X=SpaceRecvBuffer[b];
		System->AtomStorage.Expand(newna+RecvNumber[b]);

		for(int c=0;c<RecvNumber[b];c++) {
			Atom* A=System->AtomPtr(c+newna);
			A->id=X->id;
			A->xid=X->xid;
			A->nid=X->nid;
			A->flag=X->flag|FLAG_GHOST;
			A->x=X->x;
			A->y=X->y;
			A->z=X->z;
			A->vx=X->vx;
			A->vy=X->vy;
			A->vz=X->vz;
			X++;
		}
		
		newna+=RecvNumber[b];
	} // for all rubix
}

void CommunicationHandler::UnpackPosition() {
	int oldna=System->LocalAtomNumber;
	int newna=System->LocalAtomNumber;
	OMD_FLOAT xshift, yshift, zshift;
	
	for(int b=0;b<27;b++) {
		if(System->GetNeighborRank(b)<0||b==MYSELF||RecvNumber[b]==0) continue;

		if(CheckCellShift(b,xshift,yshift,zshift)){
			for(int c=0;c<RecvNumber[b];c++){
				VectorRecvBuffer[b][c].x+=xshift;
				VectorRecvBuffer[b][c].y+=yshift;
				VectorRecvBuffer[b][c].z+=zshift;
			}
		}
		
		newna+=RecvNumber[b];
		PackageVector* X=VectorRecvBuffer[b];
		for(int c=0;c<RecvNumber[b];c++) {
			Atom* A=System->AtomPtr(c+oldna);
			
			if(A->nid!=X->nid)
				die("unpacking position: atoms structure is non-contiguous! attempt to assign nid="+
					as_string(A->nid)+" with nid="+as_string(X->nid));

			A->x=X->x;
			A->y=X->y;
			A->z=X->z;
			X++;
		}
		oldna=newna;
	}

	if(newna!=GetNAtom()) die("unpacking force: missing/additional atoms total after received="+
							  as_string(newna)+" number of atoms="+as_string(GetNAtom()));
}

void CommunicationHandler::UnpackVelocity() {
	int oldna=System->LocalAtomNumber;
	int newna=System->LocalAtomNumber;
	
	// appending all received atoms from other procs...
	for(int b=0;b<27;b++) {
		
		if(System->GetNeighborRank(b)<0||b==MYSELF||RecvNumber[b]==0) continue;
		
		newna+=RecvNumber[b];
		PackageVector* X=VectorRecvBuffer[b];
		for(int c=0;c<RecvNumber[b];c++) {
			Atom* A=System->AtomPtr(c+oldna);
			
			if(A->nid!=X->nid)
				die("unpacking velocity: atoms structure is non-contiguous! attempt to assign nid="+
					as_string(A->nid)+" with nid="+as_string(X->nid));
			
			A->vx=X->x;
			A->vy=X->y;
			A->vz=X->z;
			X++;
		}
		oldna=newna;
	}
	
	if(newna!=GetNAtom()) die("unpacking force: missing/additional atoms total after received="+
							  as_string(newna)+" number of atoms="+as_string(GetNAtom()));
	
}

void CommunicationHandler::UnpackForce() {

	int oldna=System->LocalAtomNumber;
	int newna=System->LocalAtomNumber;
	
	// appending all received atoms from other procs...
	for(int b=0;b<27;b++) {
		
		if(System->GetNeighborRank(b)<0||b==MYSELF||RecvNumber[b]==0) continue;
		
		newna+=RecvNumber[b];
		PackageVector* X=VectorRecvBuffer[b];
		for(int c=0;c<RecvNumber[b];c++) {
			Atom* A=System->AtomPtr(c+oldna);
			
			if(A->nid!=X->nid)
				die("unpacking force: atoms structure is non-contiguous! attempt to assign nid="+
					as_string(A->nid)+" with nid="+as_string(X->nid));
			
			A->fx=X->x;
			A->fy=X->y;
			A->fz=X->z;
			X++;
		}
		oldna=newna;
	}
	
	if(newna!=GetNAtom()) die("unpacking force: missing/additional atoms total after received="+
							  as_string(newna)+" number of atoms="+as_string(GetNAtom()));
	
}

void CommunicationHandler::UnpackAux(int aidx) {

	int oldna=System->LocalAtomNumber;
	int newna=System->LocalAtomNumber;

	// appending all received atoms from other procs...
	for(int b=0;b<27;b++) {

		if(System->GetNeighborRank(b)<0||b==MYSELF||RecvNumber[b]==0) continue;

		newna+=RecvNumber[b];
		PackageScalar* X=ScalarRecvBuffer[b];
		for(int c=0;c<RecvNumber[b];c++) {
			Atom* A=System->AtomPtr(c+oldna);

			if(A->nid!=X->nid)
				die("unpacking aux["+as_string(aidx)+
						"]: atoms structure is non-contiguous! attempt to assign nid="+
						as_string(A->nid)+" with nid="+as_string(X->nid));

			A->aux[aidx]=X->value;
			X++;
		}

		oldna=newna;
	}

	if(newna!=GetNAtom()) die("unpacking aux["+as_string(aidx)+"] missing/additional atoms total after received="+
			as_string(newna)+" number of atoms="+as_string(GetNAtom()));

}

void CommunicationHandler::SendReceiveData(void* sendpack[], void* recvpack[], int unitlength){
	MPI_Status  stat[54];
	MPI_Request req[54];
	int nreq=0;
	
	for(int i=0;i<27;i++) {
		if((i==MYSELF)||(System->GetNeighborRank(i)<0||RecvNumber[i]<=0))continue;
		MemRealloc(recvpack[i], RecvNumber[i]*unitlength);
	}

	for(int i=0;i<27;i++){
		int proc=System->GetNeighborRank(i);
		if(i==MYSELF||proc<0)continue;
		
		if(SendNumber[i]) MPI_Send_init(sendpack[i],SendNumber[i]*unitlength,
		                       MPI_CHAR,proc,DATATAG+(26-i),
		                       MPI_COMM_WORLD,&req[nreq++]);

		if(RecvNumber[i]) MPI_Recv_init(recvpack[i],RecvNumber[i]*unitlength,
		                       MPI_CHAR,proc,DATATAG+i,
	        	               MPI_COMM_WORLD, &req[nreq++]);

	}

	MPI_Startall(nreq,req);
	MPI_Waitall(nreq,req,stat);
	for(int i=0;i<nreq;i++) MPI_Request_free(&req[i]);

}

void CommunicationHandler::SendReceive(int mode){
	timeval tstart, tend;
	gettimeofday(&tstart,NULL);

	SyncProcesses();
	CollectSendRecvNumber();

	if(mode&SYNC_SPACE) {
		PrepareSpaceBuffers();
		SendReceiveData((void**) SpaceSendBuffer, 
						(void**) SpaceRecvBuffer, sizeof(PackageSpace));
		UnpackSpace();
	}
	
	if(mode&SYNC_POSITION) { 
		PreparePositionBuffers();
		SendReceiveData((void**) VectorSendBuffer, 
						(void**) VectorRecvBuffer, sizeof(PackageVector));
		UnpackPosition();
	}
	
	if(mode&SYNC_VELOCITY) {
		PrepareVelocityBuffers();
		SendReceiveData((void**) VectorSendBuffer, 
						(void**) VectorRecvBuffer, sizeof(PackageVector));
		UnpackVelocity();
	}
		
	if(mode&SYNC_FORCE) {
		PrepareForceBuffers();
		SendReceiveData((void**) VectorSendBuffer, 
						(void**) VectorRecvBuffer, sizeof(PackageVector));
		UnpackForce();
	}
	
	if(mode&SYNC_AUX) {
		int aidx=mode&0x00FF;
		PrepareAuxBuffers(aidx);
		SendReceiveData((void**)   ScalarSendBuffer, 
						(void**)   ScalarRecvBuffer, sizeof(PackageScalar));
		UnpackAux(aidx);
	}

	gettimeofday(&tend,NULL);
	TotalComtime+=(OMD_FLOAT)((tend.tv_sec-tstart.tv_sec)+1e-6*(tend.tv_usec-tstart.tv_usec));
}

int CommunicationHandler::TakeSUM(int a){
	int mv;
	MPI_Allreduce(&a,&mv,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	return mv;
}

int CommunicationHandler::TakeMAX(int a){
	int mv;
	MPI_Allreduce(&a,&mv,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
	return mv;
}

int CommunicationHandler::TakeMIN(int a){
	int mv;
	MPI_Allreduce(&a,&mv,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
	return mv;
}


OMD_FLOAT CommunicationHandler::TakeSUM(OMD_FLOAT a){
	OMD_FLOAT mv;
	MPI_Allreduce(&a,&mv,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	return mv;
}

OMD_FLOAT CommunicationHandler::TakeMAX(OMD_FLOAT a){
	OMD_FLOAT mv;
	MPI_Allreduce(&a,&mv,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
	return mv;
}

OMD_FLOAT CommunicationHandler::TakeMIN(OMD_FLOAT a){
	OMD_FLOAT mv;
	MPI_Allreduce(&a,&mv,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
	return mv;
}

/** Gather variable values from childs. The target buffer must be
 *  allocated prior to the function call **/

void CommunicationHandler::Gather(void *data, void* buffer, int length){
	MPI_Gather(data, length, MPI_CHAR, buffer, length, MPI_CHAR, 0, MPI_COMM_WORLD);
}

void CommunicationHandler::AllGather(void *send, void* receive, int senlen, int reclen){
	MPI_Allgather(send, senlen, MPI_CHAR, receive, reclen, MPI_CHAR, MPI_COMM_WORLD);
}

void CommunicationHandler::Broadcast(void* a, int size){
	MPI_Bcast(a,size,MPI_CHAR,0,MPI_COMM_WORLD);
}

void CommunicationHandler::SyncProcesses() {MPI_Barrier(MPI_COMM_WORLD);}

void CommunicationHandler::Abort() {
	blog("Rank ("+as_string(Rank)+") sending abort message...", LOGINFO);
	MPI_Abort(MPI_COMM_WORLD, System->ExitCode);
}

void CommunicationHandler::DistributeAtomIndexSlab(int slab) {

	int center=slab+4, south=slab+1, north=slab+7, west=slab+3, east=slab+5;
	if(System->GetNeighborRank(center)<0) return;

	OMD_FLOAT rcut=System->GetMaxCutRadius()+CellRadiusTolerance;

	if(rcut>System->ProcInfo.LCellX || rcut>System->ProcInfo.LCellY || rcut>System->ProcInfo.LCellZ)
		die("processor cell size is too small cell=("+
				as_string(System->ProcInfo.LCellX)+","+
				as_string(System->ProcInfo.LCellY)+","+
				as_string(System->ProcInfo.LCellZ)+") border radius="+as_string(rcut));

	int na=System->GetLocalAtomNumber();
	SysBox Border=System->GetCellBorder();

	switch(slab) {
	case BottomSlab:

		for(int i=0;i<na;i++)
			if(Atoms(i).z<(Border.z0+rcut)) NeigRubix[center]->push(i);

		break;

	case MidSlab:

		for(int i=0;i<na;i++) NeigRubix[center]->push(i);

		break;

	case TopSlab:

		for(int i=0;i<na;i++)
			if(Atoms(i).z>(Border.z1-rcut)) NeigRubix[center]->push(i);

		break;
	}

	na=NeigRubix[center]->length;
	int a;

	if(System->GetNeighborRank(west)>=0) {
		NeigRubix[center]->reset();
		while((a=NeigRubix[center]->fetch())>=0) {
			if(Atoms(a).x<(Border.x0+rcut)) NeigRubix[west]->push(a);
		}
	}

	if(System->GetNeighborRank(east)>=0) {
		NeigRubix[center]->reset();
		while((a=NeigRubix[center]->fetch())>=0) {
			if(Atoms(a).x>(Border.x1-rcut)) NeigRubix[east]->push(a);
		}
	}

	if(System->GetNeighborRank(south)>=0) {
		NeigRubix[center]->reset();
		while((a=NeigRubix[center]->fetch())>=0) {
			if(Atoms(a).y<(Border.y0+rcut)) NeigRubix[south]->push(a);
		}
	}

	if(System->GetNeighborRank(north)>=0) {
		NeigRubix[center]->reset();
		while((a=NeigRubix[center]->fetch())>=0) {
			if(Atoms(a).y>(Border.y1-rcut)) NeigRubix[north]->push(a);
		}
	}

	na=NeigRubix[south]->length;

	if(System->GetNeighborRank(south-1)>=0) {
		NeigRubix[south]->reset();
		while((a=NeigRubix[south]->fetch())>=0) {
			if(Atoms(a).x<(Border.x0+rcut)) NeigRubix[south-1]->push(a);
		}
	}

	if(System->GetNeighborRank(south+1)>=0) {
		NeigRubix[south]->reset();
		while((a=NeigRubix[south]->fetch())>=0) {
			if(Atoms(a).x>(Border.x1-rcut)) NeigRubix[south+1]->push(a);
		}
	}

	na=NeigRubix[north]->length;

	if(System->GetNeighborRank(north-1)>=0) {
		NeigRubix[north]->reset();
		while((a=NeigRubix[north]->fetch())>=0) {
			if(Atoms(a).x<(Border.x0+rcut)) NeigRubix[north-1]->push(a);
		}
	}

	if(System->GetNeighborRank(north+1)>=0) {
		NeigRubix[north]->reset();
		while((a=NeigRubix[north]->fetch())>=0) {
			if(Atoms(a).x>(Border.x1-rcut)) NeigRubix[north+1]->push(a);
		}
	}

}


/**
 * @brief Distributes the indices of atoms to be sent to the neighboring process.
 *
 * The CommunicationHandler will assign the atom data to the buffer and exchange
 * to the neighbor procs. The index of atoms is book-keeped in IndexLists.
 */

void CommunicationHandler::DistributeAtomIndex() {
	for(int i=0;i<27;i++) NeigRubix[i]->clear();
	DistributeAtomIndexSlab(MidSlab);
	if(System->GetNeighborRank(BottomSlab+4)>=0)
		DistributeAtomIndexSlab(BottomSlab);
	if(System->GetNeighborRank(TopSlab+4)>=0)
		DistributeAtomIndexSlab(TopSlab);
}


Atom& CommunicationHandler::Atoms(int idx){return System->Atoms(idx);}
Atom* CommunicationHandler::AtomPtr(int idx){return System->AtomPtr(idx);}
int CommunicationHandler::GetNAtom(){return System->GetNAtom();}
