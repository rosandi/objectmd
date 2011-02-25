/* omdlib
 ********************************************************
 *
 * ObjectMD
 * Molecular Dynamics Class Library
 * Version 2.0 (2009)
 *
 * Released under:
 *       GNU GENERAL PUBLIC LICENSE Version 3
 *       Copyright (C) 2005  Yudi Rosandi
 *
 * Read license.txt file in the source's root directory.
 *
 ********************************************************
 * Implementation: MDSystemGrid class
 *
*/

#include <iostream>
#include <iomanip>
#include <sstream>
#include <omd/gadget.hpp>
#include <omd/systemgrid.hpp>
#include <omd/comhandler.hpp>
#include <omd/integrator.hpp>
#include <omd/iterator.hpp>

using std::ostringstream;

// ProcInfo.Box and Box <- keeps the simulation box
// Border <- keeps the processor borders (ThisProcAtoms are within this Border)

MDSystemGrid::MDSystemGrid(OMD_INT &argc,OMD_CHAR** &argv,OMD_INT nx,OMD_INT ny,OMD_INT nz)
              :MDSystem(argc,argv){
	set_name("SIMULATION_SYSTEM_GRID");
    register_class(get_name());
	// -1 == no neighbor
	for(OMD_INT i=0;i<27;i++){
		for(OMD_INT j=0;j<MAXPROC;j++) ProcInfo.NeighborList[j][i]=-1;
	}
	ClusterNX=nx;
	ClusterNY=ny;
	ClusterNZ=nz;
	CommRefreshPeriod=-1;
	Communicator=NULL;
	FirstSync=true;
	BinDirectory="pomd";
	LocalBuffer.set_logger(logger);
	LocalBuffer.set_name("localbuffer");
	GhostBuffer.set_logger(logger);
	GhostBuffer.set_name("ghostbuffer");
	LocalAtomNumber=0;
}

MDSystemGrid::MDSystemGrid(){
	set_name("SIMULATION_SYSTEM_GRID");
    register_class(get_name());
	// -1 == no neighbor
	for(OMD_INT i=0;i<27;i++){
		for(OMD_INT j=0;j<MAXPROC;j++) ProcInfo.NeighborList[j][i]=-1;
	}
	ClusterNX=ClusterNY=ClusterNZ=1;
	CommRefreshPeriod=-1;
	Communicator=NULL;
	FirstSync=true;
	BinDirectory="pomd";
	LocalBuffer.set_logger(logger);
	LocalBuffer.set_name("localbuffer");
	GhostBuffer.set_logger(logger);
	GhostBuffer.set_name("ghostbuffer");
	LocalAtomNumber=0;
}

MDSystemGrid::~MDSystemGrid() {
	delete Communicator;
}

void MDSystemGrid::SetBinDirectory(string bindir){BinDirectory.assign(bindir);}

/** 
 * @brief Get the processor number from a coordinate
 * 
 * Check to which processor an atom belongs to.
 * Return the rank of processor that ownes the atom.
 */

OMD_INT MDSystemGrid::GetCellIndex(OMD_FLOAT x, OMD_FLOAT y, OMD_FLOAT z) {
	OMD_INT u=(OMD_INT)floor((x-ProcInfo.Box.x0)/ProcInfo.LCellX);
	OMD_INT v=(OMD_INT)floor((y-ProcInfo.Box.y0)/ProcInfo.LCellY);
	OMD_INT w=(OMD_INT)floor((z-ProcInfo.Box.z0)/ProcInfo.LCellZ);
	
	if(u<0)u=IS_PERIODIC_X?ClusterNX-1:0;
	if(u>=(OMD_INT)ClusterNX)u=IS_PERIODIC_X?0:ClusterNX-1;

	if(v<0)v=IS_PERIODIC_Y?ClusterNY-1:0;
	if(v>=(OMD_INT)ClusterNY)v=IS_PERIODIC_Y?0:ClusterNY-1;

	if(w<0)w=IS_PERIODIC_Z?ClusterNZ-1:0;
	if(w>=(OMD_INT)ClusterNZ)w=IS_PERIODIC_Z?0:ClusterNZ-1;

	OMD_SIZET np=w*(ClusterNX*ClusterNY)+v*ClusterNX+u;

	if(np>=GetGridSize()){
		die("atom lying outside simulation box, probably rounding error. FIX SIMULATION BOX! ("+
			as_string(x,"%0.5f")+","+as_string(y,"%0.5f")+","+as_string(z,"%0.5f")+") box: ("+
			as_string(Box.x0)+","+as_string(Box.y0)+","+as_string(Box.z0)+") ("+
			as_string(Box.x1)+","+as_string(Box.y1)+","+as_string(Box.z1)+")"
			);
	}
	
	return np;
}

/**
 * @brief check ownership of an atom
 * 
 * The function gives the ownership of an atom position to the current cell.
 * Atoms which are not in this processor box are marked as ghosts. The atoms
 * owned by the processor are stored in ThisProcAtoms.
 * 
 */

bool MDSystemGrid::CheckOwnership(OMD_FLOAT x, OMD_FLOAT y, OMD_FLOAT z) {
	OMD_INT rank=GetRank();

	OMD_INT u=(OMD_INT)floor((x-ProcInfo.Box.x0)/ProcInfo.LCellX);
	OMD_INT v=(OMD_INT)floor((y-ProcInfo.Box.y0)/ProcInfo.LCellY);
	OMD_INT w=(OMD_INT)floor((z-ProcInfo.Box.z0)/ProcInfo.LCellZ);
	
	if(IS_PERIODIC_X){if(u<0 || u>=(OMD_INT)ClusterNX)return false;}
	if(IS_PERIODIC_Y){if(v<0 || v>=(OMD_INT)ClusterNY)return false;}
	if(IS_PERIODIC_Z){if(w<0 || w>=(OMD_INT)ClusterNZ)return false;}

	if(u<0)u=0;
	if(v<0)v=0;
	if(w<0)w=0;
	if(u>=(OMD_INT)ClusterNX)u=ClusterNX-1;
	if(v>=(OMD_INT)ClusterNY)v=ClusterNY-1;
	if(w>=(OMD_INT)ClusterNZ)w=ClusterNZ-1;

	if(u==(OMD_INT)ProcInfo.CellX[rank]&&
	   v==(OMD_INT)ProcInfo.CellY[rank]&&
	   w==(OMD_INT)ProcInfo.CellZ[rank]) return true;

	return false;
}

/**
 * @brief Arranges the processors to the Cluster architecture
 * 
 * The architecture can be linear, planar and cube.
 * (check)
 */

void MDSystemGrid::Root_ArrangeNeighbor() {
	OMD_FLOAT plx=Box.lx/(OMD_FLOAT)ClusterNX;
	OMD_FLOAT ply=Box.ly/(OMD_FLOAT)ClusterNY;
	OMD_FLOAT plz=Box.lz/(OMD_FLOAT)ClusterNZ;
	
	ProcInfo.LCellX=plx;
	ProcInfo.LCellY=ply;
	ProcInfo.LCellZ=plz;

	for(OMD_INT i=0;i<MAXPROC;i++) for(OMD_INT j=0;j<27;j++) 
		ProcInfo.NeighborList[i][j]=-1;
		
	OMD_INT nproc=GetGridSize();
	for(OMD_INT p=0;p<nproc;p++) {

		OMD_INT px=p%ClusterNX;
		OMD_INT py=(p/ClusterNX)%ClusterNY;
		OMD_INT pz=(p/(ClusterNX*ClusterNY))%ClusterNZ;
		
		ProcInfo.CellX[p]=px;
		ProcInfo.CellY[p]=py;
		ProcInfo.CellZ[p]=pz;
		ProcInfo.Border[p].x0=Box.x0+(OMD_FLOAT)px*plx;
		ProcInfo.Border[p].x1=Box.x0+(OMD_FLOAT)(px+1)*plx;
		ProcInfo.Border[p].y0=Box.y0+(OMD_FLOAT)py*ply;
		ProcInfo.Border[p].y1=Box.y0+(OMD_FLOAT)(py+1)*ply;
		
		ProcInfo.Border[p].z0=Box.z0+(OMD_FLOAT)pz*plz;
		ProcInfo.Border[p].z1=Box.z0+(OMD_FLOAT)(pz+1)*plz;
		
		//FIXME! box size must be checked_> compared with Max R_Cut
		//       especially when running on single proc using periodic
		//       boundary condition.
		
		OMD_INT n=0;
		for(OMD_INT k=-1;k<=1;k++) 
			for(OMD_INT j=-1;j<=1;j++)
				for(OMD_INT i=-1;i<=1;i++) {
					// n? => neighbor
					OMD_INT nx=px+i;
					OMD_INT ny=py+j;
					OMD_INT nz=pz+k;

					if(IS_PERIODIC_X) {
						if(nx<0)nx=ClusterNX-1;
						if(nx>=(OMD_INT)ClusterNX)nx=0;
					} else {
						if((nx<0)||(nx>=(OMD_INT)ClusterNX)){n++;continue;}
					}
					
					if(IS_PERIODIC_Y) {
						if(ny<0)ny=ClusterNY-1;
						if(ny>=(OMD_INT)ClusterNY)ny=0;
					} else {
						if((ny<0)||(ny>=(OMD_INT)ClusterNY)){n++;continue;}
					}

					if(IS_PERIODIC_Z) {
						if(nz<0)nz=ClusterNZ-1;
						if(nz>=(OMD_INT)ClusterNZ)nz=0;
					} else {
						if((nz<0)||(nz>=(OMD_INT)ClusterNZ)){n++;continue;}
					}

					ProcInfo.NeighborList[p][n++]= 
						nz*(ClusterNX*ClusterNY)+ny*ClusterNX+nx;
				}
	}
}

/**
 * @brief distribute the atoms to procs
 * 
 * 
 * The atoms are devided in files.
 * Each processor will read the file having name: p<proc-num>,
 * under "BinDirectory" directory.
 */

void MDSystemGrid::Root_DistributeAtoms() {
	
	#define NTOLE 100

	AtomContainer* ProcAtm[GetGridSize()][SystemAtoms.size()];

	for(OMD_SIZET i=0;i<SystemAtoms.size();i++){
		AtomKeeper Atm(SystemAtoms[i]->GetAtomStorage());
		OMD_INT nproc=GetGridSize();
		OMD_INT npp=Atm.GetNAtom()/nproc+NTOLE;

		for(OMD_INT a=0;a<nproc;a++){ // same atom type for all procs
			ProcAtm[a][i]=new AtomContainer;
			ProcAtm[a][i]->disable_log();
			ProcAtm[a][i]->copy_data(SystemAtoms[i]);
			ProcAtm[a][i]->Allocate(npp,false,AtomKeeper::Referral);
		}

		for(OMD_SIZET a=0;a<Atm.GetNAtom();a++) {
			OMD_INT pn=GetCellIndex(Atm[a].x,Atm[a].y,Atm[a].z);
			ProcAtm[pn][i]->GetAtomStorage().Attach(Atm[a]);
		}
	}

	// one file/proc
	OMD_SIZET cntot=0;
	for(OMD_SIZET pn=0;pn<GetGridSize();pn++) {
		string fname(BinDirectory+"/p"+as_string(pn)+"-crystal");
		OMD_SIZET cnt=0;
		for(OMD_SIZET a=0;a<SystemAtoms.size();a++){
			ProcAtm[pn][a]->Save(fname, "a");
			cnt+=ProcAtm[pn][a]->GetNAtom();
		}
		cntot+=cnt;
	}
	
	for(OMD_SIZET i=0;i<SystemAtoms.size();i++) {
		strncpy(ProcInfo.AtomNames[i],SystemAtoms[i]->get_name().c_str(),32);
		ProcInfo.AtomNames[i][31]=0x0;
	}
	
	for(OMD_SIZET a=0;a<GetGridSize();a++){
		for(OMD_SIZET i=0;i<SystemAtoms.size();i++) {
			ProcAtm[a][i]->Release();
			delete ProcAtm[a][i];
		}
	}
}

/**
 * @brief Prepare information and data
 * 
 * Construct the SystemAtoms and distribute the atoms.
 * The SystemAtoms is cleared after distributing
 */

void MDSystemGrid::Root_Prepare() {
	// In this stage Box contains simulation Box
	UnificateAtoms();
	EnumerateAtoms();
	MDSystem::AdjustSystem();	
	ProcInfo.Box=Box; // This is the valid system box after adjusted by user!
	ProcInfo.TotalAtom=GetNAtom();
	ProcInfo.NumberOfContainers=SystemAtoms.size();
	Root_ArrangeNeighbor();
	Root_DistributeAtoms();
	for(OMD_SIZET a=0;a<SystemAtoms.size();a++) {
		SystemAtoms[a]->Release();
		delete SystemAtoms[a];
	}
	AtomID=0;
	AtomStorage.Release();
	SystemAtoms.clear();
}

bool MDSystemGrid::CheckComm(){
	if(Communicator) {
		return Communicator->CheckLink();
	}
	return false;
}

bool MDSystemGrid::CheckArch() {
	string clusdim(as_string(ClusterNX)+"x"+as_string(ClusterNY)+"x"+as_string(ClusterNZ));

	if(ClusterNX<0) {
		switch(ClusterNX){
			case GRID_AUTOX:
				ClusterNX=GetGridSize();
				ClusterNY=1;
				ClusterNZ=1;
				break;
			case GRID_AUTOY:
				ClusterNX=1;
				ClusterNY=GetGridSize();
				ClusterNZ=1;
				break;
			case GRID_AUTOZ:
				ClusterNX=1;
				ClusterNY=1;
				ClusterNZ=GetGridSize();
				break;
			default:
				die("wrong architecture option");				
		}
	}

	assert(ClusterNX*ClusterNY*ClusterNZ==(OMD_INT)GetGridSize(),
	  "incompatible grid dimension: "+clusdim+
	  " on "+as_string(GetGridSize())+" processor(s)");

	assert(GetGridSize()<=MAXPROC,
	  "this program is compiled with "+as_string(MAXPROC)+" maximum processors");

	blog("grid size "+clusdim);
	  
	return true;
}

/**
 * @brief Load this processor's simulation atoms
 * 
 * search in "BinDirectory" directory
 * the structure name: p[procnum] 
 * e.g. p1
 */

void MDSystemGrid::LoadAtoms() {
	string fname(BinDirectory+"/p"+as_string(GetRank())+"-crystal");
	for(OMD_SIZET i=0;i<ProcInfo.NumberOfContainers;i++) {
		AtomContainer* A=new AtomContainer;
		A->set_logger(logger);
		A->Load(fname, ProcInfo.AtomNames[i]);
		AddAtom(A);
	}
	// don't re-enumerate!
	Enumerated=true;
}

void MDSystemGrid::CommInit(){	
	if(!Communicator){Communicator=new CommunicationHandler;}
	if(!Communicator->CheckLink()){
		Communicator->Link(this);

		sendlog("log redirected to "+BinDirectory);
		sendmemlog("memory log redirected to "+BinDirectory);
		blogstream.close();
		memlogstream.close();

  		if(GetRank()==0) {
			system(("rm -rf "+BinDirectory).c_str());
			system(("mkdir -p "+BinDirectory).c_str());
		}

		Communicator->SyncProcesses();
		memlogstream.open((BinDirectory+"/p"+as_string(GetRank())+"-mem.log").c_str());
		blogstream.open((BinDirectory+"/p"+as_string(GetRank())+"-run.log").c_str());
	}
}

/** synchronize parameters... **/
void MDSystemGrid::SyncEnvironment() {		
		OMD_CHAR parstr[4096];
		if(GetRank()==0){
			string ts;
			for(OMD_INT i=0;i<param.size();i++){
				ts.append(param[i]+" ");
			}
			memset(parstr, 0, 4096);
			ts.copy(parstr, 4095);
		}
		Communicator->Broadcast(parstr, 4096);
		param.assign(parstr);
}

void MDSystemGrid::SyncVariable() {
	Communicator->Broadcast(&SimBeginTime, sizeof(time_t));
	Communicator->Broadcast(&Step, sizeof(OMD_INT));
	Communicator->Broadcast(&PBoundary, sizeof(OMD_INT));
	Communicator->Broadcast(&Energy, sizeof(OMD_FLOAT));
	Communicator->Broadcast(&Kinetic, sizeof(OMD_FLOAT));
	Communicator->Broadcast(&Virial, sizeof(OMD_FLOAT));
	Communicator->Broadcast(&Potential, sizeof(OMD_FLOAT));
	Communicator->Broadcast(&BasePotential, sizeof(OMD_FLOAT));
	Communicator->Broadcast(&ElapsedTime, sizeof(OMD_FLOAT));
	Communicator->Broadcast(&Box, sizeof(SysBox));
	Communicator->Broadcast(&write_mode, sizeof(OMD_INT));
	
	OMD_INT varsize=RestartVars.size();
	Communicator->Broadcast(&varsize, sizeof(OMD_INT));

	for(OMD_INT i=0;i<varsize;i++) {
		struct {
			OMD_CHAR l[128];OMD_CHAR v[128];
		} rvar;
		
		if(GetRank()==ROOT) {
			RestartVars[i]->GetLabel().copy(rvar.l,128);
			RestartVars[i]->AsString().copy(rvar.v,128);
		}
		
		Communicator->Broadcast(&rvar, sizeof(rvar));
		
		if(GetRank()!=ROOT) {
			DataSlot *m=new DataSlot(rvar.l);
			m->SetDefaultData(rvar.v);
			RestartVars.push_back(m);
		}
	}

}

void MDSystemGrid::CreationFunction() {
// Here: conditional: restarting or normal.....
	if(GetRank()==0) {
		if(Mode==NORMAL_MODE)CreateSystem();
		if(Mode==RESTART_MODE)LoadSimulation();
		Root_Prepare();
	}

	SyncEnvironment();
	SyncVariable();

	Communicator->Broadcast(&ProcInfo,sizeof(StructInfo));

	LoadAtoms(); // all proc loads their own atoms, ThisProcAtomNumber is updated...

	CreateGadget();
	SyncEnvironment();
}

void MDSystemGrid::FlattenAtomBox() {
	OMD_SIZET na=GetNAtom();
	LocalBuffer.IndexBook.clear();
	GhostBuffer.IndexBook.clear();

	// scan atoms belong to this proc...
	for(OMD_SIZET i=0;i<na;i++) {
		Atom* a=AtomPtr(i);
		if(CheckOwnership(a->x, a->y, a->z)){
			LocalBuffer.IndexBook.push(i);
			a->flag&=~FLAG_GHOST; //clear ghost flag
		} else {
			GhostBuffer.IndexBook.push(i);
			a->flag|=FLAG_GHOST;  // set ghost flag
		}

	}

	LocalBuffer.AssignByIndex(AtomStorage);
	GhostBuffer.AssignByIndex(AtomStorage);

	AtomStorage.Clear();

	// head: local, tail: ghost
	AtomStorage.Append(LocalBuffer);
	AtomStorage.Append(GhostBuffer);
	LocalAtomNumber=LocalBuffer.GetNAtom();
}

void MDSystemGrid::DistributeContainers() {
	OMD_INT nc=SystemAtoms.size();

	for (OMD_INT i=0;i<nc;i++) SystemAtoms[i]->GetAtomStorage().IndexBook.clear();

	for(OMD_INT i=0;i<LocalAtomNumber;i++) {
		Atom* a=AtomPtr(i);

		if(a->id>=nc) die("undefined atom found! index="+as_string(i)+
				" id="+as_string((OMD_INT)a->id)+
				" nid="+as_string((OMD_INT)a->nid)+
				" xid="+as_string((OMD_INT)a->xid));

		SystemAtoms[(OMD_INT)(a->id)]->GetAtomStorage().IndexBook.push(i);
	}

	OMD_INT cntatom=0;
	for(OMD_INT i=0;i<nc;i++){
		SystemAtoms[i]->GetAtomStorage().AssignByIndex(AtomStorage);
		cntatom+=SystemAtoms[i]->GetNAtom();
	}

	// checking
	cntatom=Communicator->TakeSUM(cntatom);

	if(cntatom!=(OMD_INT)ProcInfo.TotalAtom) {
		DumpAtoms(LocalBuffer, "dump-err-local-"+as_string(GetRank()),WM_GHOST|WM_ID|WM_XID|WM_NID);
		DumpAtoms(GhostBuffer, "dump-err-ghost-"+as_string(GetRank()),WM_GHOST|WM_ID|WM_XID|WM_NID);
		Communicator->SyncProcesses();
		die( "missing atoms "+as_string(cntatom)+" total="+as_string(ProcInfo.TotalAtom)+
	       ". to check: boundary conditions, box size, box offset");
	}
}

/**
 * @brief Simulation synchronization
 * 
 * This function synchronizes atom data.
 *
 */

void MDSystemGrid::SyncData(OMD_INT syncmode) {
    if(syncmode&SYNC_SPACE) {
        if(Step%CommRefreshPeriod) {
            Communicator->SendReceive(syncmode);
        } else {
            Communicator->SendReceive(syncmode);
            FlattenAtomBox();
            DistributeContainers();
            UpdateRadiusTolerance();
            Communicator->DistributeAtomIndex();
            Communicator->SendReceive(syncmode);
            Iterator->SetDirty();
        }
    } else {
        Communicator->SendReceive(syncmode);
    }
}

bool MDSystemGrid::CheckRun() {
	InterruptFlag=(OMD_SIZET)Communicator->TakeMAX((OMD_INT)InterruptFlag);
	return MDSystem::CheckRun();
}

void MDSystemGrid::AdjustSystem() { // MDSystem::AdjustSystem is invoked from Root_Prepare
	
	SystemSetting();
	if(PBoundary<0) PBoundary=0;
	
	Box=ProcInfo.Box;
	TotalAtom=ProcInfo.TotalAtom;

	// assign ThisProcAtoms
	OMD_INT na=GetNAtom();
	SqrMaxVelocity=0.0;
	for(OMD_INT i=0;i<na;i++) {
		Atom* a=AtomPtr(i);
		OMD_FLOAT vv=a->vx*a->vx+a->vy*a->vy+a->vz*a->vz;
		if(SqrMaxVelocity<vv) SqrMaxVelocity=vv;
	}

	LocalAtomNumber=na;

}

void MDSystemGrid::UpdateRadiusTolerance() {
	OMD_FLOAT L;
	if(CommRefreshPeriod<1) {L=0.0; CommRefreshPeriod=1;}
	else L=sqrt(SqrMaxVelocity)*(CommRefreshPeriod)*Integrator->TimeStep;
	L+=0.1*Integrator->MaxCutRadius; // add 10% cut radius
	Communicator->SetRadiusTolerance(L);
	Iterator->SetRadiusTolerance(L);
}

void MDSystemGrid::InitGadgets() {
	MDSystem::InitGadgets();
	OMD_FLOAT L;

	UpdateRadiusTolerance();
	Iterator->SetDirty();

	// FirstSync... get ghosts from neighbors...
	Communicator->DistributeAtomIndex();
	Communicator->SendReceive(SYNC_ALL);
	FlattenAtomBox();
}

void MDSystemGrid::FirstRun() {
	MDSystem::FirstRun();
	Communicator->SyncProcesses();
}

void MDSystemGrid::ErrorHandler(const OMD_CHAR* errst) {
	MDSystem::ErrorHandler(errst);
	Communicator->Abort();
}

string MDSystemGrid::GetGridConfiguration() {
	ostringstream ost;

	ost<<"top (26..18):   ";
	for(OMD_INT i=26;i>=18;i--){
		if((i+1)%3==0)ost<<"\n";
		if(GetNeighborRank(i)>=0) ost << "["<<GetNeighborRank(i) << "]";
		else ost<<"[x]";
		
	}

	ost<<"\nmiddle (17..9): ";
	for(OMD_INT i=17;i>=9;i--){
		if((i+1)%3==0)ost<<"\n";
		if(GetNeighborRank(i)>=0) ost << "["<<GetNeighborRank(i) << "]";
		else ost<<"[x]";
	}

	ost<< "\nbottom (8..0):  ";
	for(OMD_INT i=8;i>=0;i--){
		if((i+1)%3==0)ost<<"\n";
		if(GetNeighborRank(i)>=0) ost << "["<<GetNeighborRank(i) << "]";
		else ost<<"[x]";
	}

	return ost.str();	
}


void MDSystemGrid::PrintInfo(ostream& ost) {

	OMD_CHAR buf[DEFAULT_TRANSFER_LENGTH];
	
	if(GetRank()==0) {
		PrintHeader(ost);
		PrintGadgetInfo(ost);
		PrintSystemInfo(ost);

		ost
			<< "\n\n*** Grid informations ***\n"
			<< "Computer grid size " << GetGridSize() << "\n"
			<< "Grid geometry (nx,ny,nz): ("
			<< ClusterNX << ", " << ClusterNY << ", " << ClusterNZ << ")\n"
			<< "Cell dimension (x,y,z): (" 
			<< ProcInfo.LCellX << ", "
			<< ProcInfo.LCellY << ", " 
			<< ProcInfo.LCellZ << ")\n"
			<< "Communication refresh period: " << CommRefreshPeriod << " steps\n\n"
			<< "*** Rank(0)@"<< Communicator->GetHostname() << " ***\n"
			<< "Cell borders " 
			<< "X(" <<GetCellBorder().x0<<", "<<GetCellBorder().x1<<") "
			<< "Y(" <<GetCellBorder().y0<<", "<<GetCellBorder().y1<<") "
			<< "Z(" <<GetCellBorder().z0<<", "<<GetCellBorder().z1<<")\n"
			<< "Number of atoms:\n"
			<< "  total "<<GetNAtom()<<"\n"
			<< "  belong to this cell " <<GetLocalAtomNumber()<<"\n"
			<< "  ghost atoms " << GetNAtom()-GetLocalAtomNumber()<<"\n";

			PrintContainerInfo(ost);

		ost << "\n*** Neighbor processors rank map ***\n" << GetGridConfiguration() <<"\n\n";
		
		for(OMD_SIZET i=1;i<GetGridSize();i++) {
			Communicator->RawReceive(i, buf);
			ost <<buf<<"\n";
		}
	} else {
		ostringstream outst;
		outst << std::fixed << std::setprecision(2)
			  << "*** Rank("<<GetRank()<<")@"<< Communicator->GetHostname() << " ***\n"
			  << "Cell borders " 
			  << "X(" <<GetCellBorder().x0<<", "<<GetCellBorder().x1<<") "
			  << "Y(" <<GetCellBorder().y0<<", "<<GetCellBorder().y1<<") "
			  << "Z(" <<GetCellBorder().z0<<", "<<GetCellBorder().z1<<")\n"
			  << "Number of atoms:\n"
			  << "  total "<<GetNAtom()<<"\n"
			  << "  belong to this cell: " <<GetLocalAtomNumber()<<"\n"
			  << "  ghost atoms: " << GetNAtom()-GetLocalAtomNumber()<<"\n";
		PrintContainerInfo(outst);
		outst << "\n*** Neighbor processors rank map ***\n" << GetGridConfiguration() <<"\n\n";
		strcpy(buf, outst.str().c_str());
		Communicator->RawSend(ROOT,buf,strlen(buf)+1);
	}
	ost << std::scientific << std::setprecision(4);	
}

void MDSystemGrid::PrintInfo(string fname){
	if(GetRank()==ROOT) {
		ofstream fout(fname.c_str());
		PrintInfo(fout);
		fout.close();
	} else PrintInfo(std::cerr); // non-root will not print anything...
}

void MDSystemGrid::MeasureKinetic(){
	MDSystem::MeasureKinetic();	
	Kinetic=Communicator->TakeSUM(Kinetic);
	SqrMaxVelocity=Communicator->TakeMAX(SqrMaxVelocity);
}

void MDSystemGrid::MeasurePotential(){
	MDSystem::MeasurePotential();
	Potential=Communicator->TakeSUM(Potential);
	Virial=Communicator->TakeSUM(Virial);
}

/** collect an atom keeper from all processor, and dump the content.
 * This function does semaphore to avoid overlapping write. 
 */

/* see r368 for mpi data collecting method.*/

AtomContainer* MDSystemGrid::DumpAtoms(
                             AtomKeeper& ak,
                             string fname,
	                         OMD_INT mode,
	                         bool* AuxPrintable,
	                         OMD_CHAR* AuxFormat[],
	                         string AuxNames)
{
	OMD_SIZET myrank=GetRank();
	OMD_CHAR ready=1;
	
	if(myrank==0) {
		AtomContainer::DumpAtoms(ak, fname, mode, AuxPrintable, AuxFormat, AuxNames);
	} else {
		OMD_INT cm;
		if(mode==0)cm=write_mode; else cm=mode;
		cm|=(WM_APPEND|WM_BARE);
		Communicator->RawReceive(myrank-1,&ready,1);
		AtomContainer::DumpAtoms(ak, fname, cm, AuxPrintable, AuxFormat, AuxNames);
	}

	if((myrank+1)<GetGridSize()){
		Communicator->RawSend(myrank+1,&ready,1);
	}

	return this;
}


/** Collects all atoms from all processors. Write the whole simulation box to a file **/

AtomContainer* MDSystemGrid::DumpAtoms(
                             string fname,
	                         OMD_INT mode,
	                         bool* AuxPrintable,
	                         OMD_CHAR* AuxFormat[],
	                         string AuxNames)
{	
	return DumpAtoms(AtomStorage,fname,mode,AuxPrintable,AuxFormat, AuxNames);
} 

void MDSystemGrid::Initiate() {
	CommInit();
	CheckArch();
	SyncEnvironment();
	MDSystem::Initiate();
}

void MDSystemGrid::LoadVariables(FILE* fl) {
	MDSystem::LoadVariables(fl);
}

AtomContainer* MDSystemGrid::Save(string binname, string mode) {
	// collect data
	OMD_SIZET ncon=ProcInfo.NumberOfContainers;
	for(OMD_SIZET i=0;i<ncon;i++) {
		AtomContainer ac;

		OMD_SIZET sz=SystemAtoms[i]->GetNAtom();
		OMD_SIZET asz[GetGridSize()];
		
		Communicator->Gather((void*)&sz, (void*)asz, sizeof(OMD_SIZET));

		ac.copy_data(SystemAtoms[i]);
		ac.GetAtomStorage().set_name("SAVER BUFFER");
		ac.GetAtomStorage().Copy(SystemAtoms[i]->GetAtomStorage());

		Communicator->SyncProcesses();		
		if(GetRank()==ROOT) {
			for(OMD_SIZET p=1;p<GetGridSize();p++) {
				AtomKeeper ak;
				ak.Allocate(asz[p]);
				if(asz[p]>0) {
					Communicator->RawReceive(p, ak.GetArrayPtr(), asz[p]*sizeof(Atom));
					ac.GetAtomStorage().Append(ak);
				}
			}
			ac.Save(binname.c_str(), "a");
		} else {
			if(sz>0)
				Communicator->RawSend(ROOT, ac.GetAtomStorage().GetArrayPtr(), sz*sizeof(Atom));
		}
	}
	
	Communicator->SyncProcesses();
	return this;
}

void MDSystemGrid::SaveSimulation(string binfile) {
	if(binfile=="")binfile.assign(GetRestartFilename());
	if(GetRank()==ROOT) SaveSimulationConfig(binfile);
	Save(binfile);
}

void MDSystemGrid::PrintMessages(ostream& ost) {
	if(GetRank()==ROOT) MDSystem::PrintMessages(ost);
}

void MDSystemGrid::ReadParameters() {
	MDSystem::ReadParameters();
	param.peek("comm.refresh", CommRefreshPeriod);
	if(param.exist("comm.arch")) {
		OMD_INT nx, ny, nz;
		nx=param.int_value("comm.arch",0);
		ny=param.int_value("comm.arch",1);
		nz=param.int_value("comm.arch",2);
		SetClusterArch(nx,ny,nz);
	}
	param.peek("dir.binary", BinDirectory);
}

