/* omdlib
 ********************************************************
 *
 * ObjectMD
 * Molecular Dynamics Class Library
 * Version 3.0 (2009, 2011)
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
#include <omd/omdtool.h>
#include <omd/gadget.h>
#include <omd/systemgrid.h>
#include <omd/comhandler.h>
#include <omd/integrator.h>
#include <omd/iterator.h>

using std::ostringstream;
using namespace omd;

// ProcInfo.Box and Box <- keeps the simulation box
// Border <- keeps the processor borders (ThisProcAtoms are within this Border)

MDSystemGrid::MDSystemGrid(int &argc,char** &argv,int nx,int ny,int nz)
:MDSystem(argc,argv){
	set_name("SIMULATION_SYSTEM_GRID");
  register_class(get_name());
	// -1 == no neighbor
	for(int i=0;i<27;i++){
		for(int j=0;j<MAXPROC;j++) ProcInfo.NeighborList[j][i]=-1;
	}
	ClusterNX=nx;
	ClusterNY=ny;
	ClusterNZ=nz;
	CommRefreshPeriod=-1;
  CommRadTol=0.1;  // 10% of max. cut radius
  CommAbsoluteTol=false;
	Communicator=NULL;
	FirstSync=true;
	BinDirectory="pomd";
	LocalBuffer.set_logger(logger);
	LocalBuffer.set_name("LOCAL BUFFER");
	GhostBuffer.set_logger(logger);
	GhostBuffer.set_name("GHOST BUFFER");
	LocalAtomNumber=0;
}

MDSystemGrid::MDSystemGrid(){
	set_name("SIMULATION_SYSTEM_GRID");
  register_class(get_name());
	
  // -1 == no neighbor
	for(int i=0;i<27;i++){
		for(int j=0;j<MAXPROC;j++) ProcInfo.NeighborList[j][i]=-1;
	}
  
	ClusterNX=ClusterNY=ClusterNZ=1;
	CommRefreshPeriod=-1;
	Communicator=NULL;
	FirstSync=true;
	BinDirectory="pomd";
	LocalBuffer.set_logger(logger);
	LocalBuffer.set_name("LOCAL BUFFER");
	GhostBuffer.set_logger(logger);
	GhostBuffer.set_name("GHOST BUFFER");
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

int MDSystemGrid::GetCellIndex(OMD_FLOAT x, OMD_FLOAT y, OMD_FLOAT z) {
	int u=(int)floor((x-ProcInfo.Box.x0)/ProcInfo.LCellX);
	int v=(int)floor((y-ProcInfo.Box.y0)/ProcInfo.LCellY);
	int w=(int)floor((z-ProcInfo.Box.z0)/ProcInfo.LCellZ);
	
	if(u<0)u=IS_PERIODIC_X?ClusterNX-1:0;
	if(u>=(int)ClusterNX)u=IS_PERIODIC_X?0:ClusterNX-1;

	if(v<0)v=IS_PERIODIC_Y?ClusterNY-1:0;
	if(v>=(int)ClusterNY)v=IS_PERIODIC_Y?0:ClusterNY-1;

	if(w<0)w=IS_PERIODIC_Z?ClusterNZ-1:0;
	if(w>=(int)ClusterNZ)w=IS_PERIODIC_Z?0:ClusterNZ-1;

	int np=w*(ClusterNX*ClusterNY)+v*ClusterNX+u;

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
	int rank=GetRank();

	int u=(int)floor((x-ProcInfo.Box.x0)/ProcInfo.LCellX);
	int v=(int)floor((y-ProcInfo.Box.y0)/ProcInfo.LCellY);
	int w=(int)floor((z-ProcInfo.Box.z0)/ProcInfo.LCellZ);
	
	if(IS_PERIODIC_X){if(u<0 || u>=(int)ClusterNX)return false;}
	if(IS_PERIODIC_Y){if(v<0 || v>=(int)ClusterNY)return false;}
	if(IS_PERIODIC_Z){if(w<0 || w>=(int)ClusterNZ)return false;}

	if(u<0)u=0;
	if(v<0)v=0;
	if(w<0)w=0;
	if(u>=(int)ClusterNX)u=ClusterNX-1;
	if(v>=(int)ClusterNY)v=ClusterNY-1;
	if(w>=(int)ClusterNZ)w=ClusterNZ-1;

	if(u==(int)ProcInfo.CellX[rank]&&
	   v==(int)ProcInfo.CellY[rank]&&
	   w==(int)ProcInfo.CellZ[rank]) return true;

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

	for(int i=0;i<MAXPROC;i++) for(int j=0;j<27;j++) 
		ProcInfo.NeighborList[i][j]=-1;
		
	int nproc=GetGridSize();
	for(int p=0;p<nproc;p++) {

		int px=p%ClusterNX;
		int py=(p/ClusterNX)%ClusterNY;
		int pz=(p/(ClusterNX*ClusterNY))%ClusterNZ;
		
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
		
		int n=0;
		for(int k=-1;k<=1;k++) 
			for(int j=-1;j<=1;j++)
				for(int i=-1;i<=1;i++) {
					// n? => neighbor
					int nx=px+i;
					int ny=py+j;
					int nz=pz+k;

					if(IS_PERIODIC_X) {
						if(nx<0)nx=ClusterNX-1;
						if(nx>=(int)ClusterNX)nx=0;
					} else {
						if((nx<0)||(nx>=(int)ClusterNX)){n++;continue;}
					}
					
					if(IS_PERIODIC_Y) {
						if(ny<0)ny=ClusterNY-1;
						if(ny>=(int)ClusterNY)ny=0;
					} else {
						if((ny<0)||(ny>=(int)ClusterNY)){n++;continue;}
					}

					if(IS_PERIODIC_Z) {
						if(nz<0)nz=ClusterNZ-1;
						if(nz>=(int)ClusterNZ)nz=0;
					} else {
						if((nz<0)||(nz>=(int)ClusterNZ)){n++;continue;}
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

	for(int i=0;i<(int)SystemAtoms.size();i++){
		AtomKeeper Atm(SystemAtoms[i]->GetAtomStorage());
		int nproc=GetGridSize();
		int npp=Atm.GetNAtom()/nproc+NTOLE;

		for(int a=0;a<nproc;a++){ // same atom type for all procs
			ProcAtm[a][i]=new AtomContainer;
			ProcAtm[a][i]->disable_log();
			ProcAtm[a][i]->copy_data(SystemAtoms[i]);
			ProcAtm[a][i]->Allocate(npp,false,AtomKeeper::Referral);
		}

		for(int a=0;a<Atm.GetNAtom();a++) {
			int pn=GetCellIndex(Atm[a].x,Atm[a].y,Atm[a].z);
			ProcAtm[pn][i]->GetAtomStorage().Attach(Atm[a]);
		}
	}

	// one file/proc
	int cntot=0;
	for(int pn=0;pn<GetGridSize();pn++) {
		string fname(BinDirectory+"/p"+as_string(pn)+"-crystal");
		int cnt=0;
		for(int a=0;a<(int)SystemAtoms.size();a++){
			ProcAtm[pn][a]->Save(fname, "a");
			cnt+=ProcAtm[pn][a]->GetNAtom();
		}
		cntot+=cnt;
	}
	
	for(int i=0;i<(int)SystemAtoms.size();i++) {
		strncpy(ProcInfo.AtomNames[i],SystemAtoms[i]->get_name().c_str(),32);
		ProcInfo.AtomNames[i][31]=0x0;
	}
	
	for(int a=0;a<GetGridSize();a++){
		for(int i=0;i<(int)SystemAtoms.size();i++) {
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
	MDSystem::AdjustSystem(); // box is checked here...	
	ProcInfo.Box=Box; // This is the valid system box...
	ProcInfo.TotalAtom=GetNAtom();
	ProcInfo.NumberOfContainers=SystemAtoms.size();
	Root_ArrangeNeighbor();
	Root_DistributeAtoms();
		
	for(int a=0;a<(int)SystemAtoms.size();a++) {
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

	mdassert(ClusterNX*ClusterNY*ClusterNZ==(int)GetGridSize(),
	  "incompatible grid dimension: "+clusdim+
	  " on "+as_string(GetGridSize())+" processor(s)");

	mdassert(GetGridSize()<=MAXPROC,
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
	
	for(int i=0;i<ProcInfo.NumberOfContainers;i++) {
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
		char parstr[4096];
		if(GetRank()==0){
			string ts;
			for(int i=0;i<param.size();i++){
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
	Communicator->Broadcast(&Step, sizeof(int));
	Communicator->Broadcast(&PBoundary, sizeof(int));
	Communicator->Broadcast(&Energy, sizeof(OMD_FLOAT));
	Communicator->Broadcast(&Kinetic, sizeof(OMD_FLOAT));
	Communicator->Broadcast(&Virial, sizeof(OMD_FLOAT));
	Communicator->Broadcast(&Potential, sizeof(OMD_FLOAT));
	Communicator->Broadcast(&BasePotential, sizeof(OMD_FLOAT));
	Communicator->Broadcast(&ElapsedTime, sizeof(OMD_FLOAT));
	Communicator->Broadcast(&Box, sizeof(SysBox));
	Communicator->Broadcast(&write_mode, sizeof(int));
	
	int varsize=RestartVars.size();
	Communicator->Broadcast(&varsize, sizeof(int));

	for(int i=0;i<varsize;i++) {
		struct {
			char l[128];char v[128];
		} rvar;
		
		if(GetRank()==MDROOT) {
			RestartVars[i]->GetLabel().copy(rvar.l,128);
			RestartVars[i]->AsString().copy(rvar.v,128);
		}
		
		Communicator->Broadcast(&rvar, sizeof(rvar));
		
		if(GetRank()!=MDROOT) {
			DataSlot *m=new DataSlot(rvar.l);
			m->SetDefaultData(rvar.v);
			RestartVars.push_back(m);
		}
	}

}

void MDSystemGrid::CreationFunction() {
// Here: conditional: restarting or normal.....
	if(GetRank()==0) {
		if(Mode==CONTINUE_MODE) {
			LoadSimulation(param.string_value("continue"));
		} else {
			CreateSystem();
		}
		Root_Prepare();
	}

	Communicator->SyncProcesses();
	SyncEnvironment();
	SyncVariable();
	Communicator->Broadcast(&ProcInfo,sizeof(StructInfo));

	LoadAtoms(); // all proc loads their own atoms, ThisProcAtomNumber is updated...

	for(int i=0;i<(int)SystemAtoms.size();i++) SystemAtoms[i]->set_id(i);
	CreateGadget();
	SyncEnvironment();
}

void MDSystemGrid::FlattenAtomBox() {
	int na=GetNAtom();
	LocalBuffer.IndexBook.clear();
	GhostBuffer.IndexBook.clear();

	// scan atoms belong to this proc...
	for(int i=0;i<na;i++) {
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
	int nc=SystemAtoms.size();
	int ng=SystemAtomGroups.size();
	int ing[ng];

	for (int i=0;i<nc;i++) SystemAtoms[i]->GetAtomStorage().IndexBook.clear();
	for (int i=0;i<ng;i++) {
		SystemAtomGroups[i]->GetAtomStorage().IndexBook.clear();
		ing[i]=SystemAtomGroups[i]->GetGroupMask();
	}

	for(int i=0;i<LocalAtomNumber;i++) {
		Atom* a=AtomPtr(i);

		if(a->tid>=nc) die("undefined atom found! index="+as_string(i)+
				" tid="+as_string((int)a->tid)+
				" nid="+as_string((int)a->nid)+
				" gid="+as_string((int)a->gid));

		SystemAtoms[(int)(a->tid)]->GetAtomStorage().IndexBook.push(i);

		for(int k=0;k<ng;k++) {
			if(a->gid&ing[k])
				SystemAtomGroups[k]->GetAtomStorage().IndexBook.push(i);
		}

	}

	int cntatom=0;
	for(int i=0;i<nc;i++){
		SystemAtoms[i]->GetAtomStorage().AssignByIndex(AtomStorage);
		cntatom+=SystemAtoms[i]->GetNAtom();
	}

	for(int i=0;i<ng;i++) {
		SystemAtomGroups[i]->GetAtomStorage().AssignByIndex(AtomStorage);
	}

	// checking
	cntatom=Communicator->TakeSUM(cntatom);

	if(cntatom!=(int)ProcInfo.TotalAtom) {
		DumpAtoms(LocalBuffer, "dump-err-local-"+as_string(GetRank()),WM_GHOST|WM_VELOCITY|WM_TID|WM_GID|WM_NID);
		DumpAtoms(GhostBuffer, "dump-err-ghost-"+as_string(GetRank()),WM_GHOST|WM_VELOCITY|WM_TID|WM_GID|WM_NID);
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

void MDSystemGrid::SyncData(int syncmode) {
    if(syncmode&SYNC_SPACE) {
        if(Step%CommRefreshPeriod) {
            Communicator->SendReceive(syncmode);
            DistributeContainers();
        } else {
            Communicator->SendReceive(SYNC_POSITION);
            FlattenAtomBox();
            UpdateRadiusTolerance();
            Communicator->DistributeAtomIndex();
            Communicator->SendReceive(syncmode);
            DistributeContainers();
            Iterator->SetDirty();
        }
    } else {
        Communicator->SendReceive(syncmode);
    }
}

bool MDSystemGrid::CheckRun() {
	InterruptFlag=(int)Communicator->TakeMAX((int)InterruptFlag);
	return MDSystem::CheckRun();
}

void MDSystemGrid::AdjustSystem() { // MDSystem::AdjustSystem is invoked from Root_Prepare
	
	SystemSetting();
	if(PBoundary<0) PBoundary=0;
	
	Box=ProcInfo.Box;
	TotalAtom=ProcInfo.TotalAtom;

	int na=GetNAtom();
	SqrMaxVelocity=0.0;
	for(int i=0;i<na;i++) {
		Atom* a=AtomPtr(i);
		OMD_FLOAT vv=a->vx*a->vx+a->vy*a->vy+a->vz*a->vz;
		if(SqrMaxVelocity<vv) SqrMaxVelocity=vv;
	}

	LocalAtomNumber=na;

}

// FIXME! check radius tollerance vs. verlet...
void MDSystemGrid::UpdateRadiusTolerance() {
	OMD_FLOAT L;
	if(CommRefreshPeriod<1) {L=0.0; CommRefreshPeriod=1;}
	else L=sqrt(SqrMaxVelocity)*(CommRefreshPeriod)*Integrator->TimeStep; 
  L+=(CommAbsoluteTol)?CommRadTol:CommRadTol*Integrator->MaxCutRadius;
  Communicator->SetRadiusTolerance(L);
	Iterator->SetRadiusTolerance(L);
}

void MDSystemGrid::InitGadgets() {
	MDSystem::InitGadgets();
	
	// all the gadget has been initiated..
	UpdateRadiusTolerance();
	Iterator->SetDirty();
	Iterator->SetUpdatePeriod(0); // controlled update by MDSystemGrid...


	// FirstSync... get ghosts from neighbors...
	Communicator->DistributeAtomIndex();
	Communicator->SendReceive(SYNC_SPACE|SYNC_VELOCITY|SYNC_FORCE|SYNC_AUX);
	
	FlattenAtomBox();
	DistributeContainers();
}

void MDSystemGrid::FirstRun() {
	MDSystem::FirstRun();
	Communicator->SyncProcesses();
}

void MDSystemGrid::ErrorHandler(const char* errst) {
	MDSystem::ErrorHandler(errst);
  if(Communicator) {
    if(Communicator->GetRank()==MDROOT) std::cerr<<errst<<std::endl<<std::endl;
	  if(Communicator->GetSize()>1) Communicator->Abort();
    Communicator->Close();
  } else std::cerr<<errst<<std::endl<<std::endl;
}

string MDSystemGrid::GetGridConfiguration() {
	ostringstream ost;

	ost<<"top (26..18):   ";
	for(int i=26;i>=18;i--){
		if((i+1)%3==0)ost<<"\n";
		if(GetNeighborRank(i)>=0) ost << "["<<GetNeighborRank(i) << "]";
		else ost<<"[x]";
		
	}

	ost<<"\nmiddle (17..9): ";
	for(int i=17;i>=9;i--){
		if((i+1)%3==0)ost<<"\n";
		if(GetNeighborRank(i)>=0) ost << "["<<GetNeighborRank(i) << "]";
		else ost<<"[x]";
	}

	ost<< "\nbottom (8..0):  ";
	for(int i=8;i>=0;i--){
		if((i+1)%3==0)ost<<"\n";
		if(GetNeighborRank(i)>=0) ost << "["<<GetNeighborRank(i) << "]";
		else ost<<"[x]";
	}

	return ost.str();	
}


void MDSystemGrid::PrintInfo(ostream& ost) {

	char buf[DEFAULT_TRANSFER_LENGTH];
	
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
		
		for(int i=1;i<GetGridSize();i++) {
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
		Communicator->RawSend(MDROOT,buf,strlen(buf)+1);
	}
	ost << std::scientific << std::setprecision(4);	
}

void MDSystemGrid::PrintInfo(string fname){
	if(GetRank()==MDROOT) {
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
	                         int mode,
	                         bool* AuxPrintable,
	                         char* AuxFormat[],
	                         string AuxNames)
{
	int myrank=GetRank();
	char ready=1;
	
	if(myrank==0) {
		AtomContainer::DumpAtoms(ak, fname, mode, AuxPrintable, AuxFormat, AuxNames);
	} else {
		int cm;
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
	                         int mode,
	                         bool* AuxPrintable,
	                         char* AuxFormat[],
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
	int ncon=ProcInfo.NumberOfContainers;
	for(int i=0;i<ncon;i++) {
		AtomContainer ac;

		int sz=SystemAtoms[i]->GetNAtom();
		int asz[GetGridSize()];
		Communicator->Gather((void*)&sz, (void*)asz, sizeof(int));
		ac.copy_data(SystemAtoms[i]);
		ac.GetAtomStorage().set_name("SAVER BUFFER");
		ac.GetAtomStorage().Copy(SystemAtoms[i]->GetAtomStorage());

		Communicator->SyncProcesses();		
		if(GetRank()==MDROOT) {
			for(int p=1;p<GetGridSize();p++) {
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
				Communicator->RawSend(MDROOT, ac.GetAtomStorage().GetArrayPtr(), sz*sizeof(Atom));
		}
	}
	
	Communicator->SyncProcesses();
	return this;
}

void MDSystemGrid::SaveSimulation(string binfile) {
	if(binfile=="")
		binfile.assign(replace_char(lower_case(get_name()), ' ', '_'));
	if(GetRank()==MDROOT) SaveSimulationConfig(binfile);
	Save(binfile);
}

void MDSystemGrid::PrintMessages(ostream& ost) {
	if(GetRank()==MDROOT) MDSystem::PrintMessages(ost);
}

void MDSystemGrid::ReadParameter() {
	MDSystem::ReadParameter();
	param.peek("comm.refresh", CommRefreshPeriod);
  
  if(param.exist("comm.tolerance")) {
    string str=param.string_value("comm.tolerance");
    if(char_exist(str,'%')) {
      // relative to cut radius
      CommAbsoluteTol=false;
      CommRadTol=as_double(replace_char(str,'%',0x0));
      CommRadTol/=100.0;
    } else {
      // absolute!
      CommAbsoluteTol=true;
      CommRadTol=as_double(str);
    }
  } else CommRadTol=0.1;
     
	if(param.exist("comm.arch")) {
		int nx, ny, nz;
		nx=param.int_value("comm.arch",0);
		ny=param.int_value("comm.arch",1);
		nz=param.int_value("comm.arch",2);
		SetClusterArch(nx,ny,nz);
	}
	param.peek("dir.binary", BinDirectory);
}

