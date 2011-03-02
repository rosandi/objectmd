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
 *
 * File description:
 * Implementation of MDSystem class
 *
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cfloat>
#include <signal.h>
#include <sys/types.h>
#include <unistd.h>

#include <omd/system.hpp>
#include <omd/gadget.hpp>
#include <omd/integrator.hpp>
#include <omd/iterator.hpp>
#include <omd/conditioner.hpp>
#include <omd/detector.hpp>
#include <omd/dataslot.hpp>

//-----------SIGNAL-HANDLING------------------------//
sig_atomic_t sig_int_accept=0;
void SignalHandler(OMD_INT signo){sig_int_accept=(sig_atomic_t)signo;}

// registers

vector<ForceKernel*> force_reg;

// dummy variables: needed in STATIC_MODE
class dummy_integrator: public MDIntegrator {
public:
	void Iterate(){}
	void Integrate(){}
	void Init(MDSystem *WorkSys){System=WorkSys;}
};

class dummy_force: public ForceKernel {
public:
	void CheckCompute(OMD_SIZET at, OMD_SIZET to, OMD_INT atid, OMD_INT toid){}
	void Compute(Atom &at, Atom &to){}
};

//--------------------------------------------------//

void MDSystem::SystemInit(){
	stage=stage_prepare;
	set_name("SIMULATION_SYSTEM");
	register_class(get_name());
	Argc=NULL;Argv=NULL;
	Unit=NULL;
    Potential  = BasePotential =
    Virial     = Kinetic       = 
    Energy     = ElapsedTime   = 
    Box.x0     = Box.x1        =
    Box.y0     = Box.y1        =
    Box.z0     = Box.z1        =
 	SqrMaxVelocity = 0.0;
	BoxImport=false;
	Step   = 0;
    MaxTime    = -1.0;
	SimWallTime = 0.0;
	AtomID=GroupID=ConditionerID=DetectorID=0; // Creation counter
	PBoundary= -1;
    Integrator  = NULL;
    Iterator    = NULL;
	Unificated=false;
	Enumerated=false;
	InterruptFlag=0;
	Mode=NORMAL_MODE;

	FlagBitUsed=0;
	AuxVariableUsed=0;
	ClaimFlagBit(this, "FLAG_ACTIVE");
	ClaimFlagBit(this, "FLAG_OUTSIDE");
	ClaimFlagBit(this, "FLAG_GHOST");
	ClaimFlagBit(this, "FLAG_MIRROR");
	for(OMD_INT i=0;i<MAXAUXVAR;i++)PrintableAux[i]=false;
	AuxFormat=NULL;
	AtomStorage.set_logger(this);
	ExitCode=0;
	Unit=new MDUnit;
	OutputDirectory="";
	silent_mode=false;
}

MDSystem::MDSystem(OMD_INT &argc, OMD_CHAR** &argv){
	SystemInit();
	SetArgument(argc,argv);
}

MDSystem::MDSystem(){
	SystemInit();
}

// The atoms are not freed from each atom container

MDSystem::~MDSystem() {
	for (OMD_SIZET i=0; i<Detectors.size(); i++) {
		delete Detectors[i];
	}
	for (OMD_SIZET i=0; i<Conditioners.size(); i++) {
		delete Conditioners[i];
	}
	for (OMD_SIZET i=0; i<SystemAtoms.size(); i++) {
		delete SystemAtoms[i];
	}
	delete Integrator;

	if(AuxFormat){
		for(OMD_SIZET i=0;i<SAuxFormat.size();i++) MemFree(AuxFormat[i]);
		MemFree(AuxFormat);
	}	
}

/**
 * @brief The main loop
 *
 * This function is the system main loop. This function should be executed after
 * everything is ready.
 *
 * The function also perform dynamic routine of the the system, and repeat the
 * process until the elapsed simulation time reaches the maximum value.
 *
 * The sequence in one loop is:
 *   - Call PreIntegration conditioners
 *   - Do force integration
 *   - Call PostIntegration conditioners
 *   - Measure the system by calling all Detectors
 *   - Increase simulation time
 *
*/

OMD_INT MDSystem::Run(OMD_INT mode) {

	ExitCode=EXIT_SUCCESS;
	Mode=mode;

	try {
		ReadParameters();
		switch(Mode) {
			case NORMAL_MODE: RunNormal(); break;
			case RESTART_MODE: RunRestart(); break;
			case STATIC_MODE: RunStatic(); break;
		}
	}
	catch(const OMD_CHAR* errst) {
		ErrorHandler(errst);
	}
	catch(const string &errst) {
		ErrorHandler(errst.c_str());
	}
	catch(...) {
		ErrorHandler("Undefined error...");
	}
	return ExitCode;
}


void MDSystem::ReadParameters() {
	
	if(param.exist("paramfile")) {
		string paramfilename=param.string_value("paramfile");
		assert(file_exist(paramfilename),
		       "can not find parameter file "+param.string_value("paramfile"));
		param.read(paramfilename);
	} else {
		if(file_exist(DEFAULT_CONFIG_FILENAME)) param.read(DEFAULT_CONFIG_FILENAME);
	}
	
	if(param.exist("--version")) {
		std::cerr << VERSION_INDICATOR << std::endl;
		exit(0);
	}

	if(param.exist("silent")) silent_mode=true;
	// TODO: other system settings??
	param.peek("dir.output", OutputDirectory);
	param.peek("maxtime", MaxTime);

	if(param.exist("restart")) {
		SetRestartFilename(param.string_value("restart"));
		if(file_exist(GetRestartFilename())) Mode=RESTART_MODE;
	}

	if(param.exist("boundary.periodic")) {
		string pst=lower_case(param.string_value("boundary.periodic"));
		if(pst=="no"||pst=="false") {PBoundary=0;}
		else {
			PBoundary=0;
			if(pst.find('x')!=string::npos) PBoundary|=PERIODIC_X;
			if(pst.find('y')!=string::npos) PBoundary|=PERIODIC_Y;
			if(pst.find('z')!=string::npos) PBoundary|=PERIODIC_Z;
		}
	}

	if(param.exist("dump.field")) {
		std::istringstream ist(replace_char(param.string_value("dump.field"), '+', ' '));
		while(ist.good()) {
			string st;
			ist >> st;
			if(lower_case(st)=="id") SetWriteMode(WM_ID);
			if(lower_case(st)=="nid") SetWriteMode(WM_NID);
			if(lower_case(st)=="xid") SetWriteMode(WM_XID);
			if(lower_case(st)=="velocity") SetWriteMode(WM_VELOCITY);
			if(lower_case(st)=="force") SetWriteMode(WM_FORCE);
			if(lower_case(st)=="potential") SetWriteMode(WM_POTENTIAL);
			if(lower_case(st)=="virial") SetWriteMode(WM_VIRIAL);
		}
	}

}

void MDSystem::ErrorHandler(const OMD_CHAR* errst) {
	blog(string(errst), LOGINFO);
	ExitCode=EXIT_FAILURE;
}

void MDSystem::PrintContainerInfo(ostream& ost) {
	ost << "\n*** Atom containers ***\n";
	for (OMD_SIZET i=0; i<SystemAtoms.size(); i++) {
		ost << "=>"; SystemAtoms[i]->PrintInfo(ost);
	}
}

void MDSystem::PrintGadgetInfo(ostream& ost) {

	ost << "\n*** Integrator ***\n"; 
	Integrator->PrintInfo(ost);
	
	ost << "\n*** Conditioners ***\n";
	for (OMD_SIZET i=0; i<Conditioners.size(); i++) {
		ost << "=>"; Conditioners[i]->PrintInfo(ost);
	}
	ost << "\n*** Detectors ***\n";
	for (OMD_SIZET i=0; i<Detectors.size(); i++) {
		ost << "=>"; Detectors[i]->PrintInfo(ost);
	}
	
}

void MDSystem::PrintSystemInfo(ostream& ost) {

	OMD_FLOAT volume=Box.lx*Box.ly*Box.lz;
	OMD_FLOAT tempe=Kinetic/((OMD_FLOAT)TotalAtom);
	OMD_FLOAT press=(Virial/(3.0*volume)+Kinetic/(1.5*volume));
	
	Unit->SetFormat("%0.3f");
	ost
		<< "\n*** Initial status ***"
		<< "\npotential energy= " << Unit->FormatEnergy(BasePotential)
		<< "\nkinetic energy= " << Unit->FormatEnergy(Kinetic)
		<< "\ntemperature= " << Unit->FormatTemperature(tempe)
		<< "\npressure= " << Unit->FormatPressure(press)
		<< "\nforce cut radius= "<< Unit->FormatLength(GetMaxCutRadius())
		<< "\n\n";

	ost << std::fixed << std::setprecision(2) 			
		<< "Number of simulation atoms " << TotalAtom << "\n"
		<< "Simulation boundary box " 
		<< "X(" << Box.x0 << ", " << Box.x1 << ") "
		<< "Y(" << Box.y0 << ", " << Box.y1 << ") "
		<< "Z(" << Box.z0 << ", " << Box.z1 << ")\n"
		<< "Volume =" << Unit->FormatVolume(volume) << "\n";
		
	if(IS_PERIODIC_X||IS_PERIODIC_Y||IS_PERIODIC_Z) {
		ost << "Periodic boundary: ";
		if(IS_PERIODIC_X) ost << "x ";
		if(IS_PERIODIC_Y) ost << "y ";
		if(IS_PERIODIC_Z) ost << "z ";
	} else {ost << "Non-periodic boundary\n\n";}
	Unit->SetFormat("%0.5e");

}

void MDSystem::PrintInfo(ostream& ost)
{
	ost.flush();
	PrintContainerInfo(ost);
	PrintGadgetInfo(ost);
	PrintSystemInfo(ost);
	Unit->SetFormat("%0.5e");
	ost.flush();
}

void MDSystem::PrintInfo(string fname){
	ofstream fout(fname.c_str());
	PrintInfo(fout);
	fout.close();
}

void MDSystem::PrintTime(ostream& ost)
{
	ost << "Simulation started " << ctime(&SimBeginTime);
	ost << "Simulation finished " << ctime(&SimEndTime);
	ost << "Wall time " << SimWallTime<< " seconds\n";
}


void MDSystem::EnumerateAtoms(bool force) {
	if((!Enumerated)||force) {
		OMD_SIZET na=GetNAtom();
		for(OMD_SIZET i=0;i<na;i++) {
			Atoms(i).nid=i;
		}
	}
	Enumerated=true;
}

/**
 This function will unificate all the systems atomic structures to a single
 atom array (Atoms). It will also define the boundary which is the maximum
 and minimum values of coordinates. User must then re-adjust the boundary
 in InitSystem().
 
 After the atom groups is unificated, user must readjust the boundary bos
 in InitSystem() function, since this function has no information about 
 the real box size.
 The box variable only keeps the minimum and maxsimum coordinate values.
 NAtomAll is reserved for ghost atoms. 

*/

void MDSystem::UnificateAtoms() {
    OMD_INT na=0;
    
    for (OMD_SIZET i=0;i<SystemAtoms.size();i++) na+=SystemAtoms[i]->GetNAtom();
    assert(na>0, "no atom to simulate (NAtom=0)");
   
    AtomStorage.Allocate(na);
    AtomStorage.Clear();

    for (OMD_SIZET i=0;i<SystemAtoms.size();i++) {
    	assert(SystemAtoms[i]->M>0.0&&SystemAtoms[i]->Z>0.0,
    	       "uninitialized atom properties (mass and number)",
    	       SystemAtoms[i]->get_name());

		SystemAtoms[i]->SetMaster(this);
		AtomKeeper* ak=&(SystemAtoms[i]->GetAtomStorage());
		OMD_INT n=AtomStorage.GetNAtom();
		AtomStorage.Append(*ak);
		OMD_INT na=ak->GetNAtom();
		ak->Release();
		ak->Allocate(na,AtomKeeper::Referral);
		for(OMD_SIZET i=n;i<AtomStorage.GetNAtom();i++)ak->Attach(AtomStorage[i]);
    }
    TotalAtom=GetNAtom();
	
	for (OMD_SIZET i=0; i<TotalAtom; i++) {
		Atoms(i).potential=0.0;
		Atoms(i).virial=0.0;
	}
	
    Unificated=true;
}

/**
 * This function encapsulate the SystemCreation, and make sure the atom ids are
 * consistent with the containers. The id checking is redundant with AddAtom,
 * but it is required in case of importing.
 * 
 * Task: Create the system atoms and create gadgets.
 * 
 */

void MDSystem::CreationFunction() {

    if(Mode==RESTART_MODE)
    	LoadSimulation();
    else
    	CreateSystem();
    	
    for(OMD_SIZET i=0;i<SystemAtoms.size();i++)SystemAtoms[i]->set_id(i);
    CreateGadget();
}

/** 
 * Adjust the system and give the right geometrical information to
 * user's SystemSetting() function. The box dimension is recalculated
 * regarding the possibly changed box information. Ehen the simulation
 * is reloaded (RESTART_MODE) the system box can not be changed in 
 * SystemSetting().
 * 
 * In static mode, the integrator and force is set to dummy variables
 * if none is specified in the application program.
 * 
 */

void MDSystem::AdjustSystem() {
	
		if((Mode!=RESTART_MODE)&&(!BoxImport))CalcBox();
		SysBox SaveBox=Box;
		SystemSetting();
		if(PBoundary<0) PBoundary=0;
		if(Mode==RESTART_MODE) Box=SaveBox;

		// System's boundary box maybe adjusted in InitSystem!!
		Box.lx=fabs(Box.x1-Box.x0);
		Box.ly=fabs(Box.y1-Box.y0);
		Box.lz=fabs(Box.z1-Box.z0);
		
		Box.hlx=Box.lx/2.0;
		Box.hly=Box.ly/2.0;
		Box.hlz=Box.lz/2.0;
}

/**
 * Message slots from all gadgets are collected to the main MessageSlot. In this
 * way the gadgets take care of their own slots.
 * 
 */

void MDSystem::InitGadgets() {

	// Search for the Iterator. Add if not found.
   	for (OMD_SIZET i=0; i<(Conditioners.size()); i++) {
		if(Conditioners[i]->type_of("iterator"))
			Iterator=dynamic_cast<MDIterator*>(Conditioners[i]);
	}
	
	if(!Iterator){
		warn("using default iterator: MDIterator....");
		AddConditioner(Iterator=new MDIterator());
	}

	if(!Integrator) {
		if(Mode==STATIC_MODE) SetIntegrator(new dummy_integrator);
		else {
			warn("using default integrator: MDIntegrator, with time step 1e-3");
			SetIntegrator(new MDIntegrator(1e-3));
		}
	}

	if(force_reg.empty()) {
		if(Mode==STATIC_MODE) force_reg.push_back(new dummy_force);
		else die("no interaction force....");
	}

	// insert forces and initiate integrator
	for(OMD_SIZET i=0;i<force_reg.size();i++) Integrator->AddForce(force_reg[i]);
	Integrator->Init(this);	

	// Initiate all gadgets
	for (OMD_SIZET i=0;i<Conditioners.size();i++) Conditioners[i]->Init(this);
	for (OMD_SIZET i=0;i<Detectors.size();i++) Detectors[i]->Init(this);

	// Convert aux format string to array of string
	MemAlloc(AuxFormat,SAuxFormat.size()*sizeof(OMD_CHAR*));
	for(OMD_SIZET i=0;i<SAuxFormat.size();i++){
		OMD_SIZET sz=(SAuxFormat[i].size()+1);
		MemAlloc(AuxFormat[i],sz*sizeof(OMD_CHAR));
		memset(AuxFormat[i],'\0',sz);
		SAuxFormat[i].copy(AuxFormat[i],sz);
	}

}

void MDSystem::ArrangeMessageSlots() {
	vector<DataSlot*> tslot;
	OMD_INT pri=1000;
	
	// search minimum priority number
	for(OMD_SIZET i=0; i<MessageSlots.size(); i++) {
		if(pri>MessageSlots[i]->GetPriority())
			pri=MessageSlots[i]->GetPriority();
			tslot.push_back(MessageSlots[i]);
	}
	
	OMD_SIZET adi=0;
	while(adi<MessageSlots.size()) {	
		for(OMD_SIZET i=0; i<tslot.size(); i++) {
			if(tslot[i]->GetPriority()==pri) {
				MessageSlots[adi++]=tslot[i];
				if(adi>MessageSlots.size()) break;
			}
		}
		pri++;
	}
}

/**
 This function do the main initialization of the system. First, it will 
 unificate all the atomic stuctures, which is added previously to system,
 in one array, Atoms. Next the InitSystem() function will be executed and 
 then call all initialization functions, Init(), from the contained classes, 
 in the following  sequence:
<ul>
<li> Forces initialization
<li> Conditioners initialization
<li> Detectors initialization
</ul>
*/

void MDSystem::Initiate() {
	stage=stage_create;
	PreCreation();
    CreationFunction();
   	UnificateAtoms();
   	EnumerateAtoms();
	AdjustSystem();
    PostCreation();
    stage=stage_init;
	InitGadgets();
	ArrangeMessageSlots();
	PushInfo("$ PeriodicBoundary "+as_string(PBoundary));
	
	if(write_mode&WM_ID) {
		// ID? name material_file
		for(OMD_SIZET i=0;i<SystemAtoms.size();i++){
			PushInfo("$ ID"+as_string(i)+" "+
   			   SystemAtoms[i]->get_name()+" "+
			   SystemAtoms[i]->GetMaterialFile());
		}
	}
	time(&SimBeginTime);
}

//-------------SAVING-AND-LOADING--------------

string MDSystem::GetRestartFilename(){	
	if(RestartFileName.empty()) {
		RestartFileName.assign(get_name());
		RestartFileName.append(".bin");
		RestartFileName=replace_char(lower_case(RestartFileName), ' ', '_');
	}
	return RestartFileName;
}

void MDSystem::SetRestartFilename(string filename){
	RestartFileName.assign(filename);
	RestartFileName=replace_char(lower_case(RestartFileName), ' ', '_');
}

// MaxTime excluded... to enable continuing simulation

void MDSystem::SaveVariables(FILE* fl){	
	fwrite(&SimBeginTime, sizeof(time_t), 1, fl);
	fwrite(&Step, sizeof(OMD_INT), 1, fl);
	fwrite(&PBoundary, sizeof(OMD_INT), 1, fl);
	fwrite(&Energy, sizeof(OMD_FLOAT), 1, fl);
	fwrite(&Kinetic, sizeof(OMD_FLOAT), 1, fl);
	fwrite(&Virial, sizeof(OMD_FLOAT), 1, fl);
	fwrite(&Potential, sizeof(OMD_FLOAT), 1, fl);
	fwrite(&BasePotential, sizeof(OMD_FLOAT), 1, fl);
	fwrite(&ElapsedTime, sizeof(OMD_FLOAT), 1, fl);
	fwrite(&Box, sizeof(SysBox), 1, fl);
	fwrite(&write_mode, sizeof(OMD_INT), 1, fl);	
}

void MDSystem::LoadVariables(FILE* fl){
	fread(&SimBeginTime, sizeof(time_t), 1, fl);
	fread(&Step, sizeof(OMD_INT), 1, fl);
	fread(&PBoundary, sizeof(OMD_INT), 1, fl);
	fread(&Energy, sizeof(OMD_FLOAT), 1, fl);
	fread(&Kinetic, sizeof(OMD_FLOAT), 1, fl);
	fread(&Virial, sizeof(OMD_FLOAT), 1, fl);
	fread(&Potential, sizeof(OMD_FLOAT), 1, fl);
	fread(&BasePotential, sizeof(OMD_FLOAT), 1, fl);
	fread(&ElapsedTime, sizeof(OMD_FLOAT), 1, fl);
	fread(&Box, sizeof(SysBox), 1, fl);
	fread(&write_mode, sizeof(OMD_INT), 1, fl);
}

AtomContainer* MDSystem::Save(string fname, string mode) {
	for (OMD_SIZET cr=0; cr<SystemAtoms.size(); cr++)
		SystemAtoms[cr]->Save(fname.c_str(), "a");	
	return this;
}

#define WRITEST(STR, FL) { \
	OMD_INT sz=(STR).size()+1; \
	OMD_CHAR longst[sz]; \
	memset(longst,0,sz); \
	(STR).copy(longst, sz-1); \
	fwrite(&sz, sizeof(OMD_INT), 1, FL); \
	fwrite(longst, sizeof(OMD_CHAR), sz, FL); \
}

#define READST(STR, FL) { \
	OMD_INT sz; \
	fread(&sz, sizeof(OMD_INT), 1, FL); \
	OMD_CHAR longst[sz]; \
	fread(longst, sizeof(OMD_CHAR), sz, FL); \
	(STR).assign(longst); \
}

/**
 * @brief Saving the simulation as a binary file for restarting
 * 
 * This function saves the whole simulation to a binary file. The file
 * is suitable for reloading or restarting the simulation (LoadSimulation()).
 * 
 * At the header of the file a signature and time saved is writen, as two
 * 32 character length zero terminated strings. Then a block of variables
 * follows after the first 64 bytes (see SaveVariables()).
 * The information about the container is writen after the variables, the
 * number of containers, followed by the container names all in 32 characters
 * zero terminated string. 
 * 
 * The content of every container is appended via the call to AtomContainer::Save() 
 * function. This allow reading a single container by AtomContainer::Load(),
 * disregarding the saving order.
 * 
 */

void MDSystem::SaveSimulationConfig(string binfile) {
	OMD_CHAR cname[32];
	FILE* fl = fopen(binfile.c_str(), "w");
	assert(fl, "unable to create binary file to save simulation", binfile);
	
	memset(cname,0,32);
	sprintf(cname,"OMD (c) Y ROSANDI");
	fwrite(cname, sizeof(OMD_CHAR), 32, fl);
	
	// saving time...
	time_t savetime; time(&savetime); memset(cname,0,32);
	sprintf(cname, "%s", ctime(&savetime));
	fwrite(cname, sizeof(OMD_CHAR), 32, fl);
	
	SaveVariables(fl);

	// restart variables
	OMD_SIZET nresvar=RestartVars.size();
	fwrite(&nresvar, sizeof(OMD_INT), 1, fl);		
	for(OMD_SIZET i=0; i<nresvar;i++) {
		WRITEST(RestartVars[i]->GetLabel(), fl);
		WRITEST(RestartVars[i]->AsString(), fl);
	}

	OMD_SIZET NumCont=SystemAtoms.size();
	fwrite(&NumCont, sizeof(OMD_INT), 1, fl);	

	for (OMD_SIZET i=0; i<NumCont; i++) {
		memset(cname,0,32);
		strncpy(cname, SystemAtoms[i]->get_name().c_str(), 31);
		fwrite(cname, sizeof(OMD_CHAR), 32, fl);
	}	
	
	memset(cname,0,32);
	sprintf(cname, "== END OF CONFIG HEADER ===");
	fwrite(cname, sizeof(OMD_CHAR), 32, fl);

	assert(!ferror(fl), "error writing restart file");
	fclose(fl);
}

void MDSystem::SaveSimulation(string binfile) {
	if(binfile=="") 
		binfile.assign(replace_char(lower_case(get_name()), ' ', '_'));
	SaveSimulationConfig(binfile);
	Save(binfile);
}

/**
 * @brief Load a saved simulation binary file
 * 
 * The function is called to create the system on RESTART_MODE. See
 * SaveSimulation() for the information of the order of the stored data.
 * 
 */

void MDSystem::LoadSimulation() {

	#define TOLE 100
	OMD_SIZET NumCont;
	OMD_SIZET nresvar;
	OMD_CHAR cname[32];
	
	string LoadFromFile(GetRestartFilename());

	FILE* fl = fopen(LoadFromFile.c_str(), "r");
	assert(fl, "unable to read restart file", LoadFromFile);
	
	fread(cname, sizeof(OMD_CHAR), 32, fl);
	assert(string("OMD (c) Y ROSANDI")==cname,
	       "wrong binary file"+LoadFromFile);
	       
	fread(cname, sizeof(OMD_CHAR), 32, fl);
	string stm(cname);
	stm=replace_char(stm, '\r', ' ');
	stm=replace_char(stm, '\n', ' ');
	blog("loading binary file '"+LoadFromFile+"' time stamp: "+stm, LOGINFO);

	LoadVariables(fl);
	
	fread(&nresvar, sizeof(OMD_INT), 1, fl);
	for(OMD_SIZET i=0;i<nresvar;i++) {
		string st; string vst;
		READST(st, fl); READST(vst, fl);
		DataSlot *m=new DataSlot(st);
		m->SetDefaultData(vst);
		RestartVars.push_back(m);
	}

	SystemAtoms.clear();
	fread(&NumCont, sizeof(OMD_INT), 1, fl);
	for (OMD_SIZET i=0; i<NumCont; i++) {
		fread(cname, sizeof(OMD_CHAR), 32, fl);
		AddAtom(new AtomContainer)->set_name(cname);
	}

	assert(!ferror(fl), "error loading restart file");
	fclose(fl);

	// Loads all atoms for each atom containers
	for (OMD_SIZET cr=0; cr<NumCont; cr++)
		SystemAtoms[cr]->Load(LoadFromFile);

}

//------------------------------------------------------------

void MDSystem::BorderOffset(OMD_FLOAT dx, OMD_FLOAT dy, OMD_FLOAT dz) {

	Box.x0-=dx;	Box.y0-=dy;	Box.z0-=dz;
	Box.x1+=dx;	Box.y1+=dy;	Box.z1+=dz;

	Box.lx+=2.0*dx;
	Box.ly+=2.0*dy;
	Box.lz+=2.0*dz;

	Box.hlx=Box.lx/2.0;
	Box.hly=Box.ly/2.0;
	Box.hlz=Box.lz/2.0;
}	

void MDSystem::SetBox(OMD_FLOAT x0, OMD_FLOAT y0, OMD_FLOAT z0, OMD_FLOAT x1, OMD_FLOAT y1, OMD_FLOAT z1) {
	Box.x0=x0; Box.y0=y0; Box.z0=z0;
	Box.x1=x1; Box.y1=y1; Box.z1=z1;
	Box.lx=fabs(x1-x0);
	Box.ly=fabs(y1-y0);
	Box.lz=fabs(z1-z0);
	Box.hlx=Box.lx/2.0;
	Box.hly=Box.ly/2.0;
	Box.hlz=Box.lz/2.0;
}

void MDSystem::FirstRun() {
	
	//---------------Run one---------------------
	SyncData(SYNC_ALL); // synchronize data first
	ExecuteConditioners(COND_PRE_INTEGRATION);
	Integrator->Iterate();
	ExecuteConditioners(COND_POST_INTEGRATION);
	//-------------------------------------------
	
	MeasurePotential();
	MeasureKinetic();
	BasePotential=Potential;
	for(OMD_SIZET i=0;i<GetNAtom();i++){Atoms(i).fx=Atoms(i).fy=Atoms(i).fz=0.0;}
}

void MDSystem::ExecuteConditioners(OMD_INT contype) {
	OMD_INT consize=Conditioners.size();
	for(OMD_INT i=0;i<consize;i++){
		Conditioners[i]->Execute(contype);
	}
}

void MDSystem::ExecuteDetectors() {	
	OMD_INT detsize=Detectors.size();
    for(OMD_INT i=0;i<detsize;i++) 
		if(Detectors[i]->IsActive()) Detectors[i]->Detect();
}

void MDSystem::PrintMessages(ostream& ost) {
	OMD_INT printed=0;
	for(OMD_SIZET i=0;i<MessageSlots.size();i++){
		if(MessageSlots[i]->IsPrintable()) {
			ost << MessageSlots[i]->GetFormattedText() << " ";
			printed++;
		}
	}
	if(printed!=0)ost << '\n';
	ost.flush();
}

void MDSystem::CheckBeforeRun() {
	// check output directory
	if(OutputDirectory!="") system(("mkdir -p "+OutputDirectory).c_str());
	
	if(get_type()=="simulation_system") {
		assert(PBoundary==NONPERIODIC, "use MDSystemGrid for PERIODIC boundary simulation!!");
	}

	// Check all gadgets
	assert(Integrator->Check(), "force integrator is not ready");

	for (OMD_SIZET i=0;i<Conditioners.size();i++) 
		assert(Conditioners[i]->Check(), "conditioner not ready", Conditioners[i]->get_name());

	for (OMD_SIZET i=0;i<Detectors.size();i++)
		assert(Detectors[i]->Check(), "detector not ready",Detectors[i]->get_name());

	assert(GetNAtom()>1, "insufficient number of atom", as_string(GetNAtom()));

	for(OMD_SIZET i=0;i<Detectors.size();i++){
		if(OutputDirectory!="") {
			Detectors[i]->SetFilenamePrefix(OutputDirectory+"/");
		}
	}
}

void MDSystem::PrintHeader(ostream& ost) {
	ost << "# _______________________________________________\n#\n"
	    << "# MOLECULAR-DYNAMICS SIMULATION CLASS LIBRARY\n"
	    << "# OBJECT-MD\n#\n"
	    << "# (c) Y. Rosandi, yudi@rosandi.com\n"
	    << "# _______________________________________________\n\n";
}

/**
 * The function takes command line arguments in to the parameter keeper.
 * Some environtment parameters are loaded as well:
 *  - OMD_LIB: the path of OMD library
 *  - OMD_CLASS: the path the OMD class collection
 *  - OMD_TABLE: the path to the tables (potentials, etc)
 */

void MDSystem::SetArgument(OMD_INT &argc, OMD_CHAR** &argv) {
	Argc=&argc;Argv=&argv;
	param.append(argc, argv);
	if(param.exist("load")) {
		SetRestartFilename(param.string_value("load"));
		if(file_exist(GetRestartFilename())) Mode=RESTART_MODE;
	}
}

/**
 * The detectors are called after calculation of system properties. This is to
 * guaranty that these properties are available for them.
 * 
 * The sequence of processes done in this function is the following:
 *    # Execution of scheduller function, Scheduller(). If implemented in the 
 *      descendant, this is a general purpose function reserved for schedulled processes.
 *    # The pre integration conditioners are executed.
 *    # The integration loop is invoked.
 *    # The post integration conditioners are executed.
 *    # The inline function is called. By default, this function is doing nothing.
 *      The descendant may reimplement the function when a special process need to
 *      be inserted inside the main loop.
 *    # Measuring energies, i.e. potential and kinetics.
 *    # The detectors are executed.
 *    # Checking boundary condition.
 */

void MDSystem::RunKernel() {
	stage=stage_run;
	while(CheckRun()) {
		Scheduller();
		ExecuteConditioners(COND_PRE_INTEGRATION);
		Integrator->Integrate();
		InlineFunction();
		MeasurePotential(); MeasureKinetic();
		Energy=Kinetic+Potential-BasePotential;
		ExecuteConditioners(COND_POST_INTEGRATION);
		ExecuteDetectors();
		CheckBoundary();
		if(!silent_mode) PrintMessages(std::cout);
		ElapsedTime+=Integrator->TimeStep;
		Step++;
		CheckInterruption();
	}
	stage=stage_finalize;
}

void MDSystem::RunNormal() {
	Initiate();
	PrintHeader(std::cout);
	CheckBeforeRun();
	BeforeFirstRun(); // user app
	FirstRun();
	BeforeRun();      // user app
	blog("starting simulation "+string(ctime(&SimBeginTime)));
	RunKernel();
	time(&SimEndTime);
	SimWallTime=difftime(SimEndTime,SimBeginTime);
	blog("simulation stops. step="+as_string(Step)+
		 "\nbegin time: "+ctime(&SimBeginTime)+
		 "\nend time: "+ctime(&SimEndTime)+
		 "\nwall time: "+as_string(SimWallTime)
		 );
	if (!InterruptFlag) AfterRun();
}

void MDSystem::RunRestart() {
	Initiate();
	CheckBeforeRun();
	BeforeRun();
	blog("restarting simulation "+string(ctime(&SimBeginTime)));
	RunKernel();
	time(&SimEndTime);
	SimWallTime=difftime(SimEndTime,SimBeginTime);
	blog("simulation stops. step="+as_string(Step)+
		 "\nbegin time: "+ctime(&SimBeginTime)+
		 "\nend time: "+ctime(&SimEndTime)+
		 "\nwall time: "+as_string(SimWallTime)
		 );

	if (!InterruptFlag) AfterRun();
}

/**
 * @brief Runs a molecular static simulation
 * 
 * This function executes only one run of all conditioners and detectors 
 * (no loop). If an integrator and force kernel are assigned in creation
 * time, one force iteration will be executed. The detectors takes place
 * at the end, the stage after the forces and energies are calculated.
 * 
 * If no integrator is supplied, no force iteration will be carried out.
 * 
 */

// dummy integrator and force used if non is specified @ creation stage
void MDSystem::RunStatic() {
	Initiate();
	blog("molecular static simulation\nstarts: "+string(ctime(&SimBeginTime)));
	BeforeRun();
	SyncData(SYNC_ALL);

	for(OMD_SIZET i=0; i<GetNAtom(); i++){
		Atom* a=AtomPtr(i);
		a->fx=a->fy=a->fz=a->virial=a->potential=0.0;		
	}

	ExecuteConditioners(COND_PRE_INTEGRATION);
	Integrator->Iterate();
	ExecuteConditioners(COND_POST_INTEGRATION);
	MeasurePotential();
	MeasureKinetic();
	BasePotential=Potential;
	ExecuteDetectors();
	time(&SimEndTime);
	SimWallTime=difftime(SimEndTime,SimBeginTime);	
	AfterRun();
	blog("ends: "+string(ctime(&SimBeginTime)));
}

/**
  * Compute forces, and measure the potential energy.
  * 
*/

void MDSystem::MeasurePotential()
{
	OMD_INT na=GetNAtom();
	Potential=Virial=0;
	for(OMD_INT i=0;i<na;i++){
		Atom* a=AtomPtr(i);
		if(a->flag&FLAG_GHOST) continue;
		
		if(a->flag&FLAG_ACTIVE) {
			Potential+=a->potential;
			Virial+=a->virial;
		}
	}
	// double sum in force-kernel....
	Potential*=0.5;
	Virial*=0.5;
}

void MDSystem::MeasureKinetic() {
    OMD_FLOAT ValTmp=0.0;
    OMD_FLOAT sqrv;
    Kinetic=0.0;
    SqrMaxVelocity=0.0;
    
    OMD_INT natom=GetNAtom();    
    for(OMD_INT i=0;i<natom;i++) {
    	Atom* a=AtomPtr(i);
		if(a->flag&FLAG_GHOST) continue;
		if(a->flag&FLAG_ACTIVE) {
   			sqrv=a->vx*a->vx+a->vy*a->vy+a->vz*a->vz;   			
       		ValTmp += GetMass(a->id)*sqrv;
       		if(sqrv>SqrMaxVelocity)SqrMaxVelocity=sqrv;

       	}
	}

	Kinetic=0.5*ValTmp;
}

/**
 * \brief Checking atom position on the cell border
 * 
 * This function ensures all atoms positioned inside the calculated cell. Atoms
 * those leave the cell border will be mirrored in the opposite side.
 * Applied to simulations using periodic boundary conditions.
 * 
 * The function is called at the end of run loop, after the positions are updated,
 * and all gadgets executed.
 * 
*/

void MDSystem::CheckBoundary() {
	if (PBoundary) {
		OMD_INT na=GetNAtom();
		for (OMD_INT i=0; i<na; i++) {
			Atom* a=AtomPtr(i);
			if(CheckActive(i)){
				if (PBoundary&PERIODIC_X) {
					if(a->x<Box.x0)a->x+=Box.lx;
					else if(a->x>=Box.x1)a->x-=Box.lx;
				}
				if (PBoundary&PERIODIC_Y) {
					if(a->y<Box.y0)a->y+=Box.ly;
					else if(a->y>=Box.y1)a->y-=Box.ly;
				}
				if (PBoundary&PERIODIC_Z) {
					if (a->z<Box.z0)a->z+=Box.lz;
					else if(a->z>=Box.z1)a->z-=Box.lz;
				}
			}
		}
	}
}

/**
 * \brief Distance correction 
 * 
 * This function correct distances calculation, such as the minimum distance
 * convention. The function is called from CalcSqrDistance() function of MDGadget.
 * 
 */

void MDSystem::BoundaryCorrectDistances(OMD_FLOAT& dx, OMD_FLOAT& dy, OMD_FLOAT& dz){

    	if (PBoundary) {
    		if (PBoundary&PERIODIC_X){
    			if(dx>0.0){if(dx>Box.hlx)dx-=Box.lx;} 
    			else {if(-dx>Box.hlx)dx+=Box.lx;}
    		}
    		
    		if (PBoundary&PERIODIC_Y){
    			if(dy>0.0){if(dy>Box.hly)dy-=Box.ly;} 
    			else {if(-dy>Box.hly)dy+=Box.ly;}
    		}
    		
    		if (PBoundary&PERIODIC_Z){ 
    			if(dz>0.0){if(dz>Box.hlz)dz-=Box.lz;} 
    			else {if(-dz>Box.hlz)dz+=Box.lz;}
    		}
		}

}

MDIntegrator* MDSystem::SetIntegrator(MDIntegrator* itg)
{
	assert(!Integrator,"attempt to reassign integrator");
   	itg->SetSystem(this);
   	itg->set_logger(this);
    itg->set_id(0);
	itg->SetUnit(Unit);
 	Integrator = itg;
    return Integrator;
}

ForceKernel* MDSystem::AddForce(ForceKernel* force) {
	force->set_logger(logger);
	force->SetAtomID(0,0);
	force_reg.push_back(force);
	return force;
}

ForceKernel* MDSystem::AddForce(ForceKernel* force, const OMD_CHAR* from, const OMD_CHAR* to) {
	force->set_logger(logger);
	force->SetAtomID(GetContainerID(from), GetContainerID(to));
	force_reg.push_back(force);
	return force;	
}

Detector* MDSystem::AddDetector(Detector* Detc)
{
	Detc->SetSystem(this);
	Detc->set_logger(this);
	Detc->set_id(DetectorID++); 
	Detc->SetUnit(Unit);
	Detectors.push_back(Detc); 
	return Detc;
}

Conditioner* MDSystem::AddConditioner(Conditioner* Cond)
{
	Cond->SetSystem(this);
	Cond->set_logger(this);
	Cond->set_id(ConditionerID++); 
	Cond->SetUnit(Unit);
	Conditioners.push_back(Cond); 
	return Cond;
}

AtomContainer* MDSystem::AddAtom(AtomContainer* Atm)
{	
	Atm->set_logger(this);
	if(!Atm->created)Atm->Create();
	if(Atm->get_name()=="ATOM_CONTAINER"){
		OMD_CHAR nst[32];
		sprintf(nst,"ATOM_%d",AtomID);
		Atm->set_name(nst);
	}
	
	Atm->set_id(AtomID++);
	SystemAtoms.push_back(Atm); 
	
	blog("added atoms:",
	    "name="+Atm->get_name()+
	    " id="+as_string(Atm->get_id())+
	    " n_atom="+as_string(Atm->GetNAtom()),
	    LOGCREATE);

	return Atm;
}

AtomGroup* MDSystem::AddAtomGroup(string group_name) {
	AtomGroup *ag=new AtomGroup(group_name, this, this);
	ag->set_id(GroupID++);
	SystemAtomGroups.push_back(ag);
	blog("added atom group: "+group_name, LOGCREATE);
	return ag;
}

void MDSystem::ChangeAtomID(OMD_SIZET idx, OMD_INT NewID) 
{
	assert(idx<GetNAtom(), 
	       "changing atom id out of bound index="+as_string(idx)+
	       " to id="+as_string(NewID));
	Atoms(idx).id=NewID;
}

void MDSystem::ChangeAtomID(OMD_SIZET start, OMD_SIZET end, OMD_INT NewID)
{
	assert(start>=0&&end<GetNAtom(),
	       "changing atom id out of bound start="+as_string(start)+
	       " end="+as_string(end)+" to id="+as_string(NewID));
	for (OMD_SIZET i=start;i<=end;i++) Atoms(i).id=NewID;
}

OMD_FLOAT MDSystem::GetMaxCutRadius(){return Integrator->MaxCutRadius;}
OMD_FLOAT MDSystem::GetMaxVelocity(){return sqrt(SqrMaxVelocity);}

OMD_FLOAT MDSystem::GetMass(OMD_INT atomid){
	OMD_FLOAT m;
	try {m=SystemAtoms.at(atomid)->M;}
	catch(...){throw "reading atom mass out of bound";}
	return m;
}

OMD_FLOAT MDSystem::GetNumber(OMD_INT atomid){
	OMD_FLOAT z;
	try {z=SystemAtoms.at(atomid)->Z;}
	catch(...){throw "reading atom number out of bound";}
	return z;
}

OMD_SIZET MDSystem::ClaimFlagBit(MDClass* user,string sinfo) {
	OMD_SIZET a=1<<FlagBitUsed;
	FlagBitUsed++;
	
	string scode(user->get_name());
	if(sinfo!=""){
		scode.append("<");
		scode.append(sinfo);
		scode.append(">");
	}
	
	FlagUser.push_back(scode);

	return a;
}

/**
 * \brief Claims an aux variable in the atom structure
 * 
 * This function is called via the Gadget's implementation of ClaimAuxVariable.
 * The claimed aux variables will be writen to data file and embedded in the 
 * atom structure.
 * 
 * Parameters:
 *   - the claiming user class
 *   - printable, wether the variable printed in the output data file
 *   - tag
 *   - sformat is the printing format of the variable
 * 
 * The function returns the index of the claimed aux variable.
 * 
 */

OMD_SIZET MDSystem::ClaimAuxVariable(
					MDClass* user, 
					bool printable,
					const OMD_CHAR* tag,
					const OMD_CHAR* sformat)
{
		OMD_INT AuxIdx=AuxVariableUsed;
		PrintableAux[AuxIdx]=printable;
		AuxVariableUsed++;		
		AuxUser.push_back(user->get_name());
		
		string auxtag;
		if(tag)auxtag.assign(tag);
		else {
			auxtag.assign("(");
			auxtag.append(user->get_name());
			auxtag.append(")");
		}
		
		AuxNameTag.push_back(auxtag);
		
		string auxformat;
		if(sformat)auxformat.assign(sformat);
		else auxformat.assign("%0.5e");
		SAuxFormat.push_back(auxformat);

		if(AuxVariableUsed>MAXAUXVAR) {
			string serr("Maximum aux variable exceeded. Users:");
			for(OMD_SIZET i=0;i<AuxUser.size();i++){
				serr.append(" ");
				serr.append(AuxUser[i]);
			}
			die(serr);
		}
		return AuxIdx;
}

DataSlot* MDSystem::GetMessageSlot(string slotlabel){
	for(OMD_SIZET i=0;i<MessageSlots.size();i++){
		if(slotlabel.compare(MessageSlots[i]->GetLabel())==0)
			return MessageSlots[i];
	}
	die("Can not find message slot", slotlabel);
	return NULL; // avoids warning
}

OMD_SIZET MDSystem::GetFlagBitMask(const OMD_CHAR* usagecode) {
	OMD_SIZET i;
	bool found=false;
	for(i=0;i<FlagUser.size();i++) {
		string sss("<");sss.append(usagecode);sss.append(">");
		if(FlagUser[i].find(sss)!=string::npos){found=true;break;}
	}
	assert(found, "the 'flag-bit user' mark is not found", usagecode);
	return (1<<i);
}

void MDSystem::AcceptSignal(OMD_INT signo) {
	struct sigaction act;
	memset(&act,0,sizeof(act));
	act.sa_handler=&SignalHandler;
	sigaction(signo,&act,NULL);
}

// FIXME! case and space??
MDGadget* MDSystem::SearchGadget(string name) {
	for(OMD_SIZET i=0;i<Detectors.size();i++) {
		if(Detectors[i]->get_name()==name) return Detectors[i];
	}
	
	for(OMD_SIZET i=0;i<Conditioners.size();i++) {
		if(Conditioners[i]->get_name()==name) return Conditioners[i];
	}
	
	return NULL;
}

OMD_INT MDSystem::GetContainerID(string name){
	OMD_INT cnum=SystemAtoms.size();
	for(OMD_INT i=0;i<cnum;i++) {
		if(SystemAtoms[i]->get_name()==name) return i;
	}
	die("atom "+name+" does not exist");
	return 0;
}

bool MDSystem::OnTime(OMD_FLOAT tm) {
	if((tm>(ElapsedTime-0.5*Integrator->TimeStep))
		&&(tm<(ElapsedTime+0.5*Integrator->TimeStep))) 
		return true;
	return false;
}

bool MDSystem::GadgetExist(MDGadget* gad){
	for(OMD_SIZET i=0;i<Detectors.size();i++){
		if(Detectors[i]==gad) return true;
	}
	for(OMD_SIZET i=0;i<Conditioners.size();i++){
		if(Conditioners[i]==gad) return true;
	}
	return false;
}

// FIXME! Box

/**
 * @brief Importing simulation atoms from a saved text file.
 * 
 * 
 */
AtomContainer* MDSystem::Import(string fname){
	ParamHandler p;
	p.read_pseudo(fname);
	string ids("ID0");

	OMD_INT ia=0;
	while(p.exist(ids)) {
		OMD_INT idx=p.index_of(ids);
		AddAtom(new AtomContainer(p[idx+2]))
			->Import(fname,ia)
			->SetName(p[idx+1]);
		ia++; ids.assign("ID"+as_string(ia));
	}
	
	if(p.exist("Box")) {
		blog("reading saved box geometry", LOGCREATE);
		OMD_INT idx=p.index_of("Box")+1;
		Box.x0=p.double_value(idx++);
		Box.y0=p.double_value(idx++);
		Box.z0=p.double_value(idx++);
		Box.x1=p.double_value(idx++);
		Box.y1=p.double_value(idx++);
		Box.z1=p.double_value(idx);
		BoxImport=true;
	}
	
	if(p.exist("PeriodicBoundary")) {
		PBoundary=p.int_value("PeriodicBoundary");
	}
	
	if(!ia) {
		// if only one type, 'Material' must exist 
		// or assign mat_file before import
		if(mat_file=="") {
			assert(p.exist("Material"), 
			"'Material' definition needed in import file");
			mat_file=p.string_value("Material");
		}
		
		AddAtom(new AtomContainer(mat_file))
			->Import(fname);
	}

	return this;
}

void MDSystem::CheckInterruption(){
	InterruptFlag=(OMD_INT)sig_int_accept;
	if(InterruptFlag) {
		std::cerr << "#### INTERRUPT SIGNAL: " << InterruptFlag << " ####\n";
		blog("interrupt signal caught: "+as_string(InterruptFlag));
		switch(InterruptFlag) {
			case SIGTERM:
				OnInterruptTERM();
				break;
			
			case SIGINT:
				OnInterruptINT();
				break;
			
			case SIGUSR1:
				OnInterruptUSR1();
				break;

			case SIGUSR2: 
				OnInterruptUSR2(); // non-terminating
				InterruptFlag=0;
				break;

			default:
				warn("caught: unhandled signal ("+as_string(InterruptFlag)+")");
		}
	}
}

// default interrupt handlers

void MDSystem::OnInterruptTERM(){
	blog("handling SIGTERM", LOGINFO);
	SaveSimulation();
	ExitCode=OMD_EXIT_SUSPEND;
}

void MDSystem::OnInterruptINT(){
	blog("handling SIGINT", LOGINFO);
	SaveSimulation();
	ExitCode=OMD_EXIT_SUSPEND;
}

void MDSystem::OnInterruptUSR1(){
	blog("handling SIGUSR2", LOGINFO);
	SaveSimulation();
	ExitCode=OMD_EXIT_SUSPEND;
}

// non terminating interrupt USR2
void MDSystem::OnInterruptUSR2(){
	blog("handling SIGUSR1 --- non-terminating interrupt", LOGINFO);
	SaveSimulation();
}
