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

#include <omd/omdtool.h>
#include <omd/system.h>
#include <omd/gadget.h>
#include <omd/integrator.h>
#include <omd/iterator.h>
#include <omd/modify.h>
#include <omd/detector.h>
#include <omd/dataslot.h>

using namespace omd;

//-----------SIGNAL-HANDLING------------------------//
sig_atomic_t sig_int_accept=0;
void SignalHandler(int signo){sig_int_accept=(sig_atomic_t)signo;}

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
	void CheckCompute(int at, int to, int atid, int toid){}
	void Compute(Atom &at, Atom &to){}
};

//--------------------------------------------------//

void MDSystem::SystemInit(){
	set_name("SIMULATION_SYSTEM");
	register_class(get_name());
	Argc=NULL;Argv=NULL;
	Unit=NULL;
	stopped=false;
  Potential  = BasePotential =
  Virial     = Kinetic       = 
  Energy     = ElapsedTime   = 
 	SqrMaxVelocity = 0.0;
	Step   = 0;
  MaxTime    = -1.0;
	SimWallTime = 0.0;
	AtomID=GroupID=ModifyID=DetectorID=0; // Creation counter
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
	for(int i=0;i<MAXAUXVAR;i++)PrintableAux[i]=false;
	AuxFormat=NULL;
	AtomStorage.set_logger(this);
	ExitCode=0;
	Unit=new MDUnit;
	OutputDirectory="";
	silent_mode=false;
	print_every=1;
	ParameterFilename=DEFAULT_CONFIG_FILENAME;
}

MDSystem::MDSystem(int &argc, char** &argv){
	SystemInit();
	SetArgument(argc,argv);
}

MDSystem::MDSystem(){
	SystemInit();
}

// The atoms are not freed from each atom container

MDSystem::~MDSystem() {
	for (int i=0; i<(int)Detectors.size(); i++) {
		delete Detectors[i];
	}
	for (int i=0; i<(int)Modifies.size(); i++) {
		delete Modifies[i];
	}
	for (int i=0; i<(int)SystemAtoms.size(); i++) {
		delete SystemAtoms[i];
	}
	delete Integrator;

	if(AuxFormat){
		for(int i=0;i<(int)SAuxFormat.size();i++) MemFree(AuxFormat[i]);
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
 *   - Call PreIntegration
 *   - Do force integration
 *   - Call PostIntegration
 *   - Measure the system by calling all Detectors
 *   - Increase simulation time
 *
*/

int MDSystem::Run(int mode) {

	ExitCode=EXIT_SUCCESS;
	Mode=mode;
		
	try {
    LoadParameterFile();
		ReadParameter();
		Init(); // user initialization
		switch(Mode) {
			case NORMAL_MODE: RunNormal(); break;
			case CONTINUE_MODE: RunRestart(); break;
			case STATIC_MODE: RunStatic(); break;
		}
	}
	catch(const char* errst) {
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

void MDSystem::LoadParameterFile() {
	if(param.exist("--param")) { // this is from command line
		ParameterFilename=param.string_value("--param");
		mdassert(file_exist(ParameterFilename),
             "can not find parameter file "+param.string_value("--param"));
		blog("reading parameters from ("+ParameterFilename+")");
		param.read(ParameterFilename);
	} else {
		if(file_exist(ParameterFilename)) {
			blog("reading parameters from ("+ParameterFilename+")");			
			param.read(ParameterFilename);
		}
	}
}

void MDSystem::ReadParameter() {
		
	if(param.exist("--version")) {
		std::cerr << VERSION_INDICATOR << std::endl;
		exit(0);
	}

	if(param.exist("--silent")) silent_mode=true;
	
	if(param.exist("signal")) {
		std::istringstream ist(replace_char(param.string_value("signal"),'+',' '));
		while(ist.good()) {
			string st;
			ist >> st;
			if(lower_case(st)=="int")  AcceptSignal(SIGINT);
			if(lower_case(st)=="term") AcceptSignal(SIGTERM);
			if(lower_case(st)=="usr1") AcceptSignal(SIGUSR1);
			if(lower_case(st)=="usr2") AcceptSignal(SIGUSR2);
		}
	}
	
	// TODO: other system settings??
	param.peek("dir.output", OutputDirectory);
	param.peek("time.max", MaxTime);
	
	if(param.exist("continue")) {
		SetBinaryFilename(param.string_value("continue"));
		if(file_exist(GetBinaryFilename())) Mode=CONTINUE_MODE;
		else Mode=NORMAL_MODE;
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
			if(lower_case(st)=="tid") SetWriteMode(WM_TID);
			if(lower_case(st)=="nid") SetWriteMode(WM_NID);
			if(lower_case(st)=="gid") SetWriteMode(WM_GID);
			if(lower_case(st)=="velocity") SetWriteMode(WM_VELOCITY);
			if(lower_case(st)=="force") SetWriteMode(WM_FORCE);
			if(lower_case(st)=="potential") SetWriteMode(WM_POTENTIAL);
			if(lower_case(st)=="virial") SetWriteMode(WM_VIRIAL);
		}
	}
  
  if(param.exist("log")) {
    int loglev=0;
		std::istringstream ist(replace_char(param.string_value("log"), '+', ' '));
		while(ist.good()) {
			string st;
			ist >> st;
			if(lower_case(st)=="memory") loglev|=LOGMEMORY;
			if(lower_case(st)=="create") loglev|=LOGCREATE;
			if(lower_case(st)=="destroy") loglev|=LOGDESTROY;
			if(lower_case(st)=="info") loglev|=LOGINFO;
			if(lower_case(st)=="warning") loglev|=LOGWARNING;
		}
    log_flagset(loglev);
  }
  
  if(param.exist("mode.static")) Mode=STATIC_MODE;

}

// FIXME! Untested....

/**
 This function should be reimplemented in the descendant class. If not, 
 the default system creation is defined by the program parameters.
 - import (crystal file) : import from a crystal file.
 - load (binary file) : load a atoms from a binary file. This is not continue mode, since
   the simulation times is reset to zero.
 */

void MDSystem::CreateSystem() {

	if(param.exist("import")) { // import crystal

		string fname=param.string_value("import");
		mdassert(file_exist(fname), "can not find file to import ("+fname+")");
		blog("importing file "+fname);
		if(param.peek("material", mat_file)) blog("setting material to "+mat_file);
		Import(fname);

	} else if(param.exist("load")) { // load from binary file

		string fname=param.string_value("load");
		mdassert(file_exist(fname), "can not find binary file to load ("+fname+")");
		blog("loading file "+fname);
		LoadSimulation(fname);
		EnumerateAtoms();
		ResetSimulationTime();

	}

}

void MDSystem::ErrorHandler(const char* errst) {
	blog(string(errst), LOGINFO);
	ExitCode=EXIT_FAILURE;
}

void MDSystem::PrintContainerInfo(ostream& ost) {
	ost << "\n*** Atom containers ***\n";
	for (int i=0; i<(int)SystemAtoms.size(); i++) {
		ost << "=>"; SystemAtoms[i]->PrintInfo(ost);
	}
}

void MDSystem::PrintGadgetInfo(ostream& ost) {

	ost << "\n*** Integrator ***\n"; 
	Integrator->PrintInfo(ost);
	
	ost << "\n*** Modifies ***\n";
	for (int i=0; i<(int)Modifies.size(); i++) {
		ost << "=>"; Modifies[i]->PrintInfo(ost);
	}
	ost << "\n*** Detectors ***\n";
	for (int i=0; i<(int)Detectors.size(); i++) {
		ost << "=>"; Detectors[i]->PrintInfo(ost);
	}
	
}

void MDSystem::PrintSystemInfo(ostream& ost) {

	double volume=Box.lx*Box.ly*Box.lz;
	double tempe=Kinetic/((double)TotalAtom);
	double press=(Virial+2.0*Kinetic)/(3.0*volume);
	
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
		int na=GetNAtom();
		for(int i=0;i<na;i++) {
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
  int na=0;
  
  for (int i=0;i<(int)SystemAtoms.size();i++) na+=SystemAtoms[i]->GetNAtom();
  mdassert(na>0, "no atom to simulate (NAtom=0)");
  
  AtomStorage.Allocate(na);
  AtomStorage.Clear();
  
  for (int i=0;i<(int)SystemAtoms.size();i++) {
    mdassert(SystemAtoms[i]->M>0.0&&SystemAtoms[i]->Z>0.0,
    	       "uninitialized atom properties (mass and number)",
    	       SystemAtoms[i]->get_name());
    
		SystemAtoms[i]->SetMaster(this);
		AtomKeeper* ak=&(SystemAtoms[i]->GetAtomStorage());
		int n=AtomStorage.GetNAtom();
		AtomStorage.Append(*ak);
		int na=ak->GetNAtom();
		ak->Release();
		ak->Allocate(na,AtomKeeper::Referral);
		for(int i=n;i<AtomStorage.GetNAtom();i++)ak->Attach(AtomStorage[i]);
  }
  TotalAtom=GetNAtom();
	
	for (int i=0; i<TotalAtom; i++) {
		Atoms(i).potential=0.0;
		Atoms(i).virial=0.0;
	}
	created=true;
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

    if(Mode==CONTINUE_MODE) {
    	LoadSimulation(param.string_value("continue"));
	} else {
    	CreateSystem();
	}

    for(int i=0;i<(int)SystemAtoms.size();i++) SystemAtoms[i]->set_id(i);
    CreateGadget();
}

/**
 Calculate box dimension by taking maximum and minimum extents of the
 SystemAtoms. The CalcBox() function of all containers is executed.
 */

SysBox& MDSystem::CalcBox() {
	SysBox ba,bb;

	bb.x0=bb.y0=bb.z0=DBL_MAX;
	bb.x1=bb.y1=bb.z1=-DBL_MAX;

	// find bounding box...
	for (int i=0; i<(int)SystemAtoms.size(); i++) {
		ba=SystemAtoms[i]->GetBox();
		if(bb.x0>ba.x0)bb.x0=ba.x0;
		if(bb.y0>ba.y0)bb.y0=ba.y0;
		if(bb.z0>ba.z0)bb.z0=ba.z0;
		if(bb.x1<ba.x1)bb.x1=ba.x1;
		if(bb.y1<ba.y1)bb.y1=ba.y1;
		if(bb.z1<ba.z1)bb.z1=ba.z1;
	}
	
	Box=bb;
	Box.lx=fabs(Box.x1-Box.x0);
	Box.ly=fabs(Box.y1-Box.y0);
	Box.lz=fabs(Box.z1-Box.z0);	
	Box.hlx=Box.lx/2.0;
	Box.hly=Box.ly/2.0;
	Box.hlz=Box.lz/2.0;
	
	return Box;
}

/** 
 * Adjust the system and give the right geometrical information to
 * user's SystemSetting() function. In CONTINUE_MODE
 * this function is not called.
 * 
 */

void MDSystem::AdjustSystem() {
	if(Box.undefined()) CalcBox();
	if(Mode!=CONTINUE_MODE) SystemSetting();
	if(PBoundary<0) PBoundary=0;
}

/**
 * Message slots from all gadgets are collected to the main MessageSlot. In this
 * way the gadgets take care of their own slots.
 * In static mode, the integrator and force are set to a dummy class
 * if they are not assign in CreateGadget() function.
 */

void MDSystem::InitGadgets() {
	
	// Synchronize atom's group flag...
	
	for(int i=0;i<(int)SystemAtomGroups.size();i++) {
		if(!(SystemAtomGroups[i]->created)) SystemAtomGroups[i]->Commit();
		SystemAtomGroups[i]->SyncAtomGroupMask();
    SystemAtomGroups[i]->SetMaster(this);
	}
	
	MeasureKinetic(); // needed in some initializations
	
	// Search for the Iterator. Add if not found.
   	for (int i=0; i<(int)Modifies.size(); i++) {
		if(Modifies[i]->type_of("iterator"))
			Iterator=dynamic_cast<MDIterator*>(Modifies[i]);
	}
	
	if(!Iterator){
		warn("using default iterator: MDIterator, half-loop....");
		AddModify(Iterator=new MDIterator());
		Iterator->EnableHalfLoop();
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
	for(int i=0;i<(int)force_reg.size();i++) Integrator->AddForce(force_reg[i]);
	Integrator->Init(this);	
  
  if(Mode!=STATIC_MODE)	
    mdassert(Integrator->MaxCutRadius>0.0, "bad cut radius: "+as_string(Integrator->MaxCutRadius));

	// Initiate all gadgets
	for (int i=0;i<(int)Modifies.size();i++) Modifies[i]->Init(this);
	for (int i=0;i<(int)Detectors.size();i++) Detectors[i]->Init(this);

	// Convert aux format string to array of string
	MemAlloc(AuxFormat,SAuxFormat.size()*sizeof(char*));
	for(int i=0;i<(int)SAuxFormat.size();i++){
		int sz=(SAuxFormat[i].size()+1);
		MemAlloc(AuxFormat[i],sz*sizeof(char));
		memset(AuxFormat[i],'\0',sz);
		SAuxFormat[i].copy(AuxFormat[i],sz);
	}

}

void MDSystem::ArrangeMessageSlots() {
	vector<DataSlot*> tslot;
	int pri=1000;
	
	// search minimum priority number
	for(int i=0; i<(int)MessageSlots.size(); i++) {
		if(pri>MessageSlots[i]->GetPriority())
			pri=MessageSlots[i]->GetPriority();
			tslot.push_back(MessageSlots[i]);
	}
	
	int adi=0;
	while(adi<(int)MessageSlots.size()) {
		for(int i=0; i<(int)tslot.size(); i++) {
			if(tslot[i]->GetPriority()==pri) {
				MessageSlots[adi++]=tslot[i];
				if(adi>(int)MessageSlots.size()) break;
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
<li> Force initialization
<li> Modify initialization
<li> Detector initialization
</ul>
*/

void MDSystem::Initiate() {
  PreCreation();
  CreationFunction();
  UnificateAtoms();
  EnumerateAtoms();
  AdjustSystem();
  CreateGroup();
  PostCreation();
  InitGadgets();
	ArrangeMessageSlots();
	PushInfo("$ PeriodicBoundary "+as_string(PBoundary));
	
	if(write_mode&WM_TID) {
		// ID? name material_file
		for(int i=0;i<(int)SystemAtoms.size();i++){
			PushInfo("$ ID"+as_string(i)+" "+
   			   SystemAtoms[i]->get_name()+" "+
			   SystemAtoms[i]->GetMaterialFile());
		}
	}
	time(&SimBeginTime);
}

//-------------SAVING-AND-LOADING--------------

string MDSystem::GetBinaryFilename(){	
	if(BinaryFilename.empty()) {
		BinaryFilename.assign(get_name());
		BinaryFilename.append(".bin");
		BinaryFilename=replace_char(lower_case(BinaryFilename), ' ', '_');
	}
	return BinaryFilename;
}

void MDSystem::SetBinaryFilename(string filename){
	BinaryFilename.assign(filename);
	BinaryFilename=replace_char(lower_case(BinaryFilename), ' ', '_');
}

// MaxTime excluded... to enable continuing simulation

void MDSystem::SaveVariables(FILE* fl){	
	fwrite(&SimBeginTime, sizeof(time_t), 1, fl);
	fwrite(&Step, sizeof(int), 1, fl);
	fwrite(&PBoundary, sizeof(int), 1, fl);
	fwrite(&Energy, sizeof(double), 1, fl);
	fwrite(&Kinetic, sizeof(double), 1, fl);
	fwrite(&Virial, sizeof(double), 1, fl);
	fwrite(&Potential, sizeof(double), 1, fl);
	fwrite(&BasePotential, sizeof(double), 1, fl);
	fwrite(&ElapsedTime, sizeof(double), 1, fl);
	fwrite(&Box, sizeof(SysBox), 1, fl);
	fwrite(&write_mode, sizeof(int), 1, fl);	
}

void MDSystem::LoadVariables(FILE* fl){
	fread(&SimBeginTime, sizeof(time_t), 1, fl);
	fread(&Step, sizeof(int), 1, fl);
	fread(&PBoundary, sizeof(int), 1, fl);
	fread(&Energy, sizeof(double), 1, fl);
	fread(&Kinetic, sizeof(double), 1, fl);
	fread(&Virial, sizeof(double), 1, fl);
	fread(&Potential, sizeof(double), 1, fl);
	fread(&BasePotential, sizeof(double), 1, fl);
	fread(&ElapsedTime, sizeof(double), 1, fl);
	fread(&Box, sizeof(SysBox), 1, fl);
	fread(&write_mode, sizeof(int), 1, fl);
}

AtomContainer* MDSystem::Save(string fname, string mode) {
	for (int cr=0; cr<(int)SystemAtoms.size(); cr++)
		SystemAtoms[cr]->Save(fname.c_str(), "a");	
	return this;
}

#define WRITEST(STR, FL) { \
	int sz=(STR).size()+1; \
	char longst[sz]; \
	memset(longst,0,sz); \
	(STR).copy(longst, sz-1); \
	fwrite(&sz, sizeof(int), 1, FL); \
	fwrite(longst, sizeof(char), sz, FL); \
}

#define READST(STR, FL) { \
	int sz; \
	fread(&sz, sizeof(int), 1, FL); \
	char longst[sz]; \
	fread(longst, sizeof(char), sz, FL); \
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
	char cname[32];
	FILE* fl = fopen(binfile.c_str(), "w");
	mdassert(fl, "unable to create binary file to save simulation", binfile);
	
	memset(cname,0,32);
	sprintf(cname,"OMD (c) Y ROSANDI");
	fwrite(cname, sizeof(char), 32, fl);
	
	// saving time...
	time_t savetime; time(&savetime); memset(cname,0,32);
	sprintf(cname, "%s", ctime(&savetime));
	fwrite(cname, sizeof(char), 32, fl);
	
	SaveVariables(fl);

	// continue variables
	int nresvar=RestartVars.size();
	fwrite(&nresvar, sizeof(int), 1, fl);		
	for(int i=0; i<nresvar;i++) {
		WRITEST(RestartVars[i]->GetLabel(), fl);
		WRITEST(RestartVars[i]->AsString(), fl);
	}

	int NumCont=SystemAtoms.size();
	fwrite(&NumCont, sizeof(int), 1, fl);	

	for (int i=0; i<NumCont; i++) {
		memset(cname,0,32);
		strncpy(cname, SystemAtoms[i]->get_name().c_str(), 31);
		fwrite(cname, sizeof(char), 32, fl);
	}	
	
	memset(cname,0,32);
	sprintf(cname, "== END OF CONFIG HEADER ===");
	fwrite(cname, sizeof(char), 32, fl);

	mdassert(!ferror(fl), "error writing binary file");
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
 * The function is called to create the system on CONTINUE_MODE. See
 * SaveSimulation() for the information of the order of the stored data.
 * 
 */

void MDSystem::LoadSimulation(string LoadFromFile) {

	#define TOLE 100
	int NumCont;
	int nresvar;
	char cname[32];
	mdassert(file_exist(LoadFromFile), "(LOAD) binary file doesn't exist: "+LoadFromFile);
	FILE* fl = fopen(LoadFromFile.c_str(), "r");
	mdassert(fl, "unable to read binary file", LoadFromFile);
	
	fread(cname, sizeof(char), 32, fl);
	mdassert(string("OMD (c) Y ROSANDI")==cname,
	       "wrong binary file"+LoadFromFile);
	       
	fread(cname, sizeof(char), 32, fl);
	string stm(cname);
	stm=replace_char(stm, '\r', ' ');
	stm=replace_char(stm, '\n', ' ');
	blog("loading binary file '"+LoadFromFile+"' time stamp: "+stm, LOGINFO);

	LoadVariables(fl);
	
	fread(&nresvar, sizeof(int), 1, fl);
	for(int i=0;i<nresvar;i++) {
		string st; string vst;
		READST(st, fl); READST(vst, fl);
		DataSlot *m=new DataSlot(st);
		m->SetDefaultData(vst);
		RestartVars.push_back(m);
	}

	SystemAtoms.clear();
	fread(&NumCont, sizeof(int), 1, fl);
	for (int i=0; i<NumCont; i++) {
		fread(cname, sizeof(char), 32, fl);
		AddAtom(new AtomContainer)->set_name(cname);
	}

	mdassert(!ferror(fl), "error loading binary file");
	fclose(fl);

	// Loads all atoms for each atom containers
	for (int cr=0; cr<NumCont; cr++) {
		SystemAtoms[cr]->Load(LoadFromFile);
	}
}

/**
 This function recalculate the box boundary of the system and add offsets
 to the sides of the box. The box is defined by the minimum and maximum values
 of the atom coordinates.
 */

void MDSystem::BorderOffset(double dx, double dy, double dz) {
	AtomContainer::CalcBox();
	Box.x0-=dx;	Box.y0-=dy;	Box.z0-=dz;
	Box.x1+=dx;	Box.y1+=dy;	Box.z1+=dz;

	Box.lx+=2.0*dx;
	Box.ly+=2.0*dy;
	Box.lz+=2.0*dz;

	Box.hlx=Box.lx/2.0;
	Box.hly=Box.ly/2.0;
	Box.hlz=Box.lz/2.0;
}	

void MDSystem::SetBox(double x0, double y0, double z0, double x1, double y1, double z1) {
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
	ExecuteModifies(MODIFY_PRE_INTEGRATION);
	Integrator->Iterate();
	ExecuteModifies(MODIFY_POST_INTEGRATION);
	//-------------------------------------------
	
	MeasurePotential();
	MeasureKinetic();
	BasePotential=Potential;
	for(int i=0;i<GetNAtom();i++){Atoms(i).fx=Atoms(i).fy=Atoms(i).fz=0.0;}
}

void MDSystem::ExecuteModifies(int contype) {
	int consize=Modifies.size();
	for(int i=0;i<consize;i++){
		Modifies[i]->Execute(contype);
	}
}

void MDSystem::ExecuteDetectors() {	
	int detsize=Detectors.size();
    for(int i=0;i<detsize;i++) 
		if(Detectors[i]->IsActive()) Detectors[i]->Detect();
}

void MDSystem::PrintMessages(ostream& ost) {
  if(!Step) return;
	int printed=0;
	for(int i=0;i<(int)MessageSlots.size();i++){
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
		mdassert(PBoundary==NONPERIODIC, "use MDSystemGrid for PERIODIC boundary simulation!!");
	}

	// Check all gadgets
	mdassert(Integrator->Check(), "force integrator is not ready");

	for (int i=0;i<(int)Modifies.size();i++)
		mdassert(Modifies[i]->Check(), "modify not ready", Modifies[i]->get_name());

	for (int i=0;i<(int)Detectors.size();i++)
		mdassert(Detectors[i]->Check(), "detector not ready",Detectors[i]->get_name());

	for(int i=0;i<(int)Detectors.size();i++){
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
 * Some environment parameters are loaded as well by MDClass:
 *  - OMD_LIB: the path of OMD library
 *  - OMD_CLASS: the path the OMD class collection
 *  - OMD_TABLE: the path to the tables (potentials, etc)
 */

void MDSystem::SetArgument(int &argc, char** &argv) {
	Argc=&argc;Argv=&argv;
	param.append(argc, argv);
}

/**
 * The detectors are called after calculation of system properties. This is to
 * guaranty that these properties are available for them.
 * 
 * The sequence of processes done in this function is the following:
 *    # Execution of scheduller function, Scheduler(). If implemented in the 
 *      descendant, this is a general purpose function reserved for schedulled processes.
 *    # The pre integration modifies are executed.
 *    # The integration loop is invoked.
 *    # The post integration modifies are executed.
 *    # The inline function is called. By default, this function is doing nothing.
 *      The descendant may reimplement the function when a special process need to
 *      be inserted inside the main loop.
 *    # Measuring energies, i.e. potential and kinetics.
 *    # The detectors are executed.
 *    # Checking boundary condition.
 */

void MDSystem::RunKernel() {
	while(CheckRun()) {
		Scheduler();
		ExecuteModifies(MODIFY_PRE_INTEGRATION);
		Integrator->Integrate();
		InlineFunction();
		MeasurePotential(); MeasureKinetic();
		Energy=Kinetic+Potential-BasePotential;
		ExecuteModifies(MODIFY_POST_INTEGRATION);
		ExecuteDetectors();
		CheckBoundary();
		
		if(!silent_mode && !(Step%print_every))
			PrintMessages(std::cout);
		
		ElapsedTime+=Integrator->TimeStep;
		Step++;
		CheckInterruption();
	}
}

void MDSystem::RunNormal() {
	Initiate();
	PrintHeader(std::cout);
	CheckBeforeRun();
	BeforeFirstRun(); // user app
	FirstRun();
	BeforeRun();      // user app
	blog("starting simulation "+replace_char(ctime(&SimBeginTime),'\n',' '));
	RunKernel();
	time(&SimEndTime);
	SimWallTime=difftime(SimEndTime,SimBeginTime);
	blog("simulation stops. step="+as_string(Step)+
		 "\nbegin time: "+replace_char(ctime(&SimBeginTime),'\n',' ')+
		 "\nend time: "+replace_char(ctime(&SimEndTime),'\n',' ')+
		 "\nwall time: "+as_string(SimWallTime)
		 );
	if (!InterruptFlag) AfterRun();
}

void MDSystem::RunRestart() {
	Initiate();
	CheckBeforeRun();
	BeforeRun();
	blog("restarting simulation "+replace_char(ctime(&SimBeginTime), '\n',' '));
	RunKernel();
	time(&SimEndTime);
	SimWallTime=difftime(SimEndTime,SimBeginTime);
	blog("simulation stops. step="+as_string(Step)+
		 "\nbegin time: "+replace_char(ctime(&SimBeginTime),'\n', ' ')+
		 "\nend time: "+replace_char(ctime(&SimEndTime),'\n',' ')+
		 "\nwall time: "+as_string(SimWallTime)
		 );

	if (!InterruptFlag) AfterRun();
}

/**
 * @brief Runs a molecular static simulation
 * 
 * This function executes only one run of all modifies and detectors 
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
	blog("molecular static simulation\nstarted: "+
		 replace_char(ctime(&SimBeginTime),'\n',' '));
	BeforeRun();
	SyncData(SYNC_ALL);

	for(int i=0; i<GetNAtom(); i++){
		Atom* a=AtomPtr(i);
		a->fx=a->fy=a->fz=a->virial=a->potential=0.0;		
	}

	ExecuteModifies(MODIFY_PRE_INTEGRATION);
	Integrator->Iterate();
	ExecuteModifies(MODIFY_POST_INTEGRATION);
	MeasurePotential();
	MeasureKinetic();
	BasePotential=Potential;
	ExecuteDetectors();
	time(&SimEndTime);
	SimWallTime=difftime(SimEndTime,SimBeginTime);	
	AfterRun();
	blog("finished: "+replace_char(ctime(&SimEndTime),'\n',' '));
	blog("walltime: "+as_string(SimWallTime)+" seconds");
}

void MDSystem::MeasurePotential() {
	int na=GetNAtom();
	Potential=Virial=0;
	for(int i=0;i<na;i++){
		Atom* a=AtomPtr(i);
		if(a->flag&FLAG_GHOST) continue;
		
		if(a->flag&FLAG_ACTIVE) {
			Potential+=a->potential;
			Virial+=a->virial;
		}
	}
	Virial/=2.0;
}

void MDSystem::MeasureKinetic() {
    double ValTmp=0.0;
    double sqrv;
    Kinetic=0.0;
    SqrMaxVelocity=0.0;
    
    int natom=GetNAtom();
    for(int i=0;i<natom;i++) {
    	Atom* a=AtomPtr(i);
		if(a->flag&FLAG_GHOST) continue;
		if(a->flag&FLAG_ACTIVE) {
   			sqrv=a->vx*a->vx+a->vy*a->vy+a->vz*a->vz;   			
       		ValTmp += GetMass(a)*sqrv;
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
		int na=GetNAtom();
		for (int i=0; i<na; i++) {
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

void MDSystem::BoundaryCorrectDistances(double& dx, double& dy, double& dz){

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
	mdassert(!Integrator,"attempt to reassign integrator");
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

ForceKernel* MDSystem::AddForce(ForceKernel* force, string from, string to) {
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

Modify* MDSystem::AddModify(Modify* Cond) {
	Cond->SetSystem(this);
	Cond->set_logger(this);
	Cond->set_id(ModifyID++); 
	Cond->SetUnit(Unit);
	Modifies.push_back(Cond); 
	return Cond;
}

AtomContainer* MDSystem::AddAtom(AtomContainer* Atm) {	
	Atm->set_logger(this);
	
	if(!Atm->created) Atm->Create();
	
	if(Atm->get_name()=="ATOM_CONTAINER"){
		char nst[32];
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
	ag->SetGroupMask(SystemAtomGroups.size());
	SystemAtomGroups.push_back(ag);
	blog("added atom group: "+group_name, LOGCREATE);
	return ag;
}

void MDSystem::ChangeAtomID(int idx, int NewID) 
{
	mdassert(idx<GetNAtom(), 
	       "changing atom id out of bound index="+as_string(idx)+
	       " to id="+as_string(NewID));
	Atoms(idx).tid=NewID;
}

void MDSystem::ChangeAtomID(int start, int end, int NewID)
{
	mdassert(start>=0&&end<GetNAtom(),
	       "changing atom id out of bound start="+as_string(start)+
	       " end="+as_string(end)+" to id="+as_string(NewID));
	for (int i=start;i<=end;i++) Atoms(i).tid=NewID;
}

double MDSystem::GetMaxCutRadius(){return Integrator->MaxCutRadius;}
double MDSystem::GetMaxVelocity(){return sqrt(SqrMaxVelocity);}

double MDSystem::GetMass(int idx){return SystemAtoms[Atoms(idx).tid]->M;}
double MDSystem::GetMass(Atom* a){return SystemAtoms[a->tid]->M;}
double MDSystem::GetMass(Atom& a){return SystemAtoms[a.tid]->M;}

double MDSystem::GetZ(int idx){return SystemAtoms[Atoms(idx).tid]->Z;}

int MDSystem::ClaimFlagBit(MDClass* user,string sinfo) {
	int a=1<<FlagBitUsed;
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

int MDSystem::ClaimAuxVariable(
					MDClass* user, 
					bool printable,
					const char* tag,
					const char* sformat)
{
		int AuxIdx=AuxVariableUsed;
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
			for(int i=0;i<(int)AuxUser.size();i++){
				serr.append(" ");
				serr.append(AuxUser[i]);
			}
			die(serr);
		}
		return AuxIdx;
}

DataSlot* MDSystem::GetMessageSlot(string slotlabel){
	for(int i=0;i<(int)MessageSlots.size();i++){
		if(slotlabel.compare(MessageSlots[i]->GetLabel())==0)
			return MessageSlots[i];
	}
	die("Can not find message slot", slotlabel);
	return NULL; // avoids warning
}

int MDSystem::GetFlagBitMask(const char* usagecode) {
	int i;
	bool found=false;
	for(i=0;i<(int)FlagUser.size();i++) {
		string sss("<");sss.append(usagecode);sss.append(">");
		if(FlagUser[i].find(sss)!=string::npos){found=true;break;}
	}
	mdassert(found, "the 'flag-bit user' mark is not found", usagecode);
	return (1<<i);
}

MDGadget* MDSystem::SearchGadget(string name) {
	for(int i=0;i<(int)Modifies.size();i++) {
		if(Modifies[i]->get_name()==name) return Modifies[i];
	}

	for(int i=0;i<(int)Detectors.size();i++) {
		if(Detectors[i]->get_name()==name) return Detectors[i];
	}
		
	return NULL;
}

int MDSystem::GetContainerID(string name){
	int cnum=SystemAtoms.size();
	for(int i=0;i<cnum;i++) {
		if(SystemAtoms[i]->get_name()==name) return i;
	}
	die("atom "+name+" does not exist");
	return 0;
}

bool MDSystem::OnTime(double tm) {
	if((tm>(ElapsedTime-0.5*Integrator->TimeStep))
		&&(tm<(ElapsedTime+0.5*Integrator->TimeStep))) 
		return true;
	return false;
}

bool MDSystem::OnStep(int step) {
	if(Step==step) return true;
	return false;
}

bool MDSystem::GadgetExist(MDGadget* gad){
	for(int i=0;i<(int)Detectors.size();i++){
		if(Detectors[i]==gad) return true;
	}
	for(int i=0;i<(int)Modifies.size();i++){
		if(Modifies[i]==gad) return true;
	}
	return false;
}

/**
 * @brief Importing simulation atoms from a saved text file.
 * 
 * The import file must contain material definitions (pseudo command).
 * - multi-typed crystal: IDn name material
 * - single-typed: Material material
 */

AtomContainer* MDSystem::Import(string fname){
	ParamHandler p;
	p.read_pseudo(fname);
	string ids("ID0");
	int ia=0;
	while(p.exist(ids)) { // ID* name material, can be overriden by parameter
		int idx=p.index_of(ids);
		AddAtom(new AtomContainer(p[idx+2]))
			->Import(fname,ia)
			->SetName(p[idx+1]);
		ia++;
		ids.assign("ID"+as_string(ia));
	}
	
	if(p.exist("Box")) {
		blog("reading saved box geometry", LOGCREATE);
		int idx=p.index_of("Box")+1;
		Box.x0=p.double_value(idx++);
		Box.y0=p.double_value(idx++);
		Box.z0=p.double_value(idx++);
		Box.x1=p.double_value(idx++);
		Box.y1=p.double_value(idx++);
		Box.z1=p.double_value(idx);
		Box.lx=fabs(Box.x1-Box.x0);
		Box.ly=fabs(Box.y1-Box.y0);
		Box.lz=fabs(Box.z1-Box.z0);		
		Box.hlx=Box.lx/2.0;
		Box.hly=Box.ly/2.0;
		Box.hlz=Box.lz/2.0;
	}
	
	if(p.exist("PeriodicBoundary")) {
		// leave it if it is already set!
		SetBoundaryCondition(p.int_value("PeriodicBoundary"), false);
	}
	
	if(!ia) {
		
		if(mat_file=="") {
			mdassert(p.exist("Material"), 
				   "can not find Material tag in import file ("+fname+
				   ") nor material definition in the parameter list");
			mat_file=p.string_value("Material");
		}
		
		AddAtom(new AtomContainer(mat_file))->Import(fname);
	}

	return this;
}

void MDSystem::AcceptSignal(int signo) {
	blog("accepting signal "+as_string(signo));
	struct sigaction act;
	memset(&act,0,sizeof(act));
	act.sa_handler=&SignalHandler;
	sigaction(signo,&act,NULL);
}

void MDSystem::CheckInterruption(){
	InterruptFlag=(int)sig_int_accept;
	if(InterruptFlag) {
		std::cerr << "#### INTERRUPT SIGNAL: " << InterruptFlag << " ####\n";
		blog("interrupt signal caught: "+as_string(InterruptFlag));
		switch(InterruptFlag) {
			case SIGTERM:
				OnInterruptTERM(); // default: dump&terminate
				break;
			
			case SIGINT:
				OnInterruptINT(); // default: dump&terminate
				break;
			
			case SIGUSR1:
				OnInterruptUSR1(); // default: dump&terminate
				break;

			case SIGUSR2: 
				OnInterruptUSR2(); // default: dump&continue (non-terminating)
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
	InterruptFlag=0;
}

void MDSystem::ActivateGadget(string gname) {
	MDGadget* g=SearchGadget(gname);
	if(g)g->Activate();
	blog("gadget ("+gname+","+g->get_type()+
	     ") Activated step="+as_string(Step));
}

void MDSystem::DeactivateGadget(string gname) {
	MDGadget* g=SearchGadget(gname);
	if(g) {
		g->Deactivate();
		blog("gadget ("+gname+","+g->get_type()+
		     ") Deactivated step="+as_string(Step));
	}
}

AtomContainer* MDSystem::SearchContainer(string name) {
	if(name==get_name()) return this;
	
	// Atom group has priority...
	for(int i=0;i<(int)SystemAtomGroups.size();i++)
		if(name==SystemAtomGroups[i]->get_name())
			return SystemAtomGroups[i];
	
	for(int i=0;i<(int)SystemAtoms.size();i++)
		if(name==SystemAtoms[i]->get_name())
			return SystemAtoms[i];
	
	die("can not find target atom '"+name+"'");
	return NULL; // avoids warning..
}

void MDSystem::DeleteContainer(string name) {
  for(int i=0;i<(int)SystemAtoms.size();i++) {
    if(SystemAtoms[i]->get_name()==name) {
      delete SystemAtoms[i];
      SystemAtoms.erase(SystemAtoms.begin()+i);
      return;
    }
  }
}

