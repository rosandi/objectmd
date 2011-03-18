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
 *      MDSystem
 *
 *      This header file contains the definition of
 *      The main class in the library (MDSystem).
 *
*/

#ifndef _SIMSYSTEM_H_
#define _SIMSYSTEM_H_

#include <iostream>
#include <string>
#include <vector>
#include <cstdio>
#include <ctime>
#include <signal.h>
#include <omd/config.hpp>
#include <omd/container.hpp>
#include <omd/atomgroup.hpp>
#include <omd/unit.hpp>

using std::vector;
using std::string;

#define NONPERIODIC        0
#define PERIODIC_X         1
#define PERIODIC_Y         2
#define PERIODIC_Z         4
#define PERIODIC           7

#define OMD_EXIT_FAILURE   EXIT_FAILURE
#define OMD_EXIT_SUCCESS   EXIT_SUCCESS
#define OMD_EXIT_SUSPEND   254

#define IS_PERIODIC_X (PBoundary&PERIODIC_X)
#define IS_PERIODIC_Y (PBoundary&PERIODIC_Y)
#define IS_PERIODIC_Z (PBoundary&PERIODIC_Z)

#define NORMAL_MODE  0
#define RESTART_MODE 1
#define STATIC_MODE  2

#define SYNC_ALL       0x0F00
#define SYNC_POSITION  0x0100
#define SYNC_VELOCITY  0x0200
#define SYNC_FORCE     0x0400
#define SYNC_AUX       0x0800
#define SYNC_SPACE     0x1000

#define WEST_EDGE   (Box.x0-EDGE_TOLE)
#define EAST_EDGE   (Box.x1+EDGE_TOLE)
#define SOUTH_EDGE  (Box.y0-EDGE_TOLE)
#define NORTH_EDGE  (Box.y1+EDGE_TOLE)
#define BOTTOM_EDGE (Box.z0-EDGE_TOLE)
#define TOP_EDGE    (Box.z1+EDGE_TOLE)

class MDGadget;
class DataSlot;
class MDIntegrator;
class Detector;
class Conditioner;
class MDIterator;

//----------------------MDSystem-----------------------------//

/**
 * @brief MD Simulation Manager
 *
 * MDSystem is the manager of Object-MD simulation program.
 * The class manages atoms and their containers. Here the main loop of the 
 * program takes place. In the initialization stage all gadgets must take the 
 * reference to the instance of this class.
 *
 * This class is instantiate in the application to form a simulation program.
 * The application must implement all abstract methods. The followings are
 * the abstract function:
 *     - CreateSystem(), creates the atom structures. Atom containers 
 *       are inserted to the system using AddContainer().
 *     - CreateGadget(), adds gadgets into the system.
 *     - SystemSetting(), adjustment and simulation setting. The followings are
 *       mandatory:
 *            - Boundary condition, PBoundary. Options: PERIODIC_X, 
 *              PERIODIC_Y, PERIODIC_Z, PERIODIC, NOPERIODIC, or combination
 *              between those macros (use bit-wise "or" |)
 *            - Maximum simulation time, MaxTime
 *            - etc.
 * 
 * The functions BeforeRun() and AfterRun() are reserved to define processes
 * that is executed before and after the simulation loop.
 * 
 * The main loop is invoked by calling Run() function. The function returns a
 * status code after running the simulation (OMD_EXIT_FAILURE, OMD_EXIT_SUCCESS,
 * OMD_EXIT_SUSPEND). Currently, two run modes are implemented: NORMAL_MODE and
 * RESTART_MODE. In the restart mode, the class will read a previously saved
 * Object-MD binary file, and continue the simulation.
 * 
 * MDSystem claims 3 flag bits from the Atom structure used as 
 * active, inside box, and ghost flags. The other bits may be claimed by the
 * gadgets, via the MDSystem::ClaimFlagBit(). MDClass has a function with
 * the same name, which automatically redirect the call from gadget descendants.
 * Besides flag, a gadget can also claim aux variables of atoms using
 * ClaimAuxVariable() function, implemented also in MDGadget.
 * 
 * Every step of the program loop, MDSystem prints information to stdout by 
 * PrintMessages(), one line every step. An information from gadgets can be 
 * printed in the line using message slot mechanism. A message slot (DataSlot) is
 * a data register, that points to a single data in the gadget. A gadget have
 * to register its message slot using RegisterMessageSlot() function.
 * 
 * Gadgets have their own message slot register. The registers from all gadgets
 * are gathered by MDSystem on the initialization stage (InitGadget). So, a
 * gadget may register slots before or in its Init() function. (FIXME!)
 * 
 * MDSystem will add a default iterator conditioner, if none is added in CreateGadget()
 * function. The default iterator is of type MDIterator (a slow, traditional iteration
 * loop).
 *
 * On creation, the file "omd-parameter" is loaded if exist.
*/

class MDSystem:public AtomContainer {

protected:
	friend class MDGadget;

	int FlagBitUsed;     /**<the mask of used flag bit**/
	int AuxVariableUsed; /**<the bit mask of used variables**/
    int  InterruptFlag;   /**<system interrupt flag**/

	int Mode; /**<the running mode: NORMAL_MODE|RESTART_MODE|TEST_MODE**/

	int AtomID;
	int GroupID;
	int ConditionerID;
	int DetectorID;
	bool BoxImport;
	bool Unificated;
	bool Enumerated;

    MDIterator* Iterator;

	string RestartFileName;
	string OutputDirectory;
	
// RUNTIME FLAGS
	bool silent_mode;

public:

    int* Argc;
    char*** Argv;
    int  ExitCode;
    
	vector<string>  AuxUser;
	vector<string> SAuxFormat;
	vector<string>  AuxNameTag;
	vector<string>  FlagUser;
	
	bool PrintableAux[MAXAUXVAR];
	char** AuxFormat;

	MDUnit* Unit;

    MDIntegrator*          Integrator;
    vector<Conditioner*>   Conditioners;
    vector<Detector*>      Detectors;
    vector<AtomContainer*> SystemAtoms;
    vector<AtomGroup*>     SystemAtomGroups;
	vector<DataSlot*>      MessageSlots;
	vector<DataSlot*>      RestartVars;

	time_t    SimBeginTime, SimEndTime;
	OMD_FLOAT    SimWallTime; 

    int      PBoundary;
	int    TotalAtom;
	
    OMD_FLOAT Energy;
    OMD_FLOAT Kinetic;
    OMD_FLOAT Virial;
    OMD_FLOAT Potential;
    OMD_FLOAT BasePotential;
    
    int Step;
    OMD_FLOAT ElapsedTime;
    OMD_FLOAT MaxTime;
	OMD_FLOAT SqrMaxVelocity;

protected:
	
	// Abstract functions
    virtual void CreateSystem()=0;
    virtual void CreateGadget()=0;
    virtual void SystemSetting(){} // settings may be done via parameters

	// Default=Do nothing functions. reserved for user application.
	virtual void PreCreation(){}
    virtual void PostCreation(){}
    virtual bool ConfirmResume(){return true;}
	virtual void Scheduller(){}
	virtual void InlineFunction(){}
	
	// interrupt signal handlers
	virtual void OnInterruptTERM();
	virtual void OnInterruptINT();
	virtual void OnInterruptUSR1();
	virtual void OnInterruptUSR2();
	virtual void CheckInterruption();
	
	// Normal functions
	virtual void CreationFunction();
    virtual void AdjustSystem();
    virtual void EnumerateAtoms(bool redo=false);
    virtual void UnificateAtoms();
	virtual void Initiate();
	virtual void InitGadgets();
	virtual void PrintMessages(ostream& ost);
	virtual SysBox& CalcBox();

	virtual void CheckBeforeRun();
    virtual void BeforeFirstRun(){}
	virtual void BeforeRun(){}
    virtual void AfterRun(){}
	virtual bool CheckRun(){return((ElapsedTime<=MaxTime)&&(InterruptFlag==0));}
	virtual void FirstRun();
	virtual void RunKernel();

// run modes
	virtual void RunNormal();
	virtual void RunRestart();
	virtual void RunStatic();

public:      

	MDSystem(int &argc, char** &argv);
	MDSystem();
	
    virtual ~MDSystem();
    void SystemInit();

    virtual void ReadParameters();
	virtual int Run(int mode=NORMAL_MODE);

	MDIntegrator* SetIntegrator(MDIntegrator* itg);
	Detector* AddDetector(Detector* Detc);
	Conditioner* AddConditioner(Conditioner* Cond);
	AtomContainer* AddAtom(AtomContainer* Atm);
	AtomGroup* AddAtomGroup(string group_name);
	
	class ForceKernel* AddForce(class ForceKernel* Force);
	class ForceKernel* AddForce(class ForceKernel* Force, const char* from, const char* to);

	class ForceKernel* AddForce(class ForceKernel* Force, const char* target){
		return AddForce(Force, target, target);
	}
	
	void InsertDataHeader(ofstream& fl) {} // Maybe needed later
	void ChangeAtomID(int idx, int NewID);
	void ChangeAtomID(int start, int end, int NewID);
	
	virtual void PrintContainerInfo(ostream& ost);
	virtual void PrintGadgetInfo(ostream& ost);
	virtual void PrintSystemInfo(ostream& ost);
	virtual void PrintInfo(ostream& ost);
	virtual void PrintInfo(string fname);
	
	bool OnTime(OMD_FLOAT tm);
	bool OnStep(int step);
	
	virtual void PrintTime(ostream& ost);
	virtual void PrintHeader(ostream& ost);
	virtual void ExecuteConditioners(int contype);
	virtual void ExecuteDetectors();

	virtual AtomContainer* Save(string binname, string mode="a"); // mode is not applicable...
	virtual string GetRestartFilename();
	virtual void   SetRestartFilename(string filename);
	virtual void   SaveVariables(FILE* fl);
	virtual void   LoadVariables(FILE* fl);
	virtual void   SaveSimulationConfig(string binfile);
	virtual void   SaveSimulation(string binfile="");
	virtual void   LoadSimulation();

	virtual void BorderOffset(OMD_FLOAT dx, OMD_FLOAT dy, OMD_FLOAT dz);
	virtual void BorderOffset(OMD_FLOAT dr){BorderOffset(dr,dr,dr);}
	virtual void SetBox(OMD_FLOAT x0, OMD_FLOAT y0, OMD_FLOAT z0, OMD_FLOAT x1, OMD_FLOAT y1, OMD_FLOAT z1);
	virtual MDIterator* GetIterator(){return Iterator;}
	virtual MDIntegrator* GetIntegrator(){return Integrator;}

	/** synchronizing data **/
	virtual void SyncData(int syncmode){}
		
	virtual OMD_FLOAT GetMaxCutRadius();
	virtual OMD_FLOAT GetMaxVelocity();
	virtual OMD_FLOAT GetMass(int idx);
	virtual OMD_FLOAT GetMass(Atom& a);
	virtual OMD_FLOAT GetMass(Atom* a);
	virtual OMD_FLOAT GetZ(int idx);
	
	virtual void MeasurePotential();
	virtual void MeasureKinetic();

	virtual void CheckBoundary();

	/** Corrects distances according to the boundary condition. **/
	virtual void BoundaryCorrectDistances(OMD_FLOAT& dx, OMD_FLOAT& dy, OMD_FLOAT& dz);
	
	virtual void ErrorHandler(const char* errst);
	virtual int  ClaimFlagBit(MDClass* user, string sinfo="");
	virtual int  ClaimAuxVariable(MDClass* user, bool printable=false, const char* tag=NULL, const char* sformat=NULL);
	                             
	virtual int GetFlagBitMask(const char* usagecode);
	virtual int GetTotalAtom(){return GetNAtom();}

	void SetUnit(
		const char* stime, const char* slength, const char* smass, 
        const char* sforce, const char* senergy, 
        const char* stemp, const char* spress)
	{
		if(!Unit) Unit=new MDUnit;
		Unit->SetUnit(stime,slength,smass,sforce,senergy,stemp,spress);
	}
	
	void SetUnit(MDUnit* unit){if(Unit)delete Unit; Unit=unit;}
	
	// compatibility reason...
	void SetArgument(int &argc, char** &argv);
		
	// react to this signal
	void AcceptSignal(int signo);	
		
	// output directory
	virtual void SetOutputDirectory(const string outdir) {
		OutputDirectory.assign(outdir);
	}
	
	virtual string GetOutputDirectory(){
		return OutputDirectory;
	}
	
	DataSlot* GetMessageSlot(string slotlabel);
	
	virtual void ArrangeMessageSlots();
	int GetContainerID(string name);	

	/** search gadget class by its name (Conditioner/Detector) **/
	virtual MDGadget* SearchGadget(string name);
	bool GadgetExist(MDGadget* gad);
	
	virtual AtomContainer* Import(string fname);

	int GetMode(){return Mode;}
	
	stage_type stage;
	
	void ResetSimulationTime(){Step=0; ElapsedTime=0;}

	/** sets MaxTime only when it is not already done previously **/
	void SetMaxTime(OMD_FLOAT maxtime) {if(MaxTime<0)MaxTime=maxtime;}

	/** sets boundary condition only when is not already done previously **/
	void SetBoundaryCondition(int pbc) {if(PBoundary<0)PBoundary=pbc;}

	virtual int GetLocalAtomNumber(){return GetNAtom();}
	virtual int GetTotalAtomNumber(){return GetNAtom();}
	
	void ActivateGadget(string gname);
	void DeactivateGadget(string gname);

};

#endif
