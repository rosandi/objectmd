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
#include <omd/config.h>
#include <omd/container.h>
#include <omd/atomgroup.h>
#include <omd/unit.h>

using std::vector;
using std::string;

namespace omd {

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
#define CONTINUE_MODE 1
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
class Modify;
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
 * CONTINUE_MODE. In the continue mode, the class will read a previously saved
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
 * MDSystem will add a default iterator, if none is added in CreateGadget()
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
    
    int Mode; /**<the running mode: NORMAL_MODE|CONTINUE_MODE|TEST_MODE**/
    
    int AtomID;
    int GroupID;
    int ModifyID;
    int DetectorID;
    bool Unificated;
    bool Enumerated;
    bool stopped;
    
    
    string ParameterFilename;
    string BinaryFilename;
    string OutputDirectory;
    
  public:
    
    bool silent_mode;
    int  print_every;
    
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
    
    MDIterator*            Iterator;
    MDIntegrator*          Integrator;
    vector<Modify*>        Modifies;
    vector<Detector*>      Detectors;
    vector<AtomContainer*> SystemAtoms;
    vector<AtomGroup*>     SystemAtomGroups;
    vector<DataSlot*>      MessageSlots;
    vector<DataSlot*>      RestartVars;
    
    time_t    SimBeginTime, SimEndTime;
    double    SimWallTime; 
    
    int      PBoundary;
    int    TotalAtom;
    
    double Energy;
    double Kinetic;
    double Virial;
    double Potential;
    double BasePotential;
    
    int Step;
    double ElapsedTime;
    double MaxTime;
    double SqrMaxVelocity;
    
  protected:
    
    virtual void CreateSystem();
    virtual void CreateGadget()=0;
    virtual void SystemSetting(){} // settings may be done via parameters
    virtual void CreateGroup(){} // atom grouping, called after adjust_system/system_setting
    
    // Default=Do nothing functions. reserved for user application.
    virtual void PreCreation(){}
    virtual void PostCreation(){}
    virtual bool ConfirmResume(){return true;}
    virtual void Scheduler(){}
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
    virtual bool CheckRun(){
      if(stopped) return false;
      return((ElapsedTime<=MaxTime)&&(InterruptFlag==0));
    }
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
    virtual void Init() {}; // reserved for user initialization (before Initiate())
                            // after ReadParameter
    
    virtual void ReadParameter();
    virtual int Run(int mode=NORMAL_MODE);
    
    MDIntegrator* SetIntegrator(MDIntegrator* itg);
    Detector* AddDetector(Detector* Detc);
    Modify* AddModify(Modify* Cond);
    AtomContainer* AddAtom(AtomContainer* Atm);
    AtomGroup* AddAtomGroup(string group_name);
    
    class ForceKernel* AddForce(class ForceKernel* Force);
    class ForceKernel* AddForce(class ForceKernel* Force, string from, string to);
    
    class ForceKernel* AddForce(class ForceKernel* Force, string target){
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
    
    bool OnTime(double tm);
    bool OnStep(int step);
    
    virtual void PrintTime(ostream& ost);
    virtual void PrintHeader(ostream& ost);
    virtual void ExecuteModifies(int modype);
    virtual void ExecuteDetectors();
    
    virtual AtomContainer* Save(string binname, string mode="a"); // mode is not applicable...
    virtual string GetBinaryFilename();
    virtual void   SetBinaryFilename(string filename);
    virtual void   SaveVariables(FILE* fl);
    virtual void   LoadVariables(FILE* fl);
    virtual void   SaveSimulationConfig(string binfile);
    virtual void   SaveSimulation(string binfile="");
    virtual void   LoadSimulation(string LoadFromFile);
    
    virtual void BorderOffset(double dx, double dy, double dz);
    virtual void BorderOffset(double dr){BorderOffset(dr,dr,dr);}
    virtual void SetBox(double x0, double y0, double z0, double x1, double y1, double z1);
    virtual MDIterator* GetIterator(){return Iterator;}
    virtual MDIntegrator* GetIntegrator(){return Integrator;}
    
    /** synchronizing data **/
    virtual void SyncData(int syncmode){}
		
    virtual double GetMaxCutRadius();
    virtual double GetMaxVelocity();
    virtual double GetMass(int idx);
    virtual double GetMass(Atom& a);
    virtual double GetMass(Atom* a);
    virtual double GetZ(int idx);
    
    virtual void MeasurePotential();
    virtual void MeasureKinetic();
    
    virtual void CheckBoundary();
    
    /** Corrects distances according to the boundary condition. **/
    virtual void BoundaryCorrectDistances(double& dx, double& dy, double& dz);
    
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
    
    // parameter filename. call only before run()
    void SetParameterFile(string cfgfile) {
      ParameterFilename=cfgfile;
    }
    
    void LoadParameterFile();
		
    // output directory
    void SetOutputDirectory(const string outdir) {
      OutputDirectory.assign(outdir);
    }
    
    string GetOutputDirectory(){
      return OutputDirectory;
    }
    
    DataSlot* GetMessageSlot(string slotlabel);
    
    virtual void ArrangeMessageSlots();
    int GetContainerID(string name);	
    
    /** search gadget class by its name (Modify/Detector) **/
    virtual MDGadget* SearchGadget(string name);
    bool GadgetExist(MDGadget* gad);
    
    void DeleteContainer(string name);
    virtual AtomContainer* SearchContainer(string name);
    virtual AtomContainer* Import(string fname);
    
    int GetMode(){return Mode;}
		
    void ResetSimulationTime(){Step=0; ElapsedTime=0;}
    
    /** sets MaxTime only when it is not already done previously **/
    void SetMaxTime(double maxtime) {if(MaxTime<0)MaxTime=maxtime;}
    
    /** sets boundary condition only when is not already done previously **/
    void SetBoundaryCondition(int pbc, bool force=false) {
      if(force) PBoundary=pbc;
      else if(PBoundary<0) PBoundary=pbc;
    }
    
    virtual int GetLocalAtomNumber(){return GetNAtom();}
    virtual int GetTotalAtomNumber(){return GetNAtom();}
    
    double GetBasePotential() {return BasePotential;}
    void ActivateGadget(string gname);
    void DeactivateGadget(string gname);
    void Stop(){stopped=true;}
    
    // short cuts...
    AtomContainer* AddAtom(string name, AtomContainer* Atm){return AddAtom(Atm)->SetName(name);}
    
  };
  
}

#endif
