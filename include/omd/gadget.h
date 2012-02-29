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
 *    MDGadget --> Base class of all gadgets
*/

#ifndef _MD_GADGET_H_
#define _MD_GADGET_H_

#include <cmath>
#include <omd/config.h>
#include <omd/system.h>
#include <omd/dataslot.h>

namespace omd {

//-------------MDGadget--------------//

/** 
 * @ingroup essential
 * @brief Base class for all gadgets
 *
 * MDGadget is the root class of all the plugable classes in 
 * Object-MD. These classes is called gadgets. To every gadgets is assign a Target
 * container. The access to atom structure (via Atoms()) is redirected to Target->Atoms().
 * The target container is initialized in MDGadget::Init() by defining TargetName in the
 * constructor or in ReadParameter().
 * 
 * Class type convention: lower-case, no space.
 * 
 * to explain:
 *   - System and Target
 *   - activation
 *   - claiming auxilary variable
 *   - claiming flag bit
 *   - message slot
 *   - common functions
 *   - restart variable
 * 
 */

  class MDGadget: public MDClass {
    
  protected:
    MDSystem* System;
    AtomContainer* Target;
    string TargetName;
    int  Active;
    int  ActiveCode;
    bool Ready;
    MDIterator* Iterator;
    MDIntegrator* Integrator;
    int AuxIdx;
    MDUnit* Unit;
    
    int NCalls;
    ParamHandler *SysParam;
    
  public:
    
    bool IsReady;
    
    MDGadget();
    int GetContainerID(string name);
    
    AtomContainer* SearchTarget(string name, bool strict=true);
    MDGadget* SearchGadget(string name, bool strict=true);
    MDGadget* SearchGadgetType(string type, bool strict=true);
    
    /**  
     * @brief Initialize the gadget
     * 
     * Atoms points to the same array and index of WorkSys->Atoms.
     * The iterator is set to the system's iterator. This function reads and check
     * gadgets' parameter by invoking ReadParameter and CheckParameter function if
     * reimplemented in the child gadget class.
     * 
     */
    
    virtual void Init(MDSystem* WorkSys);
    
    /**
     * Atom reference version of IterationNode()/PreIterationNode().
     * This function should be the one implemented by a descendant, if
     * it does not need the atom index information.
     */
    
    virtual bool PreIterationNode(Atom &at){return true;}
    
    /**
     The iteration node that is called by MDIterator::Iterate() function,
     in the deeper nest of the iteration loop.
     */
    
    virtual void IterationNode(Atom &at,Atom &to){}
    
    /**
     * Some classes need the index of the atoms. This version of 
     * IterationNode()/PreIterationNode()
     * provides integer parameter as index. By default the function calls
     * The atom reference version.
     */
    
    virtual bool PreIterationNode(int at){
    	return PreIterationNode(Atoms(at));
    }
    
    
    virtual void IterationNode(int at,int to){
    	IterationNode(Atoms(at), Atoms(to));
    }
    
    Atom&  Atoms(int idx);
    Atom*  AtomPtr(int idx);
    int   GetNAtom();
    double GetMass(int idx);
    double GetMass(Atom &a);
    double GetMass(Atom *a);
    double GetZ(int idx);
    double GetZ(Atom &a);
    double GetTimeStep();
    double GetElapsedTime();
    AtomContainer* GetTarget(){return Target;}
    void SetTarget(string tname){TargetName=tname;}
    
    MDGadget* Activate() {Active|=ActiveCode; return this;}
    MDGadget* Deactivate() {Active=0; return this;}
    MDGadget* SetName(string name){set_name(name);return this;}
    MDGadget* SetSystem(MDSystem* WorkSys){System=WorkSys; return this;}
    
    int ClaimFlagBit();
    int ClaimAuxVariable(bool printable, const char* tag, const char* sformat=NULL);
    int ClaimAuxVariable(){return ClaimAuxVariable(false,NULL,NULL);}
    double& AuxVariable(int i);
    
    bool OnTime(double tm);
    bool OnStep(int step);
    
    int IsActive(int Code=0);
    
    /**
     * Calculates the square distance of two atoms, identified by their index.
     * System is responsible to correct the distance, regarding the applied
     * boundary condition. The function BoundaryCorrectDistances() is called
     * after calculating distance elements.
     */
    
    double CalcSqrDistance(Atom &at, Atom &to, 
                           double &dx, double &dy, double &dz, bool check=true);
    
    double CalcSqrDistance(int at, int to, 
                           double &dx, double &dy, double &dz, bool check=true){
      return CalcSqrDistance(Atoms(at), Atoms(to), dx, dy, dz);
    }
    
    double CalcSqrDistance(Atom &at, Atom &to, bool check=true) {
      double dx, dy, dz;
      return CalcSqrDistance(at, to, dx, dy, dz, check);
    }
    
    double CalcSqrDistance(int at, int to, bool check=true) {
      return CalcSqrDistance(Atoms(at), Atoms(to), check);
    }
    
    double CalcDistance(Atom &at, Atom &to, bool check=true){
      return sqrt(CalcSqrDistance(at, to, check));
    }
    
    
    double CalcDistance(int at, int to, bool check=true) {
      return sqrt(CalcSqrDistance(Atoms(at), Atoms(to), check));
    }
    
    virtual void EvaluateForce        // inserted in force loop (at the begining of ReturnForce)
    (Atom& a, Atom& b, 
     double dx, double dy, double dz, 
     double fr, double pot, ForceKernel* fkernel) {}
    
    void SetUnit(MDUnit* unit){Unit=unit;}
    
    DataSlot* RegisterMessageSlot(DataSlot* slot);
    
    void RestartVariable(string tag, double &val);
    void RestartVariable(string tag, int &val);
    
    void PrintInfo(ostream& ost) {ost<<"id."<<id<<" "<<get_name()<<"\n";}
    
    int NumberOfCalls(){return NCalls;}
    void SyncData(int mode) {System->SyncData(mode);}
    void SyncData(int mode, int aidx) {System->SyncData(mode|aidx);}
    
    virtual void ReadParameter(){}
    virtual bool CheckParameter(){return true;}
    
  };
  
}

#endif
