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
#include <omd/config.hpp>
#include <omd/system.hpp>
#include <omd/dataslot.hpp>

//-------------MDGadget--------------//

/** 
 * @ingroup essential
 * @brief Base class for all gadgets
 *
 * MDGadget is the root class of all the plugable classes in 
 * Object-MD. These classes is called gadgets.
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
	int  Active;
	int  ActiveCode;
	string TargetName;
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
    OMD_FLOAT GetMass(int idx);
    OMD_FLOAT GetMass(Atom &a);
	OMD_FLOAT GetMass(Atom *a);
    OMD_FLOAT GetZ(int idx);
    OMD_FLOAT GetZ(Atom &a);
    OMD_FLOAT GetTimeStep();
    OMD_FLOAT GetElapsedTime();
    AtomContainer* GetTarget(){return Target;}

	MDGadget* Activate() {Active|=ActiveCode; return this;}
	MDGadget* Deactivate() {Active=0; return this;}
	MDGadget* SetSystem(MDSystem* WorkSys){System=WorkSys; return this;}
	int ClaimFlagBit();
	int ClaimAuxVariable(bool printable, const char* tag, const char* sformat=NULL);
	int ClaimAuxVariable(){return ClaimAuxVariable(false,NULL,NULL);}
	OMD_FLOAT& AuxVariable(int i);

	bool OnTime(OMD_FLOAT tm);	
	int IsActive(int Code=0);

	/**
	 * Calculates the square distance of two atoms, identified by their index.
	 * System is responsible to correct the distance, regarding the applied
	 * boundary condition. The function BoundaryCorrectDistances() is called
	 * after calculating distance elements.
	 */
	
	OMD_FLOAT CalcSqrDistance(Atom &at, Atom &to, 
						   OMD_FLOAT &dx, OMD_FLOAT &dy, OMD_FLOAT &dz, bool check=true);
	
	OMD_FLOAT CalcSqrDistance(int at, int to, 
						   OMD_FLOAT &dx, OMD_FLOAT &dy, OMD_FLOAT &dz, bool check=true){
		return CalcSqrDistance(Atoms(at), Atoms(to), dx, dy, dz);
	}

	OMD_FLOAT CalcSqrDistance(Atom &at, Atom &to, bool check=true) {
		OMD_FLOAT dx, dy, dz;
		return CalcSqrDistance(at, to, dx, dy, dz, check);
	}

	OMD_FLOAT CalcSqrDistance(int at, int to, bool check=true) {
		return CalcSqrDistance(Atoms(at), Atoms(to), check);
	}
	
	OMD_FLOAT CalcDistance(Atom &at, Atom &to, bool check=true){
		return sqrt(CalcSqrDistance(at, to, check));
	}
	

	OMD_FLOAT CalcDistance(int at, int to, bool check=true) {
		return sqrt(CalcSqrDistance(Atoms(at), Atoms(to), check));
	}

	void SetUnit(MDUnit* unit){Unit=unit;}
	
	DataSlot* RegisterMessageSlot(DataSlot* slot);
	
	void RestartVariable(string tag, OMD_FLOAT &val);
	void RestartVariable(string tag, int &val);
	
	void PrintInfo(ostream& ost) {ost<<"id."<<id<<" "<<get_name()<<"\n";}
	
	int NumberOfCalls(){return NCalls;}
	void SyncData(int mode) {System->SyncData(mode);}
	void SyncData(int mode, int aidx) {System->SyncData(mode|aidx);}

	virtual void ReadParameter(){}
	virtual bool CheckParameter(){return true;}

};

#endif
