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
 * Force kernel
 *
*/


#ifndef _FORCE_KERNEL_H_
#define _FORCE_KERNEL_H_

#include <omd/gadget.hpp>
#include <omd/conditioner.hpp>

//----------------ForceKernel---------------------//

//! The kernel of interaction force

/** 
 * 
 * This class calculates of the interaction forces between atoms
 * in the simulation system. The function Compute() is called by the Integrator
 * to compute the potential and force value between two atoms, represented
 * by "at" and "to". Two integer variables, A and B, are reserved for storing
 * the ids of interacting atom types.
 * 
 * Using this class, the interaction potential is implemented as one-to-one
 * interaction, no loop needed. The force calculation loop is taken care by 
 * the Integrator class, via its Iterator.
 * 
 * In order to make a potential class, the Compute() function (abstract) must
 * be implemented. This function have to return the calculated values, through
 * ReturnForce() function, to make sure that the virial and force components 
 * are calculated correctly.
 * 
 * Since the initialization is done by the Integrators
 * Init() function, it is guarantied that the Integrator is ready and accessible
 * in this stage, via MDSystem->Forces variable.
 * 
 * The Correction() function is called after the force-loop. This function may
 * be used to perform some operations that can be taken out from the force
 * calculation equation (e.g. multiplication factor, etc).
 *
*/
 
class ForceKernel:public MDGadget {
	friend class MDIntegrator;
	friend class MDSystem;

protected:
	int AtomTypeA, AtomTypeB;

	OMD_FLOAT CutRadius;
	OMD_FLOAT CutRadiusSqr;
	MDGadget* force_eval;

public:
	
	ForceKernel(){CutRadiusSqr=CutRadius=0.0; AtomTypeA=AtomTypeB=0;force_eval=NULL;}
	virtual ~ForceKernel(){}

	virtual void Init(MDSystem* WorkSys){
		MDGadget::Init(WorkSys);
		if(force_eval){
			assert(WorkSys->GadgetExist(force_eval), 
			       "force evaluator does not exist");
		}
	}
	
	virtual void ClearAccumulators(){}
	
	virtual void ReturnForce(Atom &at, Atom &to,
							 OMD_FLOAT dx, OMD_FLOAT dy, OMD_FLOAT dz,
							 OMD_FLOAT fr, OMD_FLOAT pot) {
		// fr=f/r
		
		if(force_eval) { // if exist: force evaluator
			force_eval->EvaluateForce(at,to,dx,dy,dz,fr,pot,this);
		}

		OMD_FLOAT fx=dx*fr,fy=dy*fr,fz=dz*fr;
		OMD_FLOAT vir=(fx*dx+fy*dy+fz*dz);
		vir*=0.5;
		pot*=0.5;
		at.virial+=vir;
		to.virial+=vir;
		at.potential+=pot;
		to.potential+=pot;
		at.fx+=fx;at.fy+=fy;at.fz+=fz;
		to.fx-=fx;to.fy-=fy;to.fz-=fz;
	}

	/*
	 * both reference and index are usefull!
	 */

	virtual void Compute(Atom &at, Atom &to) {die("Compute(Atom&,Atom&) is not implemented!!");}
	virtual void Compute(int at, int to) {Compute(Atoms(at), Atoms(to));}

	/**
	 * makes sure that the for calculate the right interacting atoms
	 */
	
	virtual void CheckCompute(int at, int to, int atid, int toid) {
		if(atid==AtomTypeA&&toid==AtomTypeB){Compute(at,to);return;}
		if(atid==AtomTypeB&&toid==AtomTypeA){Compute(at,to);return;}
		die("wrong type id of interacting atoms ids: "+as_string(atid)+" - "+as_string(toid));
	}

	virtual void   Correction() {}
	ForceKernel*   SetAtomID(int a, int b) {AtomTypeA=a; AtomTypeB=b; return this;}
	ForceKernel*   SetEvaluator(Conditioner* feval){force_eval=feval; return this;}

	bool CheckAtomID(int a, int b) {
		// commutative...
		if(a==AtomTypeA&&b==AtomTypeB) return true;
		if(a==AtomTypeB&&b==AtomTypeA) return true;
		return false;
	}
	                   
	ForceKernel*   SetCutRadius(OMD_FLOAT c) {CutRadius=c; return this;}
};

#endif
