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
 *       Atoms equation of motion integrator
 *
*/

#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

//--------------Integrator-----------------//

#include <omd/config.hpp>
#include <omd/gadget.hpp>
#include <omd/forcekernel.hpp>


/** 
 * @ingroup essential
 * @brief Base of integrator class
 * 
 * This class is the integrator base class. It is a concrete class, that
 * implements the Verlet algorithm. This class is essential in Object-MD.
 *
*/

class MDIntegrator: public MDGadget {
	OMD_INT ForcID;
	bool CurrentGhostFlag;

public:
	vector<ForceKernel*> ActForces;
	OMD_SIZET   NType;
	OMD_FLOAT TimeStep;

    //Put force pointers in matrix form
    ForceKernel* Forces[MAXATOMTYPE][MAXATOMTYPE];
    
    OMD_FLOAT MaxCutRadius;
   
                 MDIntegrator(OMD_FLOAT time_step=-1.0);
	virtual      ~MDIntegrator(); 
	virtual void Init(MDSystem* WorkSys);
	        bool Check();
	MDIntegrator* SetTimeStep(OMD_FLOAT time_step){TimeStep=time_step; return this;}
	MDIntegrator* AddForce(ForceKernel* Forc);
	
	void PrintInfo(ostream& ost);

    virtual void IterationNode(OMD_SIZET at, OMD_SIZET to);
    virtual void Integrate();
    virtual void Iterate();
	ForceKernel* GetForce(OMD_SIZET idx){return ActForces.at(idx);}	
	ForceKernel* GetForce(OMD_SIZET a, OMD_SIZET b) {return Forces[a][b];}
};

#endif
