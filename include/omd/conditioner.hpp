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
 * Base class of all conditioners
 *
*/

#ifndef _CONDITIONER_H_
#define _CONDITIONER_H_

#include <omd/gadget.hpp>

#define COND_PRE_INTEGRATION  1
#define COND_PRE_CALCULATION  2
#define COND_FORCE_MODIFIER   4
#define COND_POST_INTEGRATION 8

//------------------Conditioner-----------------//
/**
 * @ingroup conditioner
 * @brief Conditioner for modifying system atoms
 *
 * The conditioner class is used to make modification on system while the
 * simulation is running. The following are the member functions that can be 
 * implemented in the child class,
 * 	 - PreIntegration(), called before system integration. 
 *   - PreCalculation(), called by the integrator class before the force 
 *     calculation loop.
 *   - ForceModifier(), called by the integrator after the force 
 *     calculation loop.
 *   - PostIntegration(), called after the integration, before 
 *     executing Detectors.
 * 
 * The Conditioner::Active variable identifies the type of conditioner:
 * 
 *	- 1 --> Pre-integration conditioner
 *	- 2 --> Pre-calculation canditioner
 *	- 4 --> Force-modifier
 *	- 8 --> Post-integration conditioner
 * 
 * Zero value is always mean inactive. The bits of this variable is be examined
 * in the Conditioner::Execute() function.
 *
*/

class Conditioner: public MDGadget {

public:
	Conditioner() {Active=15; Target=NULL; TargetName="";}
	
	virtual ~Conditioner() {blog("cleaning up", LOGDESTROY);}
	
	void PrintInfo(ostream& ost) {
		ost<<"id."<<id<<" "<<get_name()<<" -- target: "<<Target->get_name()<<std::endl;
	}
	
	void SetConditionerType(int condtype) {Active=ActiveCode=condtype;} 

    virtual void PreIntegration()  {} // before integration
    virtual void PreCalculation()  {} // before force calculation loop
    virtual void ForceModifier()   {} // after force calculation
    virtual void PostIntegration() {} // after integration, before detection
    virtual void EvaluateForce        // inserted in force loop (before ReturnForce)
    		(Atom &a, Atom &b, OMD_FLOAT dx, OMD_FLOAT dy, OMD_FLOAT dz, OMD_FLOAT fr, OMD_FLOAT pot) {}

    virtual void Execute(int condtype) {
		if(!Active) return; // Active stors the type of conditioner...
    	if(condtype&COND_PRE_INTEGRATION&Active)  PreIntegration();
    	if(condtype&COND_PRE_CALCULATION&Active)  PreCalculation();
    	if(condtype&COND_FORCE_MODIFIER&Active)   ForceModifier();
    	if(condtype&COND_POST_INTEGRATION&Active) PostIntegration();
		NCalls++;
	}

};

/** @brief Pre integration conditioner **/
class Pre_Conditioner:public Conditioner {
	public:
		Pre_Conditioner(){Active=ActiveCode=COND_PRE_INTEGRATION;}
};

/** @brief Pre force calculation loop conditioner **/
class Calc_Conditioner:public Conditioner {
	public:
		Calc_Conditioner(){Active=ActiveCode=COND_PRE_CALCULATION;}
};

/** @brief Force modifier container **/
class Force_Conditioner:public Conditioner {
	public:
		Force_Conditioner(){Active=ActiveCode=COND_FORCE_MODIFIER;}
};

/** @brief Post integration conditioner **/
class Post_Conditioner:public Conditioner {
	public:
		Post_Conditioner(){Active=ActiveCode=COND_POST_INTEGRATION;}
};

#endif
