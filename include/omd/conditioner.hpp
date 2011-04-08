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
 * in the Conditioner::Execute() function. The calling number (NCalls) will be increased
 * when the corresponding integrator function is executed. If a conditioner has 
 * multiple-type, NCalls will be increased in every conditioner function.
 * (NCalls/step=number_of_type)
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

    virtual void Execute(int condtype) {
		if(!Active) return; // Active stors the type of conditioner...
    	if(condtype&COND_PRE_INTEGRATION&Active){PreIntegration(); NCalls++;}
    	if(condtype&COND_PRE_CALCULATION&Active){PreCalculation(); NCalls++;}
    	if(condtype&COND_FORCE_MODIFIER&Active){ForceModifier(); NCalls++;}
    	if(condtype&COND_POST_INTEGRATION&Active){PostIntegration(); NCalls++;}
	}
	
	virtual Conditioner* SetName(string newname) {
		set_name(newname);
		return this;
	}

};

/** @brief Pre integration conditioner **/
class PreConditioner:public Conditioner {
	public:
		PreConditioner(){Active=ActiveCode=COND_PRE_INTEGRATION;}
};

/** @brief Pre force calculation loop conditioner **/
class CalcConditioner:public Conditioner {
	public:
		CalcConditioner(){Active=ActiveCode=COND_PRE_CALCULATION;}
};

/** @brief Force modifier container **/
class ForceConditioner:public Conditioner {
	public:
		ForceConditioner(){Active=ActiveCode=COND_FORCE_MODIFIER;}
};

/** @brief Post integration conditioner **/
class PostConditioner:public Conditioner {
	public:
		PostConditioner(){Active=ActiveCode=COND_POST_INTEGRATION;}
};

#endif
