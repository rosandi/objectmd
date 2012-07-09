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
 *
*/

#ifndef _MODIFY_H_
#define _MODIFY_H_

#include <omd/gadget.h>

namespace omd {

#define MODIFY_PRE_INTEGRATION  1
#define MODIFY_PRE_CALCULATION  2
#define MODIFY_POST_FORCE   4
#define MODIFY_POST_INTEGRATION 8

//------------------Modify-----------------//
/**
 * @ingroup modify
 * @brief Class to modify system configurations and atoms
 *
 * The modify class is used to make modification on system while the
 * simulation is running. The following are the member functions that can be 
 * implemented in the child class,
 * 	 - PreIntegration(), called before system integration. 
 *   - PreCalculation(), called by the integrator class before the force 
 *     calculation loop.
 *   - PostForce(), called by the integrator after the force 
 *     calculation loop.
 *   - PostIntegration(), called after the integration, before 
 *     executing Detectors.
 * 
 * The Modify::Active variable identifies the type of modify:
 * 
 *	- 1 --> Pre-integration
 *	- 2 --> Pre-calculation
 *	- 4 --> Post-Force
 *	- 8 --> Post-integration
 * 
 * Zero value is always mean inactive. The bits of this variable is be examined
 * in the Modify::Execute() function. The calling number (NCalls) will be increased
 * when the corresponding integrator function is executed. If a modify has 
 * multiple-type, NCalls will be increased in every modify functions.
 * (NCalls/step=number_of_type)
 *
*/

class Modify: public MDGadget {

public:
	Modify() {Active=15; Target=NULL; TargetName="";}
	
	virtual ~Modify() {blog("cleaning up", LOGDESTROY);}
	
	void PrintInfo(ostream& ost) {
		ost<<"id."<<id<<" "<<get_name()<<" -- target: "<<Target->get_name()<<std::endl;
	}
	
	void SetModifyType(int condtype) {Active=ActiveCode=condtype;} 

    virtual void PreIntegration()  {} // before integration
    virtual void PreCalculation()  {} // before force calculation loop
    virtual void PostForce()   {} // after force calculation
    virtual void PostIntegration() {} // after integration, before detection

    virtual void Execute(int condtype) {
		if(!Active) return; // Active stors the type of modify...
    	if(condtype&MODIFY_PRE_INTEGRATION&Active){PreIntegration(); NCalls++;}
    	if(condtype&MODIFY_PRE_CALCULATION&Active){PreCalculation(); NCalls++;}
    	if(condtype&MODIFY_POST_FORCE&Active){PostForce(); NCalls++;}
    	if(condtype&MODIFY_POST_INTEGRATION&Active){PostIntegration(); NCalls++;}
	}
	
	virtual Modify* SetName(string newname) {
		set_name(newname);
		return this;
	}

};

/** @brief Pre integration modify **/
class PreModify:public Modify {
	public:
		PreModify(){Active=ActiveCode=MODIFY_PRE_INTEGRATION;}
};

/** @brief Pre force calculation loop modify **/
class CalcModify:public Modify {
	public:
		CalcModify(){Active=ActiveCode=MODIFY_PRE_CALCULATION;}
};

/** @brief Force modifier container **/
class ForceModify:public Modify {
	public:
		ForceModify(){Active=ActiveCode=MODIFY_POST_FORCE;}
};

/** @brief Post integration modify **/
class PostModify:public Modify {
	public:
		PostModify(){Active=ActiveCode=MODIFY_POST_INTEGRATION;}
};

}

#endif
