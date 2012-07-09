//-------------------------ForceDamper-------------------------//

#ifndef _FORCE_DAMPER_HPP_
#define _FORCE_DAMPER_HPP_

#include <omd/modify.h>

namespace omd {


#define FREE_SURFACE -1.0

/**
 * @ingroup modify
 * @brief Implements the force dumping functionality
 * 
 * This class damps the force in the system with dumping coefficient defined as
 * "DFactor", the 2nd constructor's parameter. The function
 * works on the chosen target system (AtomContainer). 
 * This modify is of PostForce type, operates after the forces calculation.
 *
*/

class ForceDamper:public ForceModify {
	double Factor;
	
public:		

	ForceDamper() {
		set_name("damp");
		register_class(get_name());
	}

	void ReadParameter() {
		Factor=SysParam->double_value(mytag("factor"));
    TargetName=SysParam->string_value(mytag("target"));
	}

	void PostForce() {
		double na=GetNAtom();
		for(int i=0;i<na;i++) {
			Atoms(i).fx-=Factor*Atoms(i).vx;
			Atoms(i).fy-=Factor*Atoms(i).vy;
			Atoms(i).fz-=Factor*Atoms(i).vz;		
		}
	}

	void PrintInfo(ostream& ost) {
		ost << "id." << id << " " << get_name()
			<<" -- damping ("<<Target->get_name()<<"); factor=" << Factor << "\n";
	}
};

}

#endif
