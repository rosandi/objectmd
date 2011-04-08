//-------------------------ForceDamper-------------------------//

#ifndef _FORCE_DAMPER_HPP_
#define _FORCE_DAMPER_HPP_

#include <omd/conditioner.hpp>

#define FREE_SURFACE -1.0

/**
 * @ingroup conditioner
 * @brief Implements the force dumping functionality
 * 
 * This class damps the force in the system with dumping coefficient defined as
 * "DFactor", the 2nd constructor's parameter. The function
 * works on the chosen target system (AtomContainer). 
 * This conditioner is the ForceModifier class, which operates
 * right after the forces are calculated.
 *
*/

class ForceDamper:public ForceConditioner {
	OMD_FLOAT Factor;
	
public:		

	ForceDamper(string target) {
		set_name("DAMPING");
		register_class(get_name());
		TargetName=target;
	}

	void ReadParameter() {
		string par=lower_case(replace_char(get_name(),' ','_'));
		Factor=SysParam->double_value(par+".factor");
	}

	void ForceModifier() {
		OMD_FLOAT na=Target->GetNAtom();
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

#endif
