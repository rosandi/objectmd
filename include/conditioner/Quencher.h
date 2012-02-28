//---------------------Quencher---------------------------------//

#ifndef _QUENCHER_HPP_
#define _QUENCHER_HPP_

#include <iostream>
#include <omd/conditioner.h>
using namespace omd;


/**
 * @ingroup conditioner
 * @brief Quench the system
 *
 * This class quench the velocity of all atoms in the system by
 * a quenching factor. The quenching factor is defined in the 
 * first parameter to the constructor.
 * The secong parameter is used to define the quencing desired
 * quencing period. The default values are 0.98 for quencing 
 * factor and 0 for period, which means that the system will be
 * quenched every loop.
 *
 * Parameters:
 * - @e quench.factor (factor) : quenching factor between 0 to 1
 * - @e quench.every (step) : the period of quenching in step
*/

class Quencher:public PostConditioner {
protected:
	int Period;
	double Factor;
	bool minimize;
	double last_pot;
	
public:
	Quencher() {
		set_name("QUENCH");
		register_class(get_name());
		Factor=-1.0;
		Period=1;
		minimize=false;
	}

	Quencher(string target, double quench_factor, int period=1) {
		set_name("QUENCH");
		register_class(get_name());
		Factor=quench_factor;
		Period=period;
		TargetName=target;
		minimize=false;
	}
	
	void ReadParameter() {
		// so... many quencher can be used...
		string tag=lower_case(replace_char(get_name(),' ','_'));
		SysParam->peek(tag+".factor", Factor);
		SysParam->peek(tag+".every", Period);
		SysParam->peek(tag+".minimize", minimize);
		
		mdassert(Factor>=0.0, 
			   "the quench factor ("+get_name()+
			   ".factor) must be defined in the parameter");

	}

	void PostIntegration() {
		if (!(NCalls%Period)) {
			int na=Target->GetNAtom();
			if(minimize) {
				if(last_pot<System->Potential) {
					for (int i=0; i<na; i++) {
						Atom* A=AtomPtr(i);
						A->vx*=Factor;
						A->vy*=Factor;
						A->vz*=Factor;
					}
				}
				last_pot=System->Potential;
			} else {
				for (int i=0; i<na; i++) {
					Atom* A=AtomPtr(i);
					A->vx*=Factor;
					A->vy*=Factor;
					A->vz*=Factor;
				}				
		 	}
	 	}
	}

	void PrintInfo(ostream& ost) {
		ost << "id." << id << " " << get_name() << " -- target="<<Target->get_name()
		<< "; factor="<< Factor << "; every=" << Period << " steps\n";
	}

};

#endif
