//---------------------Quencher---------------------------------//

#ifndef _POTENTIAL_MINIMIZER_HPP_
#define _POTENTIAL_MINIMIZER_HPP_

#include <iostream>
#include <omd/conditioner.h>
using namespace omd;


/**
 * @ingroup conditioner
 * @brief Minimize the potential energy of the system.
 * 
 * The atom velocity is quenched if the potential is increased.
 * 
*/

class PotentialMinimizer:public PostConditioner {
protected:
	double factor;
	double last_pot;
	
public:
	PotentialMinimizer(double QFactor) 
	{
		factor=QFactor;
		last_pot=DBL_MAX;
		set_name("potential minimizer");
	}
	
	void PrintInfo(ostream& ost) {
		ost << "ID." << id << " " << get_name()
		    << " -- quench Factor=" << factor << "\n";
	}
		
	void PostIntegration() 
	{
		if (last_pot<System->Potential) {
			int na=Target->GetNAtom();			
			for (int i=0; i<na; i++) {
				Target->Atoms(i).vx*=factor;
				Target->Atoms(i).vy*=factor;
				Target->Atoms(i).vz*=factor;
		 	}
	 	}
	 	last_pot=System->Potential;
	}
};

#endif
