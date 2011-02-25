//---------------------Quencher---------------------------------//

#ifndef _QUENCHER_HPP_
#define _QUENCHER_HPP_

#include <iostream>
#include <omd/conditioner.hpp>

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
*/

class Quencher:public Post_Conditioner {
protected:
	OMD_INT Period, Step;
	OMD_FLOAT Factor;
	
public:
	Quencher(OMD_FLOAT QFactor, OMD_INT per=0) 
	{
		Factor=QFactor; Period=per; Step=0;
		set_name("Quencher");
	}
	
	void PrintInfo(ostream& ost) {
		ost << "ID." << id << " " << Name << " -- quench Factor=" 
		    << Factor << ", period=" << Period << "steps\n";
	}
		
	void PostIntegration() 
	{
		if (Step>=Period) {
			Step=0;
			OMD_INT na=Target->GetNAtom();
			for (OMD_INT i=0; i<na; i++) {
				Target->Atoms(i).vx*=Factor;
				Target->Atoms(i).vy*=Factor;
				Target->Atoms(i).vz*=Factor;
		 	}
	 	} 
	 	Step++;
	}
};

#endif
