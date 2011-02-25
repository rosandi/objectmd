//-------------------------Range-Limiter--------------------------------//

#ifndef _RANGE_LIMITER_HPP_
#define _RANGE_LIMITER_HPP_

#include <iostream>
#include <omd/conditioner.hpp>
#include <omd/container.hpp>

/**
 * ingroup conditioner
 * @brief Range limiter
 *
 * The class to limit the movement range of atoms inside the AtomContainer, pointed by
 * LimitedAtoms parameter. The range measured respects to the corresponding container 
 * borders. The atoms which move over the limit will be deactivated,
 * and excluded from the calculation. Be aware that atom deactivation disturbs the
 * energy conservation.
 * 
*/

class RangeLimiter: public Post_Conditioner {
	OMD_FLOAT StopBorder, lW, lE, lN, lS, lT, lB;	
public:
	
	RangeLimiter(string LimitedAtoms, OMD_FLOAT Border) 
	{
		StopBorder=Border;
		TargetName=LimitedAtoms;
		set_name("Range limiter");
	}
	
	void Init(MDSystem* WorkSys) {
		Conditioner::Init(WorkSys);
		SysBox Box=WorkSys->GetBox();
		lW = Box.x0-StopBorder;
		lE = Box.x1+StopBorder;
		lS = Box.y0-StopBorder;
		lN = Box.y1+StopBorder;
		lB = Box.z0-StopBorder;
		lT = Box.z1+StopBorder;
	}

	void PrintInfo(ostream& ost){
		ost << "ID." << id << " " << get_name()
			<< "\n Target: " << Target->get_name() 
			<< "\n Limits: W=" << lW << " E=" << lE << " S=" << lS 
			<< " N=" << lN << " B=" << lB << " T=" << lT << "\n";
	}
			 
	void PostIntegration() {
		for(OMD_INT i=0; i<GetNAtom(); i++) {			
			if(Atoms(i).x<lW){SetActive(i,false);continue;}
			if(Atoms(i).x>lE){SetActive(i,false);continue;}
			if(Atoms(i).y<lS){SetActive(i,false);continue;}
			if(Atoms(i).y>lN){SetActive(i,false);continue;}
			if(Atoms(i).z<lB){SetActive(i,false);continue;}
			if(Atoms(i).z>lT){SetActive(i,false);continue;}
		}
	}
};

#endif
