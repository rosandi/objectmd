//-------------------------Range-Limiter--------------------------------//

#ifndef _RANGE_LIMITER_HPP_
#define _RANGE_LIMITER_HPP_

#include <iostream>
#include <omd/conditioner.h>
#include <omd/container.h>
using namespace omd;


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

class RangeLimiter: public PostConditioner {
	double StopBorder, lW, lE, lN, lS, lT, lB;	
public:
	
	RangeLimiter(string LimitedAtoms, double Border) 
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
		int na=GetNAtom();
		for(int i=0; i<na; i++) {
			Atom* a=AtomPtr(i);
			if(a->x<lW){SetActive(i,false);continue;}
			if(a->x>lE){SetActive(i,false);continue;}
			if(a->y<lS){SetActive(i,false);continue;}
			if(a->y>lN){SetActive(i,false);continue;}
			if(a->z<lB){SetActive(i,false);continue;}
			if(a->z>lT){SetActive(i,false);continue;}
		}
	}
};

#endif
