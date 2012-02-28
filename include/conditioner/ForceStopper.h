//------------------------ForceStopper---------------------------//
#ifndef _FORCE_STOPPER_HPP_
#define _FORCE_STOPPER_HPP_

#include <conditioner/ForceDamper.h>

/**
 * @ingroup conditioner
 * @brief Stops the atom at given region
 *
 * This class it stops atoms in a defined region.
*/

class ForceStopper: public ForceConditioner {
	double dWest, dEast, dNorth, dSouth, dTop, dBottom;
	int fixflag;

public:

	ForceStopper(string TargetAtom, 
	            double dW, double dE, 
	            double dS, double dN, 
	            double dB, double dT)
	{
		set_name("FORCE STOPPER");
		if(TargetAtom!=""||TargetAtom!="all")TargetName=TargetAtom; 
		dWest=dW;dEast=dE; 
		dSouth=dS;dNorth=dN; 
		dBottom=dB;dTop=dT;
	}

	void Init(MDSystem* WorkSys){
		Conditioner::Init(WorkSys);
		SysBox Box=System->GetBox();
		dWest  =Box.x0+dWest;
		dEast  =Box.x1-dEast;
		dSouth =Box.y0+dSouth;
		dNorth =Box.y1-dNorth;
		dBottom=Box.z0+dBottom;
		dTop   =Box.z1-dTop;
		
		fixflag=ClaimFlagBit();
		
		for(int i=0; i<GetNAtom(); i++) {
 			if(Atoms(i).x<dWest){SetFlag(i,fixflag);continue;}
			if(Atoms(i).x>dEast){SetFlag(i,fixflag);continue;}
			if(Atoms(i).y<dSouth){SetFlag(i,fixflag);continue;}
			if(Atoms(i).y>dNorth){SetFlag(i,fixflag);continue;}
			if(Atoms(i).z<dBottom){SetFlag(i,fixflag);continue;}		
			if(Atoms(i).z>dTop){SetFlag(i,fixflag);continue;}			
		}
	}

	void ForceModifier(){
		double na=GetNAtom();
		for(int i=0;i<na;i++) {
			if(CheckFlag(i,fixflag)) {
				Atoms(i).fx = Atoms(i).vx = 0.0;
				Atoms(i).fy = Atoms(i).vy = 0.0;
				Atoms(i).fz = Atoms(i).vz = 0.0;
			}
		}
	}
         
	void PrintInfo(ostream& ost){
		SysBox* Box=&(System->GetBox());
		ost << "ID." << id << " " << get_name() << "; Borders:";
		if (Box->x0<dWest)   ost << " W=" << dWest; else ost << " W=free";
		if (Box->x1>dEast)   ost << " E=" << dEast; else ost << " E=free";
		if (Box->y0<dSouth)  ost << " S=" << dSouth; else ost << " S=free";
		if (Box->y1>dNorth)  ost << " N=" << dNorth; else ost << " N=free";
		if (Box->z0<dBottom) ost << " B=" << dBottom; else ost << " B=free";
		if (Box->z1>dTop)    ost << " T=" << dTop; else ost << " T=free";
		ost << '\n';
	}

};

#endif
