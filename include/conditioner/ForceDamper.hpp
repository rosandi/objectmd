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

class ForceDamper:public Force_Conditioner {
protected:
	OMD_FLOAT Factor;
	OMD_FLOAT dWest, dEast, dNorth, dSouth, dTop, dBottom;
	SysBox Box;
	
public:		

	ForceDamper(string TargetAtom, 
	            OMD_FLOAT DFactor, 
	            OMD_FLOAT dW, OMD_FLOAT dE, 
	            OMD_FLOAT dS, OMD_FLOAT dN, 
	            OMD_FLOAT dB, OMD_FLOAT dT);
	        
	void PrintInfo(ostream&);
	void Init(MDSystem* WorkSys);
	void ForceModifier();
	virtual void Modify(int idx);
	
	
};

//---------------------------ForceDamper---------------------------//

ForceDamper::ForceDamper(string TargetAtom, 
	                     OMD_FLOAT DFactor, 
	                     OMD_FLOAT dW, OMD_FLOAT dE, 
	                     OMD_FLOAT dS, OMD_FLOAT dN, 
	                     OMD_FLOAT dB, OMD_FLOAT dT) 
{
	Factor=DFactor;
	if(TargetAtom!=""||TargetAtom!="all")TargetName=TargetAtom; 
	set_name("Force damper");
	dWest=dW;dEast=dE; 
	dSouth=dS;dNorth=dN; 
	dBottom=dB;dTop=dT;  
}
	
void ForceDamper::Init(MDSystem* WorkSys)
{	
	Conditioner::Init(WorkSys);
	Box=System->GetBox();
	dWest  =Box.x0+dWest;
    dEast  =Box.x1-dEast;
    dSouth =Box.y0+dSouth;
    dNorth =Box.y1-dNorth;
    dBottom=Box.z0+dBottom;
    dTop   =Box.z1-dTop;
}

// Factor=0 means no damper
// in "impact-code" there is 0.5*Factor is used
// here, put all numbers inside Factor!

void ForceDamper::Modify(int idx) {
	Atoms(idx).fx -= Factor*Atoms(idx).vx;
	Atoms(idx).fy -= Factor*Atoms(idx).vy;
	Atoms(idx).fz -= Factor*Atoms(idx).vz;
}    

void ForceDamper::ForceModifier()
{
	OMD_FLOAT na=Target->GetNAtom();
 	for(int i=0;i<na;i++) {
 		
 		// check also free surfaces.. beyond free surface -> free!
 		
 		if(Atoms(i).x<dWest) {
 			if(Box.x0<dWest) Modify(i); // not a free surface
 			continue;
		}
		
		if(Atoms(i).x>dEast) {
			if(Box.x1>dEast) Modify(i);
			continue;
		}
		
		if(Atoms(i).y<dSouth) {
			if(Box.y0<dSouth) Modify(i);
			continue;
		}
		
		if(Atoms(i).y>dNorth) {
			if(Box.y1>dNorth) Modify(i);
			continue;
		}
		
 		if(Atoms(i).z<dBottom) {
 			if(Box.z0<dBottom) Modify(i);
 			continue;
		}
		
		if(Atoms(i).z>dTop) {
			if(Box.z1>dTop) Modify(i);
			continue;
		}
   	}
}

void ForceDamper::PrintInfo(ostream& ost) {
		ost << "ID." << id << " " << get_name() << "; Borders:";
 	    if (Box.x0<dWest)   ost << " W=" << dWest; else ost << " W=free";
 	    if (Box.x1>dEast)   ost << " E=" << dEast; else ost << " E=free";
 	    if (Box.y0<dSouth)  ost << " S=" << dSouth; else ost << " S=free";
 	    if (Box.y1>dNorth)  ost << " N=" << dNorth; else ost << " N=free";
 	    if (Box.z0<dBottom) ost << " B=" << dBottom; else ost << " B=free";
 	    if (Box.z1>dTop)    ost << " T=" << dTop; else ost << " T=free";
 	    ost << " Factor=" << Factor << '\n';
}

#endif
