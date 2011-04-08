
#include <omd/conditioner.hpp>

/**
 * @ingroup conditioner
 * @brief Atom Coordinate Clamper
 *
 * This class clamps atoms in the chosen coordinate axis.
 * The force and velocity component of the corresponding axes is set to zero, 
 * The axes is determined my passing "coord_plane" parameter: "x", "xy", or "xyz", etc.
 * The clamped atoms must be a member of the Target container (via AtomGroup).
 * 
 * 
*/

#define CLAMP_X 1
#define CLAMP_Y 2
#define CLAMP_Z 4

class CoordClamp: public ForceConditioner {
	string saxis;
	int axis;
	
public:

	CoordClamp(string target="-", string ax="xyz") {
		set_name("CLAMP");
		register_class(get_name());
		TargetName=target;
		saxis=ax;
		axis=0;
	}
	
	void ReadParameter() {
		string tag=lower_case(replace_char(get_name(),' ','_'));
		SysParam->peek(tag+".axis", saxis);
		SysParam->peek(tag+".target", TargetName);
	}
		
	void Init(MDSystem* WorkSys){
		ForceConditioner::Init(WorkSys);
		if(saxis.find('x')) axis|=CLAMP_X;
		if(saxis.find('y')) axis|=CLAMP_Y;
		if(saxis.find('z')) axis|=CLAMP_Z;
	}
	
	void ForceModifier() {
		int na=GetNAtom();
		for (int i=0;i<na;i++) {
			if (axis|CLAMP_X) Atoms(i).fx=Atoms(i).vx=0.0;
			if (axis|CLAMP_Y) Atoms(i).fy=Atoms(i).vy=0.0;
			if (axis|CLAMP_Z) Atoms(i).fz=Atoms(i).vz=0.0;
		}
	}
	
	void PrintInfo(ostream& ost) {
		ost << "id." << id << " " << get_name()
		    <<" -- clamping ("<<Target->get_name()<<"); axis="<<saxis<<"\n";
	}
};
