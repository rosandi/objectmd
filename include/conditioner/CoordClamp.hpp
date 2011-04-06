
#include <omd/conditioner.hpp>

/**
 * @ingroup conditioner
 * @brief Coordinate clamper
 *
 * Coordinate clamper clamps atoms in the chosen coordinate axes.
 * The force and velocity component of the corresponding axes is set to zero, 
 * The axes is determined my passing "coord_plane" parameter: "x", "xy", or "xyz", etc.
 * 
 * The class claims one bit from atom flag to mark the clamped status of an atom.
 * 
*/

#define CLAMP_X 1
#define CLAMP_Y 2
#define CLAMP_Z 4

class CoordClamp: public Force_Conditioner {
	bool use_range;
	int axis;
	int gid;	
	OMD_FLOAT West, East, South, North, Bottom, Top;
	int NClamp;
	int ClampBit;
	
	public:

		CoordClamp(string clamped_name, string clamped_axis)
		{
			TargetName=clamped_name;
			axis=clamped_axis;
			set_name("coordinate clamper");
			use_range=false;
			gid=-1;
		}
		
		CoordClamp(string clamped_name, int clamped_gid, int clamped_axis) {
			TargetName=clamped_name;
			gid=clamped_gid;
			axis=clamped_axis;
			use_range=false;
			set_name("coordinate clamper");
		}

		CoordClamp* SetRange(OMD_FLOAT west, OMD_FLOAT east, 
                             OMD_FLOAT south, OMD_FLOAT north, 
                             OMD_FLOAT bottom, OMD_FLOAT top) {
                             
			West=west;
			East=east;
			South=south;
			North=north;
			Bottom=bottom;
			Top=top;
			use_range=true;
			return this;
		}
		
		void UpdateRange() {
			if(!use_range)return;
			
			NClamp=0;
			int na=GetNAtom();
			for(int i=0; i<na; i++) {
				if( (Atoms(i).x>West   && Atoms(i).x<East)  &&
					(Atoms(i).y>South  && Atoms(i).y<North) &&
					(Atoms(i).z>Bottom && Atoms(i).z<Top))
				{
					SetFlag(i,ClampBit);
					NClamp++;
				} else {
					UnsetFlag(i,ClampBit);
				}
			}
		}
		
		void Init(MDSystem* WorkSys){
			Force_Conditioner::Init(WorkSys);
			ClampBit=ClaimFlagBit();
			if(gid>=0) {
				NClamp=0;
				for(int i=0;i<GetNAtom();i++){
					if(Atoms(i).gid==gid){
						Atoms(i).flag|=ClampBit;
						NClamp++;
					}
				}
				assert(NClamp>0, "can not find atoms to clamp. gid="+as_string(gid));
			} else {
				assert(use_range, "missing: coordinate range to clamp");
				UpdateRange();
			}
		}
		
		void ForceModifier() {
			int na=Target->GetNAtom();
			for (int i=0;i<na;i++) {
				if (!CheckFlag(i,ClampBit)) continue;
				if (axis|CLAMP_X) Target->Atoms(i).fx=Target->Atoms(i).vx=0.0;
				if (axis|CLAMP_Y) Target->Atoms(i).fy=Target->Atoms(i).vy=0.0;
				if (axis|CLAMP_Z) Target->Atoms(i).fz=Target->Atoms(i).vz=0.0;
			}
		}
		
		void PrintInfo(ostream& ost) {
			Force_Conditioner::PrintInfo(ost);
			ost<<"clamping "<<NClamp<<" atoms from "<<Target->GetName()<<"\n";
			ost<<"clamping axis : "
			   <<(axis|CLAMP_X?'X':' ')
			   <<(axis|CLAMP_Y?'Y':' ')
			   <<(axis|CLAMP_Z?'Z':' ')<<"\n";
		}
};
