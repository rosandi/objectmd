
#ifndef _T_FORCE_FIELD_HPP_
#define _T_FORCE_FIELD_HPP_

#include <omd/conditioner.hpp>
#include <omd/tablereader.hpp>

enum ForceFieldType {PointSource, LineSource, PlaneSource};

/**
 * @ingroup conditioner
 * @brief Interaction with point, line, or surface potential source
 *
 * This class read potential function from a table.
 * The type of the forces is defined by the constructors.
 * The calculated potential is added to system's potential.
 *
 *
 * parameters:
 *  - (name).point ox oy oz : position of a point field
 *  - (name).line ox oy oz ax ay az : position of a line field
 *  - (name).plane ox oy oz ax ay az bx by bz : position of a plane field
 *  - (name).table tablename : potential table of the field
 *  - field.table tablename : the default table used if (name).table is not defined
 *
 * (name) represent the name of the gadget (see SetName()).
*/

class Field: public ForceConditioner {
	OMD_FLOAT ox,oy,oz;
	OMD_FLOAT ax,ay,az;
	OMD_FLOAT bx,by,bz;
	int    ivar;
	
	TableReader Phi;
	ForceFieldType Type;
	OMD_FLOAT CutRadius;
	
public:		
	Field(string TargetAtoms) {
		TargetName=TargetAtoms;
	}			

	void ReadParameter() {
		string tag=replace_char(lower_case(get_name()),' ', '_');
		
		if(SysParam->exist(tag+".table")) {
			Phi.open(SysParam->string_value(tag+".table"));
		} else {
			// default: use field.table
			Phi.open(SysParam->string_value("field.table"));
		}

		CutRadius=Phi.max_range();
				
		if(SysParam->exist(tag+".point")) {
			Type=PointSource;
			tag.append(".point");
			ox=SysParam->double_value(tag,0);
			oy=SysParam->double_value(tag,1);
			oz=SysParam->double_value(tag,2);			
		} else if(SysParam->exist(tag+".line")) {
			Type=LineSource;
			tag.append(".line");
			ox=SysParam->double_value(tag,0);
			oy=SysParam->double_value(tag,1);
			oz=SysParam->double_value(tag,2);			
			ax=SysParam->double_value(tag,3);
			ay=SysParam->double_value(tag,4);
			az=SysParam->double_value(tag,5);
		} else if(SysParam->exist(tag+".plane")) {
			Type=PlaneSource;
			tag.append(".plane");
			ox=SysParam->double_value(tag,0);
			oy=SysParam->double_value(tag,1);
			oz=SysParam->double_value(tag,2);			
			ax=SysParam->double_value(tag,3);
			ay=SysParam->double_value(tag,4);
			az=SysParam->double_value(tag,5);
			bx=SysParam->double_value(tag,6);
			by=SysParam->double_value(tag,7);
			bz=SysParam->double_value(tag,8);
		} else 
			die("parameter needed: "+tag+".[point|line|plane]");
	}
		
	OMD_FLOAT Distance(OMD_FLOAT& x, OMD_FLOAT& y, OMD_FLOAT& z) {
		OMD_FLOAT norm;
		if(Type==PointSource) {
			x=x-ox; y=y-oy; z=z-oz;
			return sqrt(x*x+y*y+z*z);
		}
		if(Type==LineSource) {
			bx=x-ox;by=y-oy;bz=z-oz;
			OMD_FLOAT cx=ax-ox,cy=ay-oy,cz=az-oz;
			OMD_FLOAT norm=sqrt(cx*cx+cy*cy+cz*cz);
			OMD_FLOAT ex=cx/norm, ey=cy/norm, ez=cz/norm;
			OMD_FLOAT ld=bx*ex+by*ey+bz*ez;
			OMD_FLOAT dx=ld*ex, dy=ld*ey,dz=ld*ez;
			x=bx-dx;y=by-dy;z=bz-dz;
			return sqrt(x*x+y*y+z*z);
		}
		if(Type==PlaneSource) {
			OMD_FLOAT cx=ax-ox,cy=ay-oy,cz=az-oz;
			OMD_FLOAT dx=bx-ox,dy=by-oy,dz=bz-oz;
			OMD_FLOAT ex=(cy*dz-cz*dy);
			OMD_FLOAT ey=(cz*dx-cx*dz);
			OMD_FLOAT ez=(cx*dy-cy*dx);
			norm=sqrt(ex*ex+ey*ey+ez*ez);
			x=(x-ox)*ex/norm;y=(y-oy)*ey/norm;z=(z-oz)*ez/norm;
			return sqrt(x*x+y*y+z*z);
		}
	}
			
	void ForceModifier(){
		OMD_FLOAT r, pot,drpot,ff,dx,dy,dz;
		int na=Target->GetNAtom();
		for(int i=0;i<na;i++) {
			dx=Atoms(i).x;dy=Atoms(i).y;dz=Atoms(i).z;
			r=Distance(dx,dy,dz);
			if(r==0.0) continue; // just in case...
			if (r<CutRadius) {
				Phi.read(r, pot, drpot);
				ff=-drpot/r;
				Atoms(i).fx+=dx*ff;
				Atoms(i).fy+=dy*ff;
				Atoms(i).fz+=dz*ff;
				Atoms(i).potential+=pot;
			}
		}		
	}
};

#endif
