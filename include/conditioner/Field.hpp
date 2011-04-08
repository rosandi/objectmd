
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
 * Not to be confused with CHARMS force field!
*/

class TField: public ForceConditioner {
	OMD_FLOAT ox,oy,oz;
	OMD_FLOAT ax,ay,az;
	OMD_FLOAT bx,by,bz;
	int    ivar;
	
	TableReader Phi;
	ForceFieldType Type;
	OMD_FLOAT CutRadius;
	
	public:		
		
		void CreateMe(const char* TargetAtoms, const char*tablefile) {
			TargetName=TargetAtoms;
			set_name("Force field");
			Phi.Open(tablefile);
			CutRadius=Phi.MaxRange();
		}

		TField(const char* TargetAtoms, const char* tablefile,
		           OMD_FLOAT originx, OMD_FLOAT originy, OMD_FLOAT originz) {
			ox=originx;oy=originy;oz=originz;
			Type=PointSource;
			CreateMe(TargetAtoms, tablefile);
		}			

		TField(const char* TargetAtoms, const char* tablefile, 
		           OMD_FLOAT originx, OMD_FLOAT originy, OMD_FLOAT originz,
		           OMD_FLOAT ariginx, OMD_FLOAT ariginy, OMD_FLOAT ariginz
		           ) {
			ox=originx;oy=originy;oz=originz;
			ax=ariginx;ay=ariginy;az=ariginz;
			Type=LineSource;
			CreateMe(TargetAtoms, tablefile);
		}
		
		TField(const char* iTargetAtoms, const char* tablefile, 
		           OMD_FLOAT originx, OMD_FLOAT originy, OMD_FLOAT originz,
		           OMD_FLOAT ariginx, OMD_FLOAT ariginy, OMD_FLOAT ariginz,
		           OMD_FLOAT briginx, OMD_FLOAT briginy, OMD_FLOAT briginz
		           ) {		
			ox=originx;oy=originy;oz=originz;
			ax=ariginx;ay=ariginy;az=ariginz;
			bx=briginx;by=briginy;bz=briginz;
			Type=PlaneSource;
			CreateMe(TargetAtoms, tablefile);
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
			OMD_FLOAT r, potential,drpot,ff,dx,dy,dz;
			for(int i=0;i<Target->NAtom;i++) {
				dx=Atoms[i].x;dy=Atoms[i].y;dz=Atoms[i].z;
				r=Distance(dx,dy,dz);
				if (r<CutRadius) {
					Phi.Read(r, potential, drpot);
					ff=-drpot/r;
					Atoms[i].fx+=dx*ff;
					Atoms[i].fy+=dy*ff;
					Atoms[i].fz+=dz*ff;
					POTENTIAL+=potential;
				}
			}
		}
};

#endif
