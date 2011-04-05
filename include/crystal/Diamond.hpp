//----------------------Diamond Structure----------------------//

#ifndef _CRYSTAL_DIAMOMD_HPP_
#define _CRYSTAL_DIAMOND_HPP_

#include <crystal/Crystalline.hpp>

class Diamond:public Crystalline {

public: 
	Diamond(string orientation, int xml, int yml, int zml, string mat_file):
	Crystalline(orientation, xml, yml, zml, mat_file){
		set_name("DIAMOND");
		register_class(get_name());
	}
	
	void Create100() {
		
		// create minimum 1 unit cell
		// crystal size will be rounded to in unit cell size
		
		int ncx=xml/4+1;
		int ncy=yml/4+1;
		int ncz=zml/4+1;
		
		Allocate(ncx*ncy*ncz*8);
		
		double a=lattice_constant;
		double bx=0.0;
		Atom* A=AtomPtr(0);
		
		for(int i=0;i<ncx;i++) {
			double by=0.0;
			for(int j=0;j<ncy;j++) {
				double bz=0.0;
				for(int k=0;k<ncz;k++) {
					A->x=bx;       A->y=by;       A->z=bz; A++;
					A->x=bx+a/2;   A->y=by+a/2;   A->z=bz; A++;
					A->x=bx+a/4;   A->y=by+a/4;   A->z=bz-a/4; A++;
					A->x=bx+3*a/4; A->y=by+3*a/4; A->z=bz-a/4; A++;
					A->x=bx;       A->y=by+a/2;   A->z=bz-a/2; A++;
					A->x=bx+a/2;   A->y=by;       A->z=bz-a/2; A++;
					A->x=bx+a/4;   A->y=by+3*a/4; A->z=bz-3*a/4; A++;
					A->x=bx+3*a/4; A->y=by+a/4;   A->z=bz-3*a/4; A++;
					bz-=a;
				}
				by+=a;
			}
			bx+=a;
		}
		
		XMLDist=YMLDist=ZMLDist=a/8.0;
		
		DumpAtoms("cry");
		
	}
	
	virtual AtomContainer* Create() {
		
		if(orientation=="100") {
			Create100();
			created=true;
			CalcBox();
			return this;
		}
		
		die("orientation:"+orientation+" is not implemented...");

	}
};

#endif
