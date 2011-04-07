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
		
		// create minimum 1 unit cell 4.xml, 4.yml, 4.zml
		// crystal size will be rounded to in unit cell size
		
		int ncx=xml/4+1;
		int ncy=yml/4+1;
		int ncz=zml/4+1;
		Allocate(ncx*ncy*ncz*8);
		
		blog("creating cell: "+
			 as_string(ncx)+"x"+as_string(ncy)+"x"+as_string(ncz)+
			 " size="+as_string(GetNAtom())
			 );
		
		double a=lattice_constant;
		double bx=0.0;
		Atom* A=AtomPtr(0);
		
		for(int i=0;i<ncx;i++) {
			double by=0.0;
			for(int j=0;j<ncy;j++) {
				double bz=0.0;
				for(int k=0;k<ncz;k++) {
					A->x=bx;         A->y=by;         A->z=bz;         A++;
					A->x=bx+a/2.;    A->y=by+a/2.;    A->z=bz;         A++;
					A->x=bx+a/4.;    A->y=by+a/4.;    A->z=bz-a/4.;    A++;
					A->x=bx+3.*a/4.; A->y=by+3.*a/4.; A->z=bz-a/4.;    A++;
					A->x=bx;         A->y=by+a/2.;    A->z=bz-a/2.;    A++;
					A->x=bx+a/2.;    A->y=by;         A->z=bz-a/2.;    A++;
					A->x=bx+a/4.;    A->y=by+3.*a/4.; A->z=bz-3.*a/4.; A++;
					A->x=bx+3.*a/4.; A->y=by+a/4.;    A->z=bz-3.*a/4.; A++;
					bz-=a;
				}
				by+=a;
			}
			bx+=a;
		}
		
		XMLDist=YMLDist=ZMLDist=a/4.0;

	}
	
	void Create100110() {
		// minimum unit cell 2.xml, 2.yml, 4.zml
		int ncx=xml/2+1;
		int ncy=yml/2+1;
		int ncz=zml/4+1;
		
		Allocate(ncx*ncy*ncz*4);
		
		blog("creating cell: "+
			 as_string(ncx)+"x"+as_string(ncy)+"x"+as_string(ncz)+
			 " size="+as_string(GetNAtom())
			 );		
		
		double a=lattice_constant;
		double dx=0.5*sqrt(2.0)*a;
		double dy=0.5*sqrt(2.0)*a;
		double dz=a;
		
		Atom* A=AtomPtr(0);
		int idx=0;
		double bx=0.0;
		for(int i=0;i<ncx;i++) {
			double by=0.0;
			for(int j=0;j<ncy;j++) {
				double bz=0.0;
				for(int k=0;k<ncz;k++) {

					A->x=bx;        A->y=by;        A->z=bz;            A++;
					A->x=bx+dx/2.0; A->y=by;        A->z=bz-dz/4.0;     A++;
					A->x=bx+dx/2.0; A->y=by+dy/2.0; A->z=bz-dz/2.0;     A++;
					A->x=bx;        A->y=by+dy/2.0; A->z=bz-3.0*dz/4.0; A++;
					
					bz-=dz;
				}
				by+=dy;
			}
			bx+=dx;
		}
		
		XMLDist=dx/2.0;
		YMLDist=dy/2.0;
		ZMLDist=dz/8.0;

	}
	
	virtual AtomContainer* Create() {
		
		if(created) return this;
		
		if(orientation=="100") { // diamond basic unit cell
			Create100();
			created=true;
			CalcBox();
			return this;
		}
		
		if(orientation=="100_110") { // z-surface=(100) x-axis=(110) diagonal
			Create100110();
			created=true;
			CalcBox();
			return this;
		}
		
		die("orientation:"+orientation+" is not implemented...");

	}
};

#endif
