
#ifndef _CRYSTAL_BCC_100_HPP_
#define _CRYSTAL_BCC_100_HPP_

#include <crystal/FCC111.hpp>


/**
  @ingroup atom
  @brief BCC(100) crystal.
**/

class CrystalBCC100: public CrystalFCC111 {
public:
	CrystalBCC100(int XMLayer, int YMLayer, int ZMLayer, string mat_file):
	CrystalFCC111(XMLayer, YMLayer, ZMLayer, mat_file){}
	
	AtomContainer* Create() {
		int eNx=(xml+1)/2;
		int eNy=(yml+1)/2;
		int oNx=xml/2;
		int oNy=yml/2;
		XMLDist=YMLDist=ZMLDist=lattice_constant/2.0;

		int na=zml*(eNx*eNy+oNx*oNy)/2 + (zml%2)*(eNx*eNy);
		Allocate(na);

		int n=0;
		for (int zm=0; zm<zml; zm++) {		
			int nx=(zm%2)?oNx:eNx;
			int ny=(zm%2)?oNy:eNy;
			OMD_FLOAT ofs=(zm%2)?lattice_constant/2.0:0.0;	
			for (int a=0; a<nx; a++)
				for (int b=0; b<ny; b++) {
					Atoms(n).x= (OMD_FLOAT)a *lattice_constant+ofs;
					Atoms(n).y= (OMD_FLOAT)b *lattice_constant+ofs;
					Atoms(n).z=-(OMD_FLOAT)zm*lattice_constant/2.0;
					n++;
				}
		}
		CalcBox();
		created=true;
		return this;
	}
};

#endif
