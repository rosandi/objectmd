
#ifndef _CRYSTAL_BCC_100_HPP_
#define _CRYSTAL_BCC_100_HPP_

#include <crystal/FCC111.hpp>


/**
  @ingroup atom
  @brief BCC(100) crystal.
**/

class CrystalBCC100: public CrystalFCC111 {
public:
	CrystalBCC100(OMD_INT XMLayer, OMD_INT YMLayer, OMD_INT ZMLayer, string mat_file):
	CrystalFCC111(XMLayer, YMLayer, ZMLayer, mat_file){}
	
	AtomContainer* Create() {
		OMD_INT eNx=(xml+1)/2;
		OMD_INT eNy=(yml+1)/2;
		OMD_INT oNx=xml/2;
		OMD_INT oNy=yml/2;
		XMLDist=YMLDist=ZMLDist=lattice_constant/2.0;

		OMD_INT na=zml*(eNx*eNy+oNx*oNy)/2 + (zml%2)*(eNx*eNy);
		Allocate(na);

		OMD_INT n=0;
		for (OMD_INT zm=0; zm<zml; zm++) {		
			OMD_INT nx=(zm%2)?oNx:eNx;
			OMD_INT ny=(zm%2)?oNy:eNy;
			OMD_FLOAT ofs=(zm%2)?lattice_constant/2.0:0.0;	
			for (OMD_INT a=0; a<nx; a++)
				for (OMD_INT b=0; b<ny; b++) {
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
