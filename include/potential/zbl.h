/* lib
 ***************************************************
 * Object-MD
 * Object Oriented Molecular Dynamics Program
 * 
 * (c)Yudi Rosandi, 2005
 *
 * Version 1.0 (07.10.07)
 *
 * Project started on July 2005.
 *
 ***************************************************
 *
 *
*/
 
#include <omd/forcekernel.h>

NOT FINISHED!

class zbl: public ForceKernel {	
	
	double zc,ze0,ze1,ze2,ze3;

	public:
	
	zbl(double rcut) {
		set_name("ZBL");
		register_class(get_name());
		CutRadius=rcut;
		CutRadiusSqr=CutRadius*CutRadius;
	}

	void Init(MDSystem* WorkSys) {
		ForceKernel::Init(WorkSys);
		double znum1=System->SystemAtoms[AtomTypeA]->Z;
		double znum2=System->SystemAtoms[AtomTypeB]->Z;
		double znor=0.468521928/(pow(znum1,0.23)+pow(znum2, 0.23));
		zc =14.39964415*znum1*znum2;
		ze0=3.2000/znor;
		ze1=0.9423/znor;
		ze2=0.4029/znor;
		ze3=0.2016/znor;		
	}

	void ComputeHalf(Atom& at, Atom& to) {
		
		double fex[]= {
			0.1818 *exp(-ze0*r),
			0.5099 *exp(-ze1*r),
			0.2802 *exp(-ze2*r),
			0.02817*exp(-ze3*r)
		};
		
		double fpot=(fex[0]+fex[1]+fex[2]+fex[3])*zc/r;
		double fr=((1.0/r+ze0)*fex[0]+
					  (1.0/r+ze1)*fex[1]+
					  (1.0/r+ze2)*fex[2]+
					  (1.0/r+ze3)*fex[3])*zc/(r*r);
		
		ReturnForce(at,to,dx,dy,dz,fr,fpot);

	}

};
