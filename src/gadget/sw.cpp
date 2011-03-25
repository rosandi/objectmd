/* omdlib
 ********************************************************
 *
 * ObjectMD
 * Molecular Dynamics Class Library
 * Version 3.1 (2009,2011)
 *
 * Released under:
 *       GNU GENERAL PUBLIC LICENSE Version 3
 *       Copyright (C) 2005  Yudi Rosandi
 *
 * Read license.txt file in the source's root directory.
 *
 ********************************************************
 * Stillinger-Weber Potential
 *
*/

#include <omd/system.hpp>
#include <conditioner/VerletList.hpp>
#include <potential/sw.hpp>

StillingerWeber::StillingerWeber(string material) {
	paramfile=material;
	set_name("STILLINGER WEBER");
	register_class(get_name());
}

void StillingerWeber::ReadParameter() {
	param.read(search_path("$OMD_TABLE", "sw."+paramfile));
	eps=param.double_value("epsilon");
	A=param.double_value("A");
	B=param.double_value("B");
	p=param.double_value("p");
	q=param.double_value("q");
	sigma=param.double_value("sigma");
	alpha=param.double_value("a");
	lambda=param.double_value("lambda");
	gamma=param.double_value("gamma");
	lc=param.double_value("lattice_constant");
}

void StillingerWeber::Init(MDSystem* WorkSys) {
	ForceKernel::Init(WorkSys);
	
	Verlet=dynamic_cast<VerletList*>(System->GetIterator());
	assert(Verlet->type_of("verlet list"), 
		   "StillingerWeber potential kernel needs (VERLET LIST) iterator");
	
	CutRadius=alpha*sigma;
	CutRadiusSqr=CutRadius*CutRadius;
    K=new OMD_FLOAT[10];
    K[0]=A*eps*p*B*pow(sigma,p);
    K[1]=A*eps*q*pow(sigma,q);
    K[2]=A*eps*B*pow(sigma,p+1);
    K[3]=A*eps*pow(sigma,q+1);
    K[4]=A*eps*B*pow(sigma,p);
    K[5]=A*eps*pow(sigma,q);
	
    sigga=sigma*gamma;//*2./lc;
    rsisi=sigma*alpha;//*2./lc;
    onethird=1.0/3.0;
    epsla=eps*lambda;
}

StillingerWeber::~StillingerWeber() {
	delete[] K;
}

void StillingerWeber::PrintInfo(ostream& ost) {
	ost << "id."<<id<<" "<<get_name()<<"\n"
		<< "atoms "<<System->SystemAtoms[AtomTypeA]->get_name()
		<< "<->"<< System->SystemAtoms[AtomTypeB]->get_name()<< "\n"
		<<" parameter="<<paramfile;
}

void StillingerWeber::TwoBodyTerm(Atom& at, Atom& to,
								  OMD_FLOAT r,
								  OMD_FLOAT dx,
								  OMD_FLOAT dy,
								  OMD_FLOAT dz)
{
	OMD_FLOAT rp=pow(r,-p);
	OMD_FLOAT rq=pow(r,-q);
	OMD_FLOAT irr=1./(r*r);
	OMD_FLOAT ira=(1.0/(r-CutRadius));
	OMD_FLOAT irra=ira*ira*r;
	OMD_FLOAT esira=exp(sigma*ira);
	OMD_FLOAT fr=(K[0]*rp-K[1]*rq+(K[2]*rp-K[3]*rq)*irra)*esira*irr;
	OMD_FLOAT ep=0.5*(K[4]*rp-K[5]*rq)*esira;

	// return force...
	at.fx+=dx*fr;
	at.fy+=dy*fr;
	at.fz+=dz*fr;
	to.fx-=dx*fr;
	to.fy-=dy*fr;
	to.fz-=dz*fr;
	at.potential+=ep;
	to.potential+=ep;
}

void StillingerWeber::ThreeBodyTerm(Atom& at0, Atom& at1, Atom& at2,
                                    OMD_FLOAT rjl, // at <- first_nb
                                    OMD_FLOAT dxjl,
                                    OMD_FLOAT dyjl,
                                    OMD_FLOAT dzjl,
                                    OMD_FLOAT ril, // at <- second_nb (iterated
                                    OMD_FLOAT dxil,
                                    OMD_FLOAT dyil,
                                    OMD_FLOAT dzil)
{
    double edil=1.0/ril;
    double edjl=1.0/rjl;
    double f1il=exp(sigga/(ril-rsisi));
    double f3il=sigga/((ril-rsisi)*(ril-rsisi)*ril);
    double f1jl=exp(sigga/(rjl-rsisi));
    double f3jl=sigga/((rjl-rsisi)*(rjl-rsisi)*rjl);
	
    double splilj=(dxil*dxjl+dyil*dyjl+dzil*dzjl)*edil*edjl;
    double f2lilj=splilj+onethird;
    double f4lilj=splilj*pow(edjl,2.0);
    double f4ljli=splilj*pow(edil,2.0);
    double f5lilj=edil*edjl;
    double ep=epsla*f1il*f1jl*f2lilj*f2lilj/3.0;
	
    double vorfak=epsla*f1il*f1jl*f2lilj;
    double fril=f2lilj*f3il-2.0*(f5lilj-f4ljli);
    double frjl=f2lilj*f3jl-2.0*(f5lilj-f4lilj);

	double frx=-vorfak*(fril*dxil+frjl*dxjl);
    double fry=-vorfak*(fril*dyil+frjl*dyjl);
    double frz=-vorfak*(fril*dzil+frjl*dzjl);
	
    double fr2il=2.0*f5lilj;
    double fr2jl=-(f2lilj*f3jl+2.0*f4lilj);
	
    double fr2x=-vorfak*(fr2il*dxil+fr2jl*dxjl);
    double fr2y=-vorfak*(fr2il*dyil+fr2jl*dyjl);
    double fr2z=-vorfak*(fr2il*dzil+fr2jl*dzjl);
	
    double fr3il=-(f2lilj*f3il+2.0*f4ljli);
    double fr3jl=2.0*f5lilj;
	
    double fr3x=-vorfak*(fr3il*dxil+fr3jl*dxjl);
    double fr3y=-vorfak*(fr3il*dyil+fr3jl*dyjl);
    double fr3z=-vorfak*(fr3il*dzil+fr3jl*dzjl);

	
    // return force...
    at0.fx+=frx;
    at0.fy+=fry;
    at0.fz+=frz;
    at1.fx+=fr2x;
    at1.fy+=fr2y;
    at1.fz+=fr2z;
    at2.fx+=fr3x;
    at2.fy+=fr3y;
    at2.fz+=fr3z;
	
    at0.potential+=ep;
    at1.potential+=ep;
    at2.potential+=ep;

}


// FIXME! no virial calculation...
void StillingerWeber::Compute(Atom& at, Atom& to) {
    OMD_FLOAT dxjl, dyjl, dzjl;
    OMD_FLOAT dxil, dyil, dzil;
    int iat,ito,nidx,lstart,lend;
	
    OMD_FLOAT RRJL=CalcSqrDistance(at,to,dxjl,dyjl,dzjl);
    if(RRJL<CutRadiusSqr) {
        OMD_FLOAT rjl=sqrt(RRJL);
        TwoBodyTerm(at,to,rjl,dxjl,dyjl,dzjl);
		
        Verlet->GetIterationVariables(iat,ito,nidx,lstart,lend);
        for(int i=nidx+1;i<lend;i++) {
            int ka=Verlet->GetNeighbor(i);
            OMD_FLOAT RRIL=CalcSqrDistance(Atoms(ka),at,dxil,dyil,dzil);
            if(RRIL<CutRadiusSqr) {
                ThreeBodyTerm(at,to,Atoms(ka),
                              rjl,-dxjl,-dyjl,-dzjl,
                              sqrt(RRIL),dxil,dyil,dzil);
            }
        }
 
		if(force_eval) force_eval->EvaluateForce(at,to,dxjl,dyjl,dzjl,0.0,0.0,this); 
	}
}
