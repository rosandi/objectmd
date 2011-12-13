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
#include <conditioner/VerletListFull.hpp>
#include <potential/sw.hpp>

OMD_FLOAT onethird;

int nzz, nspl, nsw;

StillingerWeber::StillingerWeber(string material) {
	paramfile=material;
	set_name("STILLINGER WEBER");
	register_class(get_name());
}

void StillingerWeber::ReadParameter() {
	param.read(search_path("$OMD_TABLE", "sw."+paramfile));
	eps=param.double_value("sw.epsilon");
	bigA=param.double_value("sw.A");
	bigB=param.double_value("sw.B");
	powerp=param.double_value("sw.p");
	powerq=param.double_value("sw.q");
	sigma=param.double_value("sw.sigma");
	alpha=param.double_value("sw.a");
	lambda=param.double_value("sw.lambda");
	gamma=param.double_value("sw.gamma");
	icos0=param.double_value("sw.invers_cos0");
	
	if(param.exist("zbl")) {
		param.peek("zbl", usezbl);
		if(usezbl) {
			sx0=param.double_value("spline.range", 0);
			sx1=param.double_value("spline.range", 1);
			sa0=param.double_value("spline.coeff", 0);
			sa1=param.double_value("spline.coeff", 1);
			sa2=param.double_value("spline.coeff", 2);
			sa3=param.double_value("spline.coeff", 3);			
		}
	}

}

void StillingerWeber::Init(MDSystem* WorkSys) {
	ForceKernel::Init(WorkSys);
	
	Verlet=dynamic_cast<VerletListFull*>(System->GetIterator());
	assert(Verlet->type_of("verlet list full neighbor"), 
		   "StillingerWeber potential kernel needs (VERLET LIST FULL NEIGHBOR) iterator");

	CutRadius=alpha*sigma;
	CutRadiusSqr=CutRadius*CutRadius;
	c1 = bigA*eps*powerp*bigB*pow(sigma,powerp);
  c2 = bigA*eps*powerq*pow(sigma,powerq);
	c3 = bigA*eps*bigB*pow(sigma,powerp+1.0);
	c4 = bigA*eps*pow(sigma,powerq+1.0);
  c5 = bigA*eps*bigB*pow(sigma,powerp);
  c6 = bigA*eps*pow(sigma,powerq);
	c7 = sigma*gamma;
	c8 = lambda*eps;
	c9 = 2.0*lambda*eps;
	c10= 1.0/icos0;
	
	onethird=1./3.;
	
	OMD_FLOAT znum1=System->SystemAtoms[AtomTypeA]->Z;
	OMD_FLOAT znum2=System->SystemAtoms[AtomTypeB]->Z;
	OMD_FLOAT znor=0.468521928/(pow(znum1,0.23)+pow(znum2, 0.23));
	zc =14.39964415*znum1*znum2;
	ze0=3.2000/znor;
	ze1=0.9423/znor;
	ze2=0.4029/znor;
	ze3=0.2016/znor;
	
}

void StillingerWeber::PrintInfo(ostream& ost) {
	ost << "id."<<id<<" "<<get_name()<<" -- "
		<< "atoms ("<<System->SystemAtoms[AtomTypeA]->get_name()
		<< "<->"<< System->SystemAtoms[AtomTypeB]->get_name()
		<<"); epsilon="<<eps
		<<"; A="<<bigA
		<<"; B="<<bigB
		<<"; p="<<powerp
		<<"; q="<<powerq
		<<"; sigma="<<sigma
		<<"; a="<<alpha
		<<"; lambda="<<lambda
		<<"; gamma="<<gamma
		<<"; cos0="<<1./icos0;
	
	if(usezbl) {
		ost << "; spline to ZBL at ("<<sx0<<","<<sx1
		<< ") coeff ("<<sa0<<" "<<sa1<<" "<<sa2<<" "<<sa3<<")";
	}
	
	ost << std::endl;
}

void StillingerWeber::ZBL(Atom& at, Atom& to, 
						  OMD_FLOAT r, OMD_FLOAT dx, OMD_FLOAT dy, OMD_FLOAT dz)
{
	
	OMD_FLOAT fex[]= {
		0.1818 *exp(-ze0*r),
		0.5099 *exp(-ze1*r),
		0.2802 *exp(-ze2*r),
		0.02817*exp(-ze3*r)
	};
	
	OMD_FLOAT fpot=(fex[0]+fex[1]+fex[2]+fex[3])*zc/r;
	OMD_FLOAT fr=((1.0/r+ze0)*fex[0]+
				  (1.0/r+ze1)*fex[1]+
				  (1.0/r+ze2)*fex[2]+
				  (1.0/r+ze3)*fex[3])*zc/(r*r); // -del fpot
	
	// return force...
	at.fx+=dx*fr;
	at.fy+=dy*fr;
	at.fz+=dz*fr;
	to.fx-=dx*fr;
	to.fy-=dy*fr;
	to.fz-=dz*fr;
	fpot*=0.5;
	at.potential+=fpot;
	to.potential+=fpot;	
	
}

void StillingerWeber::Spline(Atom& at, Atom& to,
							 OMD_FLOAT r, OMD_FLOAT dx, OMD_FLOAT dy, OMD_FLOAT dz)
{
	OMD_FLOAT fr=sa1+sa2*r+sa3*r*r;
	OMD_FLOAT fpot=0.5*(sa0+fr*r);

	// return force...
	at.fx-=dx*fr;
	at.fy-=dy*fr;
	at.fz-=dz*fr;
	to.fx+=dx*fr;
	to.fy+=dy*fr;
	to.fz+=dz*fr;
	at.potential+=fpot;
	to.potential+=fpot;	
	
}

void StillingerWeber::TwoBodyTerm(Atom& at, Atom& to,
								  OMD_FLOAT r,
								  OMD_FLOAT dx,
								  OMD_FLOAT dy,
								  OMD_FLOAT dz)
{
	OMD_FLOAT ira=1.0/(r-CutRadius);
	OMD_FLOAT rp=pow(r,-powerp);
	OMD_FLOAT rq=pow(r,-powerq);
	OMD_FLOAT fex=exp(sigma*ira);
	OMD_FLOAT fr=(c1*rp-c2*rq + (c3*rp-c4*rq)*ira*ira*r)*fex/(r*r);
	OMD_FLOAT fpot=0.5*(c5*rp-c6*rq)*fex;

	// return force...
	at.fx+=dx*fr;
	at.fy+=dy*fr;
	at.fz+=dz*fr;
	to.fx-=dx*fr;
	to.fy-=dy*fr;
	to.fz-=dz*fr;
	at.potential+=fpot;
	to.potential+=fpot;

}

void StillingerWeber::ThreeBodyTerm(Atom& at0, Atom& at1, Atom& at2,
									OMD_FLOAT r1, // at <- first_nb
									OMD_FLOAT dx1,
									OMD_FLOAT dy1,
									OMD_FLOAT dz1,
									OMD_FLOAT r2, // at <- second_nb (iterated)
									OMD_FLOAT dx2,
									OMD_FLOAT dy2,
									OMD_FLOAT dz2)
{
		
	OMD_FLOAT ira[]={1.0/(r1-CutRadius), 1.0/(r2-CutRadius)};
	OMD_FLOAT fira[]= {c7*ira[0], c7*ira[1]};
	OMD_FLOAT ir1r2 = 1.0/(r1*r2);
	OMD_FLOAT fex=exp(fira[0])*exp(fira[1]);	
	OMD_FLOAT costheta=(dx1*dx2+dy1*dy2+dz1*dz2)*ir1r2;
	OMD_FLOAT dcos = costheta-c10;
	OMD_FLOAT dccos = dcos*dcos;
	OMD_FLOAT ftheta = c9*fex*dcos;
	OMD_FLOAT ft12=ir1r2*ftheta;
	
	OMD_FLOAT fpot=c8*fex*dccos;
	OMD_FLOAT fcos=costheta*ftheta;
	OMD_FLOAT fr[]={fpot*fira[0]*ira[0]/r1,fpot*fira[1]*ira[1]/r2};
	OMD_FLOAT ft[]={fcos/(r1*r1), fcos/(r2*r2)};
	OMD_FLOAT fx[]={dx1*(fr[0]+ft[0])-dx2*ft12,dx2*(fr[1]+ft[1])-dx1*ft12};
	OMD_FLOAT fy[]={dy1*(fr[0]+ft[0])-dy2*ft12,dy2*(fr[1]+ft[1])-dy1*ft12};
	OMD_FLOAT fz[]={dz1*(fr[0]+ft[0])-dz2*ft12,dz2*(fr[1]+ft[1])-dz1*ft12};
	
	// return force...
	at0.fx-=fx[0]+fx[1];
	at0.fy-=fy[0]+fy[1];
	at0.fz-=fz[0]+fz[1];
	at1.fx+=fx[0];
	at1.fy+=fy[0];
	at1.fz+=fz[0];
	at2.fx+=fx[1];
	at2.fy+=fy[1];
	at2.fz+=fz[1];

	fpot*=onethird;
	at0.potential+=fpot;
	at1.potential+=fpot;
	at2.potential+=fpot;	
	
}

// FIXME! no virial calculation...

void StillingerWeber::ComputeHalf(Atom& at, Atom& to) {
	OMD_FLOAT  dx, dy, dz;
	OMD_FLOAT RR=CalcSqrDistance(at,to,dx,dy,dz);
	if(RR<CutRadiusSqr) {
		OMD_FLOAT r=sqrt(RR);
		if(usezbl) {
			if(r<sx0)      ZBL(at,to,r,dx,dy,dz);  // ZBL region
			else if(r<sx1) Spline(at,to,r,dx,dy,dz); // spline region
			else           TwoBodyTerm(at,to,r,dx,dy,dz); // sw region
		}
		else TwoBodyTerm(at,to,r,dx,dy,dz);
	}
}

void StillingerWeber::ComputeFull(Atom& at, Atom& to) {
	OMD_FLOAT  dx, dy, dz;
	OMD_FLOAT kdx,kdy,kdz;
	
	int iat,ito,nidx;
	NeighborList* nlist;
	Verlet->GetIterationVariables(iat,ito,nidx,nlist);
	if(ito==nlist->end-1) return;

	OMD_FLOAT RR=CalcSqrDistance(to,at,dx,dy,dz);
	if(RR<CutRadiusSqr) {
		for(int i=nidx+1;i<nlist->end;i++) {
			int ka=nlist->list[i];
			OMD_FLOAT KRR=CalcSqrDistance(Atoms(ka),at,kdx,kdy,kdz);
			if(KRR<CutRadiusSqr)
				ThreeBodyTerm(at,to,Atoms(ka),sqrt(RR),dx,dy,dz,sqrt(KRR),kdx,kdy,kdz);
		}
		if(force_eval) force_eval->EvaluateForce(at,to,dx,dy,dz,0.0,0.0,this);
	}
}
