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
}

void StillingerWeber::Init(MDSystem* WorkSys) {
	ForceKernel::Init(WorkSys);
	
	Verlet=dynamic_cast<VerletList*>(System->GetIterator());
	assert(Verlet->type_of("verlet list"), 
		   "StillingerWeber potential kernel needs (VERLET LIST) iterator");
	
	CutRadius=alpha*sigma;
	CutRadiusSqr=CutRadius*CutRadius;
	K=new OMD_FLOAT[9];
	K[0]=A*eps*p*B*pow(sigma,p);
	K[1]=A*eps*q*pow(sigma,q);
	K[2]=A*eps*B*pow(sigma,p+1);
	K[3]=A*eps*pow(sigma,q+1);
	K[4]=A*eps*B*pow(sigma,p);
	K[5]=A*eps*pow(sigma,q);
	K[6]=sigma*gamma;
	K[7]=lambda*eps;
	K[8]=2.0*K[7];
}

StillingerWeber::~StillingerWeber() {
	delete[] K;
}

void StillingerWeber::PrintInfo(ostream& ost) {
	ost << "id."<<id<<" "<<get_name()<<"\n"
		<< "atoms "<<System->SystemAtoms[A]->get_name()
		<< "<->"<< System->SystemAtoms[B]->get_name()<< "\n"
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
									OMD_FLOAT r1, // at <- first_nb
									OMD_FLOAT dx1,
									OMD_FLOAT dy1,
									OMD_FLOAT dz1,
									OMD_FLOAT r2, // at <- second_nb (iterated)
									OMD_FLOAT dx2,
									OMD_FLOAT dy2,
									OMD_FLOAT dz2)
{
	OMD_FLOAT irr[]   = {1.0/(r1*r1), 1.0/(r2*r2)};
	OMD_FLOAT ira[]   = {1.0/(r1-CutRadius), 1.0/(r2-CutRadius)};
	OMD_FLOAT sgira[] = {K[6]*ira[0], K[6]*ira[1]};
	OMD_FLOAT sgirra[]= {sgira[0]*ira[0]/r1, sgira[1]*ira[1]/r2};
	OMD_FLOAT esgira[]= {exp(sgira[0]), exp(sgira[1])};	
	OMD_FLOAT irr12   = 1.0/(r1*r2);
	OMD_FLOAT costh   = (dx1*dx2+dy1*dy2+dz1*dz2)*irr12 + 1.0/3.0;
	OMD_FLOAT efac    = esgira[0]*esgira[1];
	OMD_FLOAT radfac  = K[7]*efac*costh*costh;
	OMD_FLOAT fac[]   = {radfac*sgirra[0],radfac*sgirra[1]};
	OMD_FLOAT fang    = K[8]*efac*costh;
	OMD_FLOAT fang12  = fang*irr12;
	OMD_FLOAT cosfang = fang*costh;
	OMD_FLOAT cosf[]  = {irr[0]*cosfang,irr[1]*cosfang};
	OMD_FLOAT fx[]    = {dx1*(fac[0]+cosf[0])-dx2*fang12, dx1*(fac[1]+cosf[1])-dx2*fang12};
	OMD_FLOAT fy[]    = {dy1*(fac[0]+cosf[0])-dy2*fang12, dy1*(fac[1]+cosf[1])-dy2*fang12};
	OMD_FLOAT fz[]    = {dz1*(fac[0]+cosf[0])-dz2*fang12, dz1*(fac[1]+cosf[1])-dz2*fang12};	

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
	radfac/=3.0;
	at0.potential+=radfac;
	at1.potential+=radfac;
	at2.potential+=radfac;	
}


// FIXME! no virial calculation...
void StillingerWeber::Compute(Atom& at, Atom& to) {
	OMD_FLOAT  dx, dy, dz;
	OMD_FLOAT kdx,kdy,kdz;
	int iat,ito,nidx,lstart,lend;

	OMD_FLOAT RR=CalcSqrDistance(at,to,dx,dy,dz);
	if(RR<=CutRadiusSqr) {
		OMD_FLOAT r=sqrt(RR);
		TwoBodyTerm(at,to,r,dx,dy,dz);
		Verlet->GetIterationVariables(iat,ito,nidx,lstart,lend);
		for(int i=nidx+1;i<lend;i++) {
			int ka=Verlet->GetNeighbor(i);
			OMD_FLOAT KRR=CalcSqrDistance(Atoms(ka),at,kdx,kdy,kdz); // sign inverted!
			if(KRR<=CutRadius) {
				dx=-dx;dy=-dy;dz=-dz;
				ThreeBodyTerm(at,to,Atoms(ka),r,-dx,-dy,-dz,sqrt(KRR),kdx,kdy,kdz);
			}
		}
		if(force_eval) force_eval->EvaluateForce(at,to,dx,dy,dz,0.0,0.0,this);
	}
}
