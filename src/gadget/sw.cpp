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
	bigA=param.double_value("A");
	bigB=param.double_value("B");
	powerp=param.double_value("p");
	powerq=param.double_value("q");
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
	c1 = bigA*eps*powerp*bigB*pow(sigma,powerp);
    c2 = bigA*eps*powerq*pow(sigma,powerq);
	c3 = bigA*eps*bigB*pow(sigma,powerp+1.0);
	c4 = bigA*eps*pow(sigma,powerq+1.0);
    c5 = bigA*eps*bigB*pow(sigma,powerp);
    c6 = bigA*eps*pow(sigma,powerq);
	
	sigma_gamma=sigma*gamma;
	lambda_eps=lambda*eps;
	lambda_eps2=2.0*lambda*eps;
}

StillingerWeber::~StillingerWeber() {
}

void StillingerWeber::PrintInfo(ostream& ost) {
	ost << "id."<<id<<" "<<get_name()<<"\n"
		<< "atoms "<<System->SystemAtoms[AtomTypeA]->get_name()
		<< "<->"<< System->SystemAtoms[AtomTypeB]->get_name()<< "\n"
		<<" parameter:"
		<<"\n epsilon="<<eps
		<<"\n A="<<bigA
		<<"\n B="<<bigB
		<<"\n p="<<powerp
		<<"\n q="<<powerq
		<<"\n sigma="<<sigma
		<<"\n a="<<alpha
		<<"\n lambda="<<lambda
		<<"\n gamma="<<gamma<<"\n";
}

void StillingerWeber::TwoBodyTerm(Atom& at, Atom& to,
								  OMD_FLOAT r,
								  OMD_FLOAT dx,
								  OMD_FLOAT dy,
								  OMD_FLOAT dz)
{
	double rinvsq = 1.0/(r*r);
	double rp = pow(r,-powerp);
	double rq = pow(r,-powerq);
	double rainv = 1.0 / (r - CutRadius);
	double rainvsq = rainv*rainv*r;
	double expsrainv = exp(sigma * rainv);
	double fr = (c1*rp - c2*rq +
			  (c3*rp - c4*rq) * rainvsq) * expsrainv * rinvsq;
	double ep = 0.5*(c5*rp - c6*rq) * expsrainv;

	// return force...
	at.fx+=dx*fr;
	at.fy+=dy*fr;
	at.fz+=dz*fr;
	to.fx-=dx*fr;
	to.fy-=dy*fr;
	to.fz-=dz*fr;
 
	at.potential+=ep;
	to.potential+=ep;
	
	std::cerr << "pair("<<at.nid<<","<<to.nid<<")\n";
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
	double rinvsq1 = 1.0/(r1*r1);
	double rainv1 = 1.0/(r1 - CutRadius);
	double gsrainv1 = sigma_gamma * rainv1;
	double gsrainvsq1 = gsrainv1*rainv1/r1;
	double expgsrainv1 = exp(gsrainv1);
	
	double rinvsq2 = 1.0/(r2*r2);
	double rainv2 = 1.0/(r2 - CutRadius);
	double gsrainv2 = sigma_gamma * rainv2;
	double gsrainvsq2 = gsrainv2*rainv2/r2;
	double expgsrainv2 = exp(gsrainv2);
	
	double rinv12 = 1.0/(r1*r2);
	double cs = (dx1*dx2 + dy1*dy2 + dz1*dz2) * rinv12;
	double delcs = cs + 1.0/3.0;
	double delcssq = delcs*delcs;

	double facexp = expgsrainv1*expgsrainv2;
	
	double facrad = lambda_eps * facexp*delcssq;
	double frad1 = facrad*gsrainvsq1;
	double frad2 = facrad*gsrainvsq2;
	double facang = lambda_eps2 * facexp*delcs;
	double facang12 = rinv12*facang;
	double csfacang = cs*facang;
	double csfac1 = rinvsq1*csfacang;
	double fj[3], fk[3];

	fj[0] = dx1*(frad1+csfac1)-dx2*facang12;
	fj[1] = dy1*(frad1+csfac1)-dy2*facang12;
	fj[2] = dz1*(frad1+csfac1)-dz2*facang12;
	
	double csfac2 = rinvsq2*csfacang;
	fk[0] = dx2*(frad2+csfac2)-dx1*facang12;
	fk[1] = dy2*(frad2+csfac2)-dy1*facang12;
	fk[2] = dz2*(frad2+csfac2)-dz1*facang12;

	// return force...
	at0.fx-=fj[0]+fk[0];
	at0.fy-=fj[1]+fk[1];
	at0.fz-=fj[2]+fk[2];
	at1.fx+=fj[0];
	at1.fy+=fj[1];
	at1.fz+=fj[2];
	at2.fx+=fk[0];
	at2.fy+=fk[1];
	at2.fz+=fk[2];

	facrad/=3.0;
	at0.potential+=facrad;
	at1.potential+=facrad;
	at2.potential+=facrad;	
	
	std::cerr << "ang("<<at0.nid<<","<<at1.nid<<","<<at2.nid<<")\n";
}

void StillingerWeber::ComputeHalf(Atom& at, Atom& to) {
	OMD_FLOAT  dx, dy, dz;
	OMD_FLOAT kdx,kdy,kdz;
	int iat,ito,nidx,lstart,lend;

	OMD_FLOAT RR=CalcSqrDistance(at,to,dx,dy,dz);
	if(RR<CutRadiusSqr) {
		TwoBodyTerm(at,to,sqrt(RR),dx,dy,dz);
		OMD_FLOAT RR=CalcSqrDistance(to,at,dx,dy,dz);
		if(RR<CutRadiusSqr) {
			Verlet->GetIterationVariables(iat,ito,nidx,lstart,lend);
			for(int i=nidx+1;i<lend;i++) {
				int ka=Verlet->GetNeighbor(i);
				OMD_FLOAT KRR=CalcSqrDistance(Atoms(ka),at,kdx,kdy,kdz);
				if(KRR<CutRadiusSqr)
					ThreeBodyTerm(at,to,Atoms(ka),sqrt(RR),dx,dy,dz,sqrt(KRR),kdx,kdy,kdz);
			}
		}
		if(force_eval) force_eval->EvaluateForce(at,to,dx,dy,dz,0.0,0.0,this);
	}
}

// FIXME! no virial calculation...
/*
void StillingerWeber::ComputeHalf(Atom& at, Atom& to) {
	OMD_FLOAT  dx, dy, dz;
	OMD_FLOAT RR=CalcSqrDistance(at,to,dx,dy,dz);
	if(RR<CutRadiusSqr) TwoBodyTerm(at,to,sqrt(RR),dx,dy,dz);
}

void StillingerWeber::ComputeFull(Atom& at, Atom& to) {
	OMD_FLOAT  dx, dy, dz;
	OMD_FLOAT kdx,kdy,kdz;
	int iat,ito,nidx,lstart,lend;
	
	OMD_FLOAT RR=CalcSqrDistance(to,at,dx,dy,dz);
	if(RR<CutRadiusSqr) {
		Verlet->GetIterationVariables(iat,ito,nidx,lstart,lend);
		for(int i=nidx+1;i<lend;i++) {
			int ka=Verlet->GetNeighbor(i);
			OMD_FLOAT KRR=CalcSqrDistance(Atoms(ka),at,kdx,kdy,kdz);
			if(KRR<CutRadiusSqr)
				ThreeBodyTerm(at,to,Atoms(ka),sqrt(RR),dx,dy,dz,sqrt(KRR),kdx,kdy,kdz);
		}
		if(force_eval) force_eval->EvaluateForce(at,to,dx,dy,dz,0.0,0.0,this);
	}
}
*/
