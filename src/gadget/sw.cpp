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
	OMD_FLOAT rinvsq = 1.0/(r*r);
	OMD_FLOAT rp = pow(r,-powerp);
	OMD_FLOAT rq = pow(r,-powerq);
	OMD_FLOAT rainv = 1.0 / (r - CutRadius);
	OMD_FLOAT rainvsq = rainv*rainv*r;
	OMD_FLOAT expsrainv = exp(sigma * rainv);
	OMD_FLOAT fr = (c1*rp - c2*rq +
			  (c3*rp - c4*rq) * rainvsq) * expsrainv * rinvsq;
	OMD_FLOAT ep = 0.5*(c5*rp - c6*rq) * expsrainv;

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
	OMD_FLOAT rinvsq1 = 1.0/(r1*r1);
	OMD_FLOAT rainv1 = 1.0/(r1 - CutRadius);
	OMD_FLOAT gsrainv1 = c7 * rainv1;
	OMD_FLOAT gsrainvsq1 = gsrainv1*rainv1/r1;
	OMD_FLOAT expgsrainv1 = exp(gsrainv1);
	
	OMD_FLOAT rinvsq2 = 1.0/(r2*r2);
	OMD_FLOAT rainv2 = 1.0/(r2 - CutRadius);
	OMD_FLOAT gsrainv2 = c7 * rainv2;
	OMD_FLOAT gsrainvsq2 = gsrainv2*rainv2/r2;
	OMD_FLOAT expgsrainv2 = exp(gsrainv2);
	
	OMD_FLOAT rinv12 = 1.0/(r1*r2);
	OMD_FLOAT cs = (dx1*dx2 + dy1*dy2 + dz1*dz2) * rinv12;
	OMD_FLOAT delcs = cs + 1.0/3.0;
	OMD_FLOAT delcssq = delcs*delcs;

	OMD_FLOAT facexp = expgsrainv1*expgsrainv2;
	
	OMD_FLOAT facrad = c8 * facexp*delcssq;
	OMD_FLOAT frad1 = facrad*gsrainvsq1;
	OMD_FLOAT frad2 = facrad*gsrainvsq2;
	OMD_FLOAT facang = c9 * facexp*delcs;
	OMD_FLOAT facang12 = rinv12*facang;
	OMD_FLOAT csfacang = cs*facang;
	OMD_FLOAT csfac1 = rinvsq1*csfacang;
	OMD_FLOAT fj[3], fk[3];

	fj[0] = dx1*(frad1+csfac1)-dx2*facang12;
	fj[1] = dy1*(frad1+csfac1)-dy2*facang12;
	fj[2] = dz1*(frad1+csfac1)-dz2*facang12;
	
	OMD_FLOAT csfac2 = rinvsq2*csfacang;
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
	
}

// FIXME! no virial calculation...

void StillingerWeber::ComputeHalf(Atom& at, Atom& to) {
	OMD_FLOAT  dx, dy, dz;
	OMD_FLOAT RR=CalcSqrDistance(at,to,dx,dy,dz);
	if(RR<CutRadiusSqr) TwoBodyTerm(at,to,sqrt(RR),dx,dy,dz);
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
