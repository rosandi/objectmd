/* omdlib
 ********************************************************
 *
 * ObjectMD
 * Molecular Dynamics Class Library
 * Version 2.0 (2009)
 *
 * Released under:
 *       GNU GENERAL PUBLIC LICENSE Version 3
 *       Copyright (C) 2005  Yudi Rosandi
 *
 * Read license.txt file in the source's root directory.
 *
 ********************************************************
 * ObjectMD header file
 *
 * Morse potential implementation
 *
 *
*/

// FIXME! update time?ps

#include <cstdio>
#include <cmath>
#include <omd/forcekernel.hpp>
#include <omd/conditioner.hpp>
#include <omd/iterator.hpp>

class ClassNeighborList: public Calc_Conditioner {
friend class Tersoff;

	OMD_SIZET **list;
	OMD_INT*  length;
	OMD_INT*  alloc_size;
	OMD_SIZET num_atom;
	OMD_FLOAT cut_radius;

public:

	ClassNeighborList(OMD_FLOAT cut) {
		list=NULL;
		length=NULL;
		alloc_size=NULL;
		num_atom=0;
		cut_radius=cut;
		register_class("neighborlist_generator");
		set_name("NEIGHBORLIST_GENERATOR");
	}

	virtual ~ClassNeighborList() {
		for(OMD_SIZET i=0;i<num_atom;i++) MemFree(list[i]);
		MemFree(list);
		MemFree(length);
		MemFree(alloc_size);
	}

	void Init(MDSystem* WorkSys) {
		Calc_Conditioner::Init(WorkSys);
	}

	void realloc_array(OMD_SIZET size) {

		if(size<=num_atom) return;
		MemRealloc(list, size*sizeof(int*));
		MemRealloc(length, size*sizeof(int));
		MemRealloc(alloc_size, size*sizeof(int));

		for(OMD_SIZET i=num_atom; i<size; i++){
			list[i]=NULL;
			alloc_size[i]=0;
		}

		for(OMD_INT i=0;i<size;i++) length[i]=0;
		num_atom=size;

	}

	void grow(OMD_SIZET idx, OMD_SIZET size) {
		if(size<=alloc_size[idx]) return;
		MemRealloc(list[idx], size*sizeof(int));
	}

	void IterationNode(OMD_SIZET at, OMD_SIZET to) {
		if(CalcDistance(at,to)<System->GetMaxCutRadius()) {
			grow(at,length[at]+1);
//			grow(to,length[to]+1);
			list[at][length[at]++]=to;
//			list[to][length[to]++]=at;
		}
	}

	void PreCalculation() {
		realloc_array(GetNAtom());
		Iterator->Iterate(this);
	}

};

/**
 * @f V_{ij} = f_C (r_{ij}) \left[ f_R(r_{ij} + b_{ij} f_A(r_{ij}) \right]
*/

class Tersoff: public ForceKernel {

	string material_name;

	//-- tersoff parameters --
	OMD_FLOAT a;
	OMD_FLOAT b;
	OMD_FLOAT c2;
	OMD_FLOAT d2;
	OMD_FLOAT h;
	OMD_FLOAT n;
	OMD_FLOAT De;
	OMD_FLOAT Re;
	OMD_FLOAT alpha;
	OMD_FLOAT beta;
	OMD_FLOAT lambda_1;
	OMD_FLOAT lambda_2;
	OMD_FLOAT lambda_3;

	//-- zbl parameters ------
	OMD_FLOAT zbl_factor;
	OMD_FLOAT zbl_norm;
	OMD_FLOAT zbl_cut;
	OMD_FLOAT sharpness;
	OMD_FLOAT ze0, ze1, ze2, ze3;
	OMD_FLOAT spline_range[2];

	bool connect_zbl;
	ClassNeighborList *NeighborList;

public:	

	Tersoff(string s_material="-") {
		register_class("tersoff_bond_order");
		set_name("TERSOFF_BOND_ORDER");
		material_name.assign(s_material);
		connect_zbl=false;
	}

	void Init(MDSystem *WorkSys) {
		ForceKernel::Init(WorkSys);

		SysParam->peek("tersoff.material", material_name);
		string paramfile=System->search_path("$OMD_TABLE", "tersoff."+material_name);
		assert(file_exist(paramfile), "can not find potenial parameters for "+material_name);
		param.read(paramfile);

		// put file path information in parameter
		param.append("table_file_path "+paramfile);

		// read parameters...
		a=param.double_value("tersoff.A");
		b=param.double_value("tersoff.B");
		c2=param.double_value("tersoff.c");
		d2=param.double_value("tersoff.d");
		h=param.double_value("tersoff.h");
		n=param.double_value("tersoff.n");
		De=param.double_value("tersoff.D");
		Re=param.double_value("tersoff.R");
		alpha=param.double_value("tersoff.alpha");
		beta=param.double_value("tersoff.beta");
		lambda_1=param.double_value("tersoff.lambda_1");
		lambda_2=param.double_value("tersoff.lambda_2");
		lambda_3=param.double_value("tersoff.lambda_3");

		if(param.exist("tersoff.zbl")) {
			if(param.string_value("tersoff.zbl")=="true") {
				connect_zbl=true;
				zbl_cut=param.double_value("tersoff.zbl_cut");
				sharpness=param.double_value("tersoff.zbl_sharpness");
				spline_range[0]=param.double_value("tersoff.spline_range", 0);
				spline_range[1]=param.double_value("tersoff.spline_range", 1);
			}
		}

		// cancel zbl if stated in system's parameter...
		SysParam->peek("tersoff.zbl",connect_zbl);

		c2*=c2;
		d2*=d2;
		CutRadius=De+Re;
		NeighborList=new ClassNeighborList(CutRadius);
		System->AddConditioner(NeighborList);

		OMD_FLOAT z1=System->SystemAtoms.at(A)->Z;
		OMD_FLOAT z2=System->SystemAtoms.at(B)->Z;

		zbl_factor=14.40*z1*z2;
		zbl_norm=0.46852/(pow(z1,0.23)+pow(z2,0.23));

        ze0=3.2/zbl_norm;
        ze1=0.9423/zbl_norm;
        ze2=0.4029/zbl_norm;
        ze3=0.2016/zbl_norm;

	}

	void ComputeConnect(OMD_SIZET at, OMD_SIZET to, OMD_FLOAT &r,
			            OMD_FLOAT &dx, OMD_FLOAT &dy, OMD_FLOAT &dz,
			            OMD_FLOAT f_at[], OMD_FLOAT &pat, OMD_FLOAT f_to[], OMD_FLOAT &pto) {
		// good numbers: sharpness=30, zbl_cut=1
		OMD_FLOAT fexp=-sharpness*(r-zbl_cut);
		OMD_FLOAT ff=1.0/(1.0+fexp);
		OMD_FLOAT dff_dr=-sharpness*ff*ff*fexp;
		OMD_FLOAT zpat, zpto, zfat[3], zfto[3], dpdr;
		OMD_FLOAT r_hat[]={dx/r, dy/r, dz/r};
/*
		ComputeZBL(r,dx,dy,dz,zfat,zpat,zfto,zpto);
		ComputeTersoff(at,to,r,dx,dy,dz,f_at,pat,f_to,pto);

		// ------- zbl --------------- tersoff ---
		f_at[0]=(-dff_dr*zpat*r_hat[0] + (1-ff)*zfat[0]) + (-dff_dr*pat*r_hat[0] + ff*f_at[0]);
		f_at[1]=(-dff_dr*zpat*r_hat[1] + (1-ff)*zfat[1]) + (-dff_dr*pat*r_hat[1] + ff*f_at[1]);
		f_at[2]=(-dff_dr*zpat*r_hat[2] + (1-ff)*zfat[2]) + (-dff_dr*pat*r_hat[2] + ff*f_at[2]);

		f_to[0]=(dff_dr*zpto*r_hat[0] + (1-ff)*zfto[0]) - (dff_dr*pat*r_hat[0] - ff*f_at[0]);
		f_to[1]=(dff_dr*zpto*r_hat[1] + (1-ff)*zfto[1]) - (dff_dr*pat*r_hat[1] - ff*f_at[1]);
		f_to[2]=(dff_dr*zpto*r_hat[2] + (1-ff)*zfto[2]) - (dff_dr*pat*r_hat[2] - ff*f_at[2]);

		pat =(1-ff)*zpat+ff*pat;
		pto =(1-ff)*zpto+ff*pto;
*/
	}

	void ComputeZBL(OMD_FLOAT &r, OMD_FLOAT &dx, OMD_FLOAT &dy, OMD_FLOAT &dz,
			        OMD_FLOAT f_at[], OMD_FLOAT &pat, OMD_FLOAT f_to[], OMD_FLOAT &pto) {
		OMD_FLOAT zexp0=0.1818 *exp(-ze0*r);
		OMD_FLOAT zexp1=0.5099 *exp(-ze1*r);
		OMD_FLOAT zexp2=0.2802 *exp(-ze2*r);
		OMD_FLOAT zexp3=0.02817*exp(-ze3*r);
		OMD_FLOAT zr=zbl_factor/r;
		OMD_FLOAT inr=1.0/r;
		OMD_FLOAT fr=(inr+ze0)*zexp0+(inr+ze1)*zexp1+(inr+ze2)*zexp2+(inr+ze3)*zexp3;
		pto=pat=(zexp0+zexp1+zexp2+zexp3)*zr;

		f_at[0]= fr*zr*dx;
		f_at[1]= fr*zr*dy;
		f_at[2]= fr*zr*dz;
		f_to[0]=-fr*zr*dx;
		f_to[1]=-fr*zr*dy;
		f_to[2]=-fr*zr*dz;
	}

	OMD_FLOAT CutFunction(OMD_FLOAT r) {
		if(r<(Re-De)) return 1.0;
		if(r<(Re+De)) return 0.5-0.5*sin(M_PI*(r-Re)/(2.0*De));
		return 0.0;
	}

	OMD_FLOAT dCutFunction(OMD_FLOAT r){
		if(r<(Re-De)) return 0.0;
		if(r<(Re+De)) return -M_PI/(4.0*De)*cos(M_PI*(r-Re)/(2.0*De));
		return 0.0;
	}

	void ComputeBondOrder(OMD_INT at, OMD_INT to, OMD_FLOAT rij,
						  OMD_FLOAT dx, OMD_FLOAT dy, OMD_FLOAT dz) {

		OMD_SIZET* list=NeighborList->list[at];
		OMD_SIZET  nn=NeighborList->length[at];;
		OMD_FLOAT g_cos_theta, dg_cos_theta;
		OMD_FLOAT fc_rij, dfc_rij;

		OMD_FLOAT dr2=dx*dx+dy*dy+dz*dz;
		OMD_FLOAT rij_hat[]={dx/rij,dy/rij,dz/rij};
		OMD_FLOAT zeta=0.0, dzeta[]={0.,0,.0};
		OMD_FLOAT dtmpx[nn], dtmpy[nn], dtmpz[nn];
		OMD_SIZET nj;

		for(OMD_INT i=0;i<nn;i++) {
			dtmpx[i]=dtmpy[i]=dtmpz[i]=0.0;
			if(list[i]==to) nj=i;
		}

		for(OMD_INT i=0; i<nn; i++) {

			if(list[i]==to) continue;

			OMD_FLOAT ndx,ndy,ndz;
			OMD_FLOAT rik=sqrt(CalcSqrDistance(at, list[i], ndx, ndy, ndz));
			OMD_FLOAT rik_hat[]={ndx/rik,ndy/rik,ndz/rik};
			OMD_FLOAT rijxrik=rij*rik;
			OMD_FLOAT costheta=(dx*ndx+dy*ndy+dz*ndz)/rijxrik;
			OMD_FLOAT g_denom=(d2+(h-costheta)*(h-costheta));

			OMD_FLOAT fc=CutFunction(rik);
			OMD_FLOAT dfc_dr=dCutFunction(rik);
			OMD_FLOAT gtheta=1+c2/d2-c2/g_denom;
			OMD_FLOAT dgtheta_dcostheta=-2.*c2*(h-costheta)/(g_denom*g_denom);
			OMD_FLOAT fexp=exp(pow(lambda_3*(rij-rik),3.0));
			OMD_FLOAT dfexp_dr=3.*lambda_3*lambda_3*lambda_3*(rij-rik)*(rij-rik)*fexp;

			OMD_FLOAT dxj_costheta=(rik_hat[0]-costheta*rij_hat[0])/rij;
			OMD_FLOAT dyj_costheta=(rik_hat[1]-costheta*rij_hat[1])/rij;
			OMD_FLOAT dzj_costheta=(rik_hat[2]-costheta*rij_hat[2])/rij;
			OMD_FLOAT dxk_costheta=(rij_hat[0]-costheta*rik_hat[0])/rik;
			OMD_FLOAT dyk_costheta=(rij_hat[1]-costheta*rik_hat[1])/rik;
			OMD_FLOAT dzk_costheta=(rij_hat[2]-costheta*rik_hat[2])/rik;

			OMD_FLOAT dcostheta[]= {
					dxj_costheta+dxk_costheta,
					dyj_costheta+dyk_costheta,
					dzj_costheta+dzk_costheta
			};

			OMD_FLOAT eq_dfc=dfc_dr*gtheta*fexp;
			OMD_FLOAT eq_dgtheta=fc*dgtheta_dcostheta*fexp;
			OMD_FLOAT eq_dfexp=fc*gtheta*dfexp_dr;

			dzeta[0]+=eq_dfc*rik_hat[0]+eq_dgtheta*dcostheta[0]+eq_dfexp*(rij_hat[0]-rik_hat[0]);
			dzeta[1]+=eq_dfc*rik_hat[1]+eq_dgtheta*dcostheta[1]+eq_dfexp*(rij_hat[1]-rik_hat[1]);
			dzeta[2]+=eq_dfc*rik_hat[2]+eq_dgtheta*dcostheta[2]+eq_dfexp*(rij_hat[2]-rik_hat[2]);
			dtmpx[nj]+=(-eq_dgtheta*dxj_costheta-eq_dfexp*rij_hat[0]);
			dtmpy[nj]+=(-eq_dgtheta*dyj_costheta-eq_dfexp*rij_hat[1]);
			dtmpz[nj]+=(-eq_dgtheta*dzj_costheta-eq_dfexp*rij_hat[2]);
			dtmpx[i]+=(-eq_dfc*rik_hat[0]-eq_dgtheta*dxk_costheta+eq_dfexp*rik_hat[0]);
			dtmpy[i]+=(-eq_dfc*rik_hat[1]-eq_dgtheta*dyk_costheta+eq_dfexp*rik_hat[1]);
			dtmpz[i]+=(-eq_dfc*rik_hat[2]-eq_dgtheta*dzk_costheta+eq_dfexp*rik_hat[2]);
			zeta += fc*gtheta*fexp;
		}

		OMD_FLOAT bij, att;

		if(zeta !=0.0) {
			OMD_FLOAT b_kern=1+pow(beta*zeta,n);
			bij=pow(b_kern,-1.0/(2.*n));
			OMD_FLOAT db_dzeta=-0.5*pow(b_kern,-(1./(2.*n) + 1.0))*pow(beta,n)*pow(zeta,n-1.0);

			OMD_FLOAT dbij[]={
				db_dzeta*dzeta[0],
				db_dzeta*dzeta[1],
				db_dzeta*dzeta[2]
			};

			// attractive part...
			att=b*exp(-lambda_2*rij);

			Atoms(at).fx+=att*dbij[0];
			Atoms(at).fy+=att*dbij[1];
			Atoms(at).fz+=att*dbij[2];

			// for the neighbors...
			for(int i=0;i<nn;i++){
				Atoms(list[i]).fx+=att*db_dzeta*dtmpx[i];
				Atoms(list[i]).fy+=att*db_dzeta*dtmpy[i];
				Atoms(list[i]).fz+=att*db_dzeta*dtmpz[i];
			}
		} else {
			bij=1.0;
		}

		// repulsive part...
		OMD_FLOAT rep= a*exp(-lambda_1*rij);
		OMD_FLOAT dv_dr=lambda_1*rep-lambda_2*bij*att;
		OMD_FLOAT dvc_dr=CutFunction(rij)*dv_dr+dCutFunction(rij)*(rep-bij*att);

		Atoms(at).fx+=dvc_dr*rij_hat[0];
		Atoms(at).fy+=dvc_dr*rij_hat[1];
		Atoms(at).fz+=dvc_dr*rij_hat[2];
		Atoms(at).potential=Atoms(to).potential=CutFunction(rij)*(rep+bij*att);

	}

/*
	void ComputeTersoff(OMD_SIZET at, OMD_SIZET to, OMD_FLOAT &r,
			            OMD_FLOAT &dx, OMD_FLOAT &dy, OMD_FLOAT &dz,
			            OMD_FLOAT f_at[], OMD_FLOAT &pat, OMD_FLOAT f_to[], OMD_FLOAT &pto) {


		ComputeBondOrder(at,to,r, dx, dy, dz, att, bat, dbat);

		pat=rep+bat*att;
		pto=rep+bto*att;

		f_at[0]=-cut*(drep_dr*r_hat[0] + datt_dr*bat*r_hat[0] + dbat[0]*att)-dcut*pat*r_hat[0];
		f_at[1]=-cut*(drep_dr*r_hat[1] + datt_dr*bat*r_hat[1] + dbat[1]*att)-dcut*pat*r_hat[1];
		f_at[2]=-cut*(drep_dr*r_hat[2] + datt_dr*bat*r_hat[2] + dbat[2]*att)-dcut*pat*r_hat[2];
		f_to[0]= cut*(drep_dr*r_hat[0] + datt_dr*bto*r_hat[0] - dbto[0]*att)+dcut*pat*r_hat[0];
		f_to[1]= cut*(drep_dr*r_hat[1] + datt_dr*bto*r_hat[1] - dbto[1]*att)+dcut*pat*r_hat[1];
		f_to[2]= cut*(drep_dr*r_hat[2] + datt_dr*bto*r_hat[2] - dbto[2]*att)+dcut*pat*r_hat[2];

		pat*=cut;
		pto*=cut;
	}
*/

	void Compute(OMD_SIZET at, OMD_SIZET to) {
		OMD_FLOAT dx, dy, dz;
		OMD_FLOAT r = sqrt(CalcSqrDistance(at, to, dx, dy, dz));

		if(r<CutRadius) {
			OMD_FLOAT pat,pto;
			OMD_FLOAT f_at[3], f_to[3];
			if(connect_zbl) {/*
				if (r<spline_range[0]) ComputeZBL(r,dx,dy,dz,f_at,pat,f_to,pto);
				else if (r<spline_range[1]) ComputeConnect(at,to,r,dx,dy,dz,f_at,pat,f_to,pto);
				else ComputeTersoff(at,to,r,dx,dy,dz,f_at,pat,f_to,pto);
				*/
				die("diabled...");
			} else {
				ComputeBondOrder(at,to,r,dx,dy,dz);
			}
		}
	}

	void   PrintInfo(ostream& ost) {
		ost << "ID." << id << " " << get_name() << "\n";
		ost << "(cut radius = " << CutRadius << ", Table=" << param.string_value("table_file_path") << ")\n";
		ost << "A=" << a << "; B=" << b
			<< "\nlambda_1=" << lambda_1 << "; lambda_2="<<lambda_2
			<< "\nalpha=" << alpha << "; beta="<<beta<<"; n="<<n
			<< "\nc="<<sqrt(c2)<<"; d="<<sqrt(d2)<<"; h="<<h
			<< "\nlambda_3="<<lambda_3<<"; R="<<Re<<"; D="<<De;
		ost << "\n***\n";
	}

};
