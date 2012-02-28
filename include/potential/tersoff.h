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
#include <omd/forcekernel.h>
#include <omd/conditioner.h>
#include <omd/iterator.h>

class ClassNeighborList: public Calc_Conditioner {
friend class Tersoff;

	int **list;
	int*  length;
	int*  alloc_size;
	int num_atom;
	double cut_radius;

public:

	ClassNeighborList(double cut) {
		list=NULL;
		length=NULL;
		alloc_size=NULL;
		num_atom=0;
		cut_radius=cut;
		register_class("neighborlist_generator");
		set_name("NEIGHBORLIST_GENERATOR");
	}

	virtual ~ClassNeighborList() {
		for(int i=0;i<num_atom;i++) MemFree(list[i]);
		MemFree(list);
		MemFree(length);
		MemFree(alloc_size);
	}

	void Init(MDSystem* WorkSys) {
		Calc_Conditioner::Init(WorkSys);
	}

	void realloc_array(int size) {

		if(size<=num_atom) return;
		MemRealloc(list, size*sizeof(int*));
		MemRealloc(length, size*sizeof(int));
		MemRealloc(alloc_size, size*sizeof(int));

		for(int i=num_atom; i<size; i++){
			list[i]=NULL;
			alloc_size[i]=0;
		}

		for(int i=0;i<size;i++) length[i]=0;
		num_atom=size;

	}

	void grow(int idx, int size) {
		if(size<=alloc_size[idx]) return;
		MemRealloc(list[idx], size*sizeof(int));
	}

	void IterationNode(int at, int to) {
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
	double a;
	double b;
	double c2;
	double d2;
	double h;
	double n;
	double De;
	double Re;
	double alpha;
	double beta;
	double lambda_1;
	double lambda_2;
	double lambda_3;

	//-- zbl parameters ------
	double zbl_factor;
	double zbl_norm;
	double zbl_cut;
	double sharpness;
	double ze0, ze1, ze2, ze3;
	double spline_range[2];

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

		double z1=System->SystemAtoms.at(A)->Z;
		double z2=System->SystemAtoms.at(B)->Z;

		zbl_factor=14.40*z1*z2;
		zbl_norm=0.46852/(pow(z1,0.23)+pow(z2,0.23));

        ze0=3.2/zbl_norm;
        ze1=0.9423/zbl_norm;
        ze2=0.4029/zbl_norm;
        ze3=0.2016/zbl_norm;

	}

	void ComputeConnect(int at, int to, double &r,
			            double &dx, double &dy, double &dz,
			            double f_at[], double &pat, double f_to[], double &pto) {
		// good numbers: sharpness=30, zbl_cut=1
		double fexp=-sharpness*(r-zbl_cut);
		double ff=1.0/(1.0+fexp);
		double dff_dr=-sharpness*ff*ff*fexp;
		double zpat, zpto, zfat[3], zfto[3], dpdr;
		double r_hat[]={dx/r, dy/r, dz/r};
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

	void ComputeZBL(double &r, double &dx, double &dy, double &dz,
			        double f_at[], double &pat, double f_to[], double &pto) {
		double zexp0=0.1818 *exp(-ze0*r);
		double zexp1=0.5099 *exp(-ze1*r);
		double zexp2=0.2802 *exp(-ze2*r);
		double zexp3=0.02817*exp(-ze3*r);
		double zr=zbl_factor/r;
		double inr=1.0/r;
		double fr=(inr+ze0)*zexp0+(inr+ze1)*zexp1+(inr+ze2)*zexp2+(inr+ze3)*zexp3;
		pto=pat=(zexp0+zexp1+zexp2+zexp3)*zr;

		f_at[0]= fr*zr*dx;
		f_at[1]= fr*zr*dy;
		f_at[2]= fr*zr*dz;
		f_to[0]=-fr*zr*dx;
		f_to[1]=-fr*zr*dy;
		f_to[2]=-fr*zr*dz;
	}

	double CutFunction(double r) {
		if(r<(Re-De)) return 1.0;
		if(r<(Re+De)) return 0.5-0.5*sin(M_PI*(r-Re)/(2.0*De));
		return 0.0;
	}

	double dCutFunction(double r){
		if(r<(Re-De)) return 0.0;
		if(r<(Re+De)) return -M_PI/(4.0*De)*cos(M_PI*(r-Re)/(2.0*De));
		return 0.0;
	}

	void ComputeBondOrder(int at, int to, double rij,
						  double dx, double dy, double dz) {

		int* list=NeighborList->list[at];
		int  nn=NeighborList->length[at];;
		double g_cos_theta, dg_cos_theta;
		double fc_rij, dfc_rij;

		double dr2=dx*dx+dy*dy+dz*dz;
		double rij_hat[]={dx/rij,dy/rij,dz/rij};
		double zeta=0.0, dzeta[]={0.,0,.0};
		double dtmpx[nn], dtmpy[nn], dtmpz[nn];
		int nj;

		for(int i=0;i<nn;i++) {
			dtmpx[i]=dtmpy[i]=dtmpz[i]=0.0;
			if(list[i]==to) nj=i;
		}

		for(int i=0; i<nn; i++) {

			if(list[i]==to) continue;

			double ndx,ndy,ndz;
			double rik=sqrt(CalcSqrDistance(at, list[i], ndx, ndy, ndz));
			double rik_hat[]={ndx/rik,ndy/rik,ndz/rik};
			double rijxrik=rij*rik;
			double costheta=(dx*ndx+dy*ndy+dz*ndz)/rijxrik;
			double g_denom=(d2+(h-costheta)*(h-costheta));

			double fc=CutFunction(rik);
			double dfc_dr=dCutFunction(rik);
			double gtheta=1+c2/d2-c2/g_denom;
			double dgtheta_dcostheta=-2.*c2*(h-costheta)/(g_denom*g_denom);
			double fexp=exp(pow(lambda_3*(rij-rik),3.0));
			double dfexp_dr=3.*lambda_3*lambda_3*lambda_3*(rij-rik)*(rij-rik)*fexp;

			double dxj_costheta=(rik_hat[0]-costheta*rij_hat[0])/rij;
			double dyj_costheta=(rik_hat[1]-costheta*rij_hat[1])/rij;
			double dzj_costheta=(rik_hat[2]-costheta*rij_hat[2])/rij;
			double dxk_costheta=(rij_hat[0]-costheta*rik_hat[0])/rik;
			double dyk_costheta=(rij_hat[1]-costheta*rik_hat[1])/rik;
			double dzk_costheta=(rij_hat[2]-costheta*rik_hat[2])/rik;

			double dcostheta[]= {
					dxj_costheta+dxk_costheta,
					dyj_costheta+dyk_costheta,
					dzj_costheta+dzk_costheta
			};

			double eq_dfc=dfc_dr*gtheta*fexp;
			double eq_dgtheta=fc*dgtheta_dcostheta*fexp;
			double eq_dfexp=fc*gtheta*dfexp_dr;

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

		double bij, att;

		if(zeta !=0.0) {
			double b_kern=1+pow(beta*zeta,n);
			bij=pow(b_kern,-1.0/(2.*n));
			double db_dzeta=-0.5*pow(b_kern,-(1./(2.*n) + 1.0))*pow(beta,n)*pow(zeta,n-1.0);

			double dbij[]={
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
		double rep= a*exp(-lambda_1*rij);
		double dv_dr=lambda_1*rep-lambda_2*bij*att;
		double dvc_dr=CutFunction(rij)*dv_dr+dCutFunction(rij)*(rep-bij*att);

		Atoms(at).fx+=dvc_dr*rij_hat[0];
		Atoms(at).fy+=dvc_dr*rij_hat[1];
		Atoms(at).fz+=dvc_dr*rij_hat[2];
		Atoms(at).potential=Atoms(to).potential=CutFunction(rij)*(rep+bij*att);

	}

/*
	void ComputeTersoff(int at, int to, double &r,
			            double &dx, double &dy, double &dz,
			            double f_at[], double &pat, double f_to[], double &pto) {


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

	void Compute(int at, int to) {
		double dx, dy, dz;
		double r = sqrt(CalcSqrDistance(at, to, dx, dy, dz));

		if(r<CutRadius) {
			double pat,pto;
			double f_at[3], f_to[3];
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
