/*
 ********************************************************
 * 
 *   (c) 01-00, Cemal Engin                               
 *          09, Yudi Rosandi
 *
 *   source: Interatomic potentials from first-principles
 *           calculations: the force-matching method
 * 
 *           F. Ercolessi and J. B. Adams
 *           Europhys. Lett, 26 (8), pp. 583-588 (1994)
 * 
 *           The potential is splined to ZBL at high energy region.
 *
 * 
*/

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <omd/conversion.hpp>
#include <omd/paramhandler.hpp>

using std::istringstream;
using std::ifstream;
using std::string;
using std::ios;
using std::cout;
using std::cerr;


#define POTENTIAL_NAME          "Aluminium, EAM-Potential"

#define REL_MASS_Al             26.982
#define LATTICE_CONSTANT         4.032
#define R_CUT                    5.5580
#define Z1                      13.0
#define Z2                      13.0 /* Atomic number of Al */
#define PI M_PI

/*----------------------------------------*/
/* general stuff                          */
/*----------------------------------------*/

#define ZBL_COEFF0               0.1818
#define ZBL_COEFF1               0.5099
#define ZBL_COEFF2               0.2802
#define ZBL_COEFF3               0.02817
#define ZBL_EXP0                 3.2
#define ZBL_EXP1                 0.9423
#define ZBL_EXP2                 0.4029
#define ZBL_EXP3                 0.2016

/*----------------------------------------*/
/* FermiFunktion Parameter                */
/*----------------------------------------*/

#define RM                      1.75
#define delta                   0.1

/*----------------------------------------*/


#define ELEMENTARY_CHARGE        1.6022e-19
#define FIELD_CONSTANT           8.8542e-12

#define true   1
#define false  0
#define sq(x) ((x)*(x))
#define cub(x) ((sq(x))*(x))
#define tiefe  4
#define Rmin   2

double ZBLFactor_Al, ZBLFactor_AlXe, ZBLCoeff[4], ZBLExp_Al[4],ZBLExp_AlXe[4];
double knot[15],a_emb[13],b_emb[13],c_emb[13],d_emb[13];
double pot[15],b_pot[15],c_pot[15], d_pot[15];
double dent[15],b_dent[15],c_dent[15],d_dent[15];

double Embed_energy(double Rho){
	double dRho, Embe;
	int pos;
	if (Rho >= 1.2999 ) {
		Embe= 0.712814*(Rho-1.3)-2.5515022;
	} else {
		pos = (int)(10.0*Rho);
		dRho = Rho - .1*pos;
		Embe = a_emb[pos] + b_emb[pos]*dRho + c_emb[pos]*sq(dRho)+d_emb[pos]*cub(dRho);
	}
	return Embe;
}

double dEmbed_energy (double Rho) { 
	double dRho, dEmbe;
	int pos;
	if (Rho >= 1.2999 ) {
		dEmbe = 0.712814;
	} else {
		pos = (int)(10*Rho);
		dRho = Rho - .1*pos;
		dEmbe = b_emb[pos] + 2*c_emb[pos]*dRho + 3*d_emb[pos]*sq(dRho);
	}
	return dEmbe;
}

/* Ableitung der Elektronendichte     */
double CalcDrElecDens (double R) {
	double dx;
	int pos,fertig;
	if ( R >= 2.021 ) {
		pos= ((int)R - Rmin)*tiefe;
		fertig = false;
		if ((R < knot[pos]) && (pos > 0)) { pos--; fertig = true;}
		if (R_CUT < R ) { pos=14;fertig = true;}
		while (!fertig && (pos <= 13))
		if (!((knot[pos] <= R) && (R < knot[pos+1]))) pos++; else fertig = true;
		dx = R - knot[pos];
		if (pos==14) return 0.;
		return b_dent[pos] + 2*c_dent[pos]*dx + 3*d_dent[pos]*sq(dx);
	} /* end if (R >=2.021) */
	if( R > 0.9004871795 ) return 0.078;
	return 0.;
}

double CalcElecDens (double R){  
	double dx;
	int    pos,fertig;
	if ( R >= 2.021 ) {
		pos=((int)R - Rmin)*tiefe;
		fertig = false;
		if ((R < knot[pos]) && (pos > 0)) { pos--; fertig = true;}
		if (R_CUT < R ) { pos=14;fertig = true;}
		while (!fertig && (pos <= 13))
	 		if (!((knot[pos] <= R) && (R < knot[pos+1]))) pos++; else fertig = true;		
		dx = R - knot[pos];
		if (pos==14) return  0.0;
		return dent[pos]+b_dent[pos]*dx+c_dent[pos]*sq(dx)+d_dent[pos]*cub(dx);
	}
	if(R>0.9004871795)return 0.078*(R-2.021)+0.0874; 
	return 0.;
}

void CalcPairContrib (double R, double& Phi, double& drPhi) {
	int i,pos,fertig;
	double dx, ElDens, drElDens, Embed, drEmbed,Phi_ZBL,drPhi_ZBL;
	double Term, Factor, Pot, drPot,Phi_EAM,drPhi_EAM,T_INV;
	T_INV=1+exp((R-RM)/delta); /* 1/Fermifkt. */
	if (R >= 2.021) {
		pos= ((int)R - Rmin)*tiefe;
		fertig = false;
		if ((R < knot[pos])&& (pos > 0)) { pos--; fertig = true;}
		if (R_CUT < R ) { pos = 14; fertig = true; }
		while (!fertig && (pos <= 13)) {
			if (!((knot[pos] <= R)&&(R < knot[pos+1])))pos++;
			else fertig = true;
		}
		dx = R - knot[pos];
		if(pos==14){Phi_EAM=0. ; drPhi_EAM=0.; }
		else {
			Phi_EAM = pot[pos]+b_pot[pos]*dx+c_pot[pos]*sq(dx)+d_pot[pos]*cub(dx);
			drPhi_EAM = b_pot[pos] + 2*c_pot[pos]*dx +3*d_pot[pos]*sq(dx);
		}
	} else { Phi_EAM = -6.9356*R + 15.9759476;drPhi_EAM=-6.9356;}
	Factor = ZBLFactor_Al / R;
	Pot = drPot = 0.;
	for (i = 0; i < 4; i++) {
		Pot   += (Term = ZBLCoeff[i] * exp(- ZBLExp_Al[i] * R));
		drPot += - Term * (ZBLExp_Al[i]+ 1. / R); 
	}
	drElDens = CalcDrElecDens (R);
	ElDens   = CalcElecDens (R);
	
	// CHECK THIS!
	Embed=Embed_energy( ElDens);
	drEmbed=dEmbed_energy(ElDens);
	
	Phi_ZBL    =  Factor * Pot - 2. * Embed;
	drPhi_ZBL  = Factor * drPot - 2. * drEmbed * drElDens;
	drPhi_EAM  = drPhi_EAM*(1-1/T_INV) + Phi_EAM*(T_INV - 1)/delta/sq(T_INV);
	Phi_EAM    = Phi_EAM*(1-1/T_INV);
	drPhi_ZBL  = drPhi_ZBL/T_INV - Phi_ZBL*(T_INV-1)/delta/sq(T_INV);
	Phi_ZBL    = Phi_ZBL/T_INV;
	Phi        = Phi_ZBL + Phi_EAM;
	drPhi      = drPhi_ZBL + drPhi_EAM;
}

void assert(bool cond, string msg){
	if(!cond) {
		cerr<<msg<<"\n";
		exit(1);
	}
}

void InitPotential(string tablefile){
	double ZBLNorm_Al,ZBLNorm_AlXe;

	double a,b,c,d,e;
	
	ifstream finp(tablefile.c_str());
	string s;

	char sln[256]="";	

	while(s!="#?table:pair_potential") {
		finp.getline(sln, 256);
		istringstream ss(sln);ss>>s;
		assert(finp.good(), "bad table format: "+tablefile);
	}

	for(int i=0;i<15;i++){
		finp >>a>>b>>c>>d>>e;
		assert(finp.good(), "error reading pair potential data from "+tablefile);
		knot[i]=a; pot[i]=b; b_pot[i]=c; c_pot[i]=d; d_pot[i]=e;
	}
	
	finp.seekg(0);
	while(s!="#?table:electron_density") {
		finp.getline(sln, 256);
		istringstream ss(sln);ss>>s;
		assert(finp.good(), "bad table format: "+tablefile);
	}

	for(int i=0;i<15;i++) {
		finp >>a>>b>>c>>d;
		assert(finp.good(), "error reading electron density data from "+tablefile);
		dent[i]=a; b_dent[i]=b; c_dent[i]=c; d_dent[i]=d;
	}
	
	finp.seekg(0);
	while(s!="#?table:embedding_function") {
		finp.getline(sln, 256);
		istringstream ss(sln);ss>>s;
		assert(finp.good(), "bad table format: "+tablefile);
	}

	for(int i=0;i<13;i++) {
		finp >>a>>b>>c>>d;
		assert(finp.good(), "error reading electron density data from "+tablefile);
		a_emb[i]=a; b_emb[i]=b; c_emb[i]=c; d_emb[i]=d;
	}

	finp.close(); 

	ZBLCoeff[0] = ZBL_COEFF0;
	ZBLCoeff[1] = ZBL_COEFF1;
	ZBLCoeff[2] = ZBL_COEFF2;
	ZBLCoeff[3] = ZBL_COEFF3;
	ZBLNorm_Al = 0.88534 * 0.5292 / (pow(Z1,0.23) + pow(Z1,0.23));
	ZBLFactor_Al = Z1*Z1*ELEMENTARY_CHARGE / 4./PI/FIELD_CONSTANT * 1e10;
	ZBLExp_Al[0] = ZBL_EXP0 / ZBLNorm_Al;
	ZBLExp_Al[1] = ZBL_EXP1 / ZBLNorm_Al;
	ZBLExp_Al[2] = ZBL_EXP2 / ZBLNorm_Al;
	ZBLExp_Al[3] = ZBL_EXP3 / ZBLNorm_Al;
}

int main(int argc, char* argv[]) {

    if(argc==1){
        cerr << "syntax: "<<argv[0]<<" parameter_file\n";
        exit(1);
    }

    InitPotential(argv[1]);

    int ndat=1000;
    double dr=R_CUT/(double)ndat;
    cout.setf(ios::scientific);
    cout.precision(15);

    //------------------- pair potential -------------------

    cout<< "#################################################################\n"
        << "# Force Matching Method\n#"
        << "#     Ercolessi, F. & Adams, J. B., Europhys. Lett., 1994, 26, 583\n"
        << "# Fitted to ZBL at high energy:\n"
        << "#     Cemal Engin, 2004\n#\n"
        << "# (c) Rosandi, 2007\n#\n"
        << "#$PAIR NumberOfData " << ndat << '\n'
        << "#$PAIR Spacing " << dr << '\n'
        << "#$PAIR Offset " << dr << '\n'
        << "#$PAIR Format RMULT\n"
        << "#$PAIR --\n";

    for (double r=dr; r<=R_CUT; r+=dr) {
		double phi, drphi;
		CalcPairContrib(r, phi, drphi);
		cout<<phi*r<<"\n";
    }

    //------------------ electron density --------------------

    cout << "#$EDENS NumberOfData " << ndat << '\n'
         << "#$EDENS Spacing " << dr << '\n'
         << "#$EDENS Offset " << dr << '\n'
         << "#$EDENS --\n";

    double edmax=-1e6, ed;
    for (double r=dr; r<=R_CUT; r+=dr) {
        ed=CalcElecDens(r);
        if (edmax<ed) edmax=ed;
        cout << ed << '\n';
    }

    //------------------- embedding function ------------------
    // 20 times edmax
    dr=20.0*edmax/(double)ndat;

    cout<< "#$EMBED NumberOfData " << ndat << '\n'
        << "#$EMBED Spacing " << dr << '\n'
        << "#$EMBED Offset 0.0\n"
        << "#$EMBED --\n";

    for (double r=0.0; r<20.0*edmax; r+=dr) {cout << Embed_energy(r) << '\n';}
    return 0;
}
