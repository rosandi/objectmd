/*========================================================*/
/*                                                        */
/*  "potential.h"                                         */
/*                                                        */
/*   part of an MD-simulation                             */
/*                                                        */
/*   (c) 01-00, Cemal Engin                               */
/*                                                        */
/*    deutschsprachiger Kommentar: R. Aderjan             */
/*                                                        */
/*========================================================*/


/*============================================================================*/

/*----------------------------------------*/
/*          EAM - Potential               */
/*----------------------------------------*/

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
#define Z1                      13
#define Z2                      13 /* Atomic number of Al */


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
#define Rmin   2   /* 2=int(x[]) wobei x[] kleinster wert in der tabelle */


/*============================================================================*/

/*============================================================================*/

double ZBLFactor_Al, ZBLFactor_AlXe, ZBLCoeff[4], ZBLExp_Al[4],ZBLExp_AlXe[4];
double knot[15],emb[13],b_emb[13],c_emb[13],d_emb[13];
double pot[15],b_pot[15],c_pot[15], d_pot[15];
double dent[15],b_dent[15],c_dent[15],d_dent[15];

/*============================================================================*/



void InitPotential (void)
{
  double ZBLNorm_Al,ZBLNorm_AlXe;
  
  int n;
  double a,b,c,d,e;
  FILE *finp1,*finp2,*finp3;


  /*----------------------------------------*/
  /* EAM stuff                              */
  /*----------------------------------------*/

  
  finp1=fopen("potential.dat", "r");
  finp2=fopen("atomden.dat","r");
  finp3=fopen("embeddfct.dat","r");

n=0;
  while(fscanf(finp1, " %lf %lf %lf %lf %lf \n ", &a, &b, &c, &d, &e) != EOF) {
   knot[n]=a; pot[n]=b; b_pot[n]=c; c_pot[n]=d; d_pot[n]=e;
   /* printf(" %lf %lf %lf %lf %lf \n ",a, b,c,d,e); */
   n++;
  };
fclose(finp1);


n=0;
  while(fscanf(finp2, " %lf %lf %lf %lf \n ", &a, &b, &c, &d ) != EOF) {
 dent[n]=a; b_dent[n]=b; c_dent[n]=c; d_dent[n]=d;
 /* printf(" %lf %lf %lf %lf \n ",a, b,c,d); */
  n++;
  };
fclose(finp2);


n=0;
  while(fscanf(finp3, " %lf %lf %lf %lf \n ", &a, &b, &c, &d ) != EOF) {
 emb[n]=a; b_emb[n]=b; c_emb[n]=c; d_emb[n]=d;
 /* printf(" %lf %lf %lf %lf \n ",a, b,c,d); */
 n++; 
};
fclose(finp3); 
  /*----------------------------------------*/
  /* ZBL stuff                              */
  /*----------------------------------------*/

  ZBLCoeff[0] = ZBL_COEFF0;
  ZBLCoeff[1] = ZBL_COEFF1;
  ZBLCoeff[2] = ZBL_COEFF2;
  ZBLCoeff[3] = ZBL_COEFF3;
 
 
  ZBLNorm_Al = 0.88534 * 0.5292 / ( exp(0.23*log(Z1)) + exp(0.23*log(Z1)) );
   ZBLFactor_Al = Z1*Z1*ELEMENTARY_CHARGE / 4./M_PI/FIELD_CONSTANT * 1e10;
 
  ZBLExp_Al[0] = ZBL_EXP0 / ZBLNorm_Al;
  ZBLExp_Al[1] = ZBL_EXP1 / ZBLNorm_Al;
  ZBLExp_Al[2] = ZBL_EXP2 / ZBLNorm_Al;
  ZBLExp_Al[3] = ZBL_EXP3 / ZBLNorm_Al;

}


void Embedd_energy (double Rho,double *Embedding, double *drEmbedding)
{ 
  double dRho;
  int pos;
 if (Rho >= 1.2999 ) {
    *Embedding= 0.712814*(Rho-1.3)-2.5515022;
    *drEmbedding = 0.712814;
  }
else {
  pos = (int)(10*Rho);
  dRho = Rho - .1*pos;
  *Embedding = emb[pos] + b_emb[pos]*dRho + c_emb[pos]*sq(dRho)+d_emb[pos]*cub(dRho);
  *drEmbedding = b_emb[pos] + 2*c_emb[pos]*dRho + 3*d_emb[pos]*sq(dRho);
  }
 
}


/* Ableitung der Elektronendichte     */
double CalcDrElecDens (double R) {

  double dx;
  int    pos,fertig;

  if ( R >= 2.021 ) {

    pos= ((int)R - Rmin)*tiefe;
   fertig = false;
    if ((R < knot[pos]) && (pos > 0)) { pos--; fertig = true;}
   if (R_CUT < R ) { pos=14;fertig = true;}

   while (!fertig && (pos <= 13))
    if (!((knot[pos] <= R) && (R < knot[pos+1]))) pos++; else fertig = true;

   dx = R - knot[pos];
   if (pos==14)
    return 0.;
   else 
     return b_dent[pos] + 2*c_dent[pos]*dx + 3*d_dent[pos]*sq(dx);
 
 } /* end if (R >=2.021) */

 else 
   if ( R > 0.9004871795 )  return 0.078; else return 0.;
   

}
   
double CalcElecDens (double R)
{  
   double dx;
   int    pos,fertig;

 	 if ( R >= 2.021 ) {

	   pos= ((int)R - Rmin)*tiefe;
	   fertig = false;

	   if ((R < knot[pos]) && (pos > 0)) { pos--; fertig = true;}
	   if (R_CUT < R ) { pos=14;fertig = true;}

	   while (!fertig && (pos <= 13))
	     if (!((knot[pos] <= R) && (R < knot[pos+1]))) pos++; else fertig = true;

	   dx = R - knot[pos];
	   if (pos==14)
              return  0.;
	    else
	      return dent[pos] + b_dent[pos]*dx + c_dent[pos]*sq(dx) + d_dent[pos]*cub(dx);
	 }

	 else 
	   if ( R > 0.9004871795 ) 
	     return 0.078*(R-2.021)+0.0874;
	   else 
	     return 0.;
}
 
/*================================================================================*/




void CalcPairContrib (double R, double *Phi, double *drPhi)

{

  int i,pos,fertig;
  double dx, ElDens, drElDens, Embed, drEmbed,Phi_ZBL,drPhi_ZBL;
  double Term, Factor, Pot, drPot,Phi_EAM,drPhi_EAM,T_INV;

  T_INV=1+exp((R-RM)/delta); /* 1/Fermifkt. */


/*------------------------------------------------*/
/*   EAM Potential                                */
/*------------------------------------------------*/
  if (R >= 2.021) {
      pos= ((int)R - Rmin)*tiefe;
      fertig = false;

    if ((R < knot[pos]) && (pos > 0)) { pos--; fertig = true;}
    if (R_CUT < R ) { pos = 14; fertig = true; }
     
   while (!fertig && (pos <= 13))
     if (!((knot[pos] <= R) && (R < knot[pos+1])))  pos++; else fertig = true;
   
   dx = R - knot[pos];
   if (pos==14) 
     {Phi_EAM=0. ; drPhi_EAM=0.; }
   else 
     {Phi_EAM = pot[pos]+b_pot[pos]*dx+c_pot[pos]*sq(dx)+d_pot[pos]*cub(dx);
      drPhi_EAM = b_pot[pos] + 2*c_pot[pos]*dx +3*d_pot[pos]*sq(dx);
     }

  }
  else { Phi_EAM = -6.9356*R + 15.9759476;drPhi_EAM=-6.9356;}


/*----------------------------------------------------*/
/* ZBL Potential                                      */
/*----------------------------------------------------*/

  Factor = ZBLFactor_Al / R;

    Pot = drPot = 0.;
    for (i = 0; i < 4; i++) {
     
      Pot   += (Term = ZBLCoeff[i] * exp(- ZBLExp_Al[i] * R));
      drPot += - Term * (ZBLExp_Al[i]+ 1. / R); 
      
    }
 
   
    drElDens = CalcDrElecDens (R);
    ElDens   = CalcElecDens (R);
    Embedd_energy ( ElDens ,&Embed,&drEmbed);

    Phi_ZBL   =  Factor * Pot - 2. * Embed;
    drPhi_ZBL = Factor * drPot - 2. * drEmbed * drElDens;



/*------------------------------------------------------------*/
/*  connecting ZBL and EAM mit Hilfe der Fermifunktion        */
/*------------------------------------------------------------*/

 drPhi_EAM  = drPhi_EAM*(1-1/T_INV) + Phi_EAM*(T_INV - 1)/delta/sq(T_INV);
 Phi_EAM    = Phi_EAM*(1-1/T_INV);

 drPhi_ZBL  = drPhi_ZBL/T_INV - Phi_ZBL*(T_INV-1)/delta/sq(T_INV);
 Phi_ZBL    = Phi_ZBL/T_INV;
    
 *Phi       = Phi_ZBL + Phi_EAM;
 *drPhi     = drPhi_ZBL + drPhi_EAM;


}   /*   CalcPairContrib   */


int main() {

    InitPotential();

    int ndat=500;
    double dr=R_CUT/(double)ndat;
    cout.setf(ios::scientific);
    cout.precision(15);

    //------------------- pair potential -------------------

    cout<< "##################################\n"
        << "# Force Matching Method, LEA\n#\n"
        << "# (c)Rosandi, 2007\n#\n"
        << "# BEGIN-PART PAIR\n"
        << "#$ NumberOfData " << ndat << '\n'
        << "#$ Spacing " << dr << '\n'
        << "#$ Offset " << dr << '\n';

    for (double r=dr; r<=R_CUT; r+=dr) {
		double phi, drphi;
		CalcPairContrib(r, &phi, &drphi);
		cout<<phi*r<<"\n";
    }
    cout << "# END-PART PAIR\n";

    //------------------ electron density --------------------

    cout << "# BEGIN-PART EDENS\n"
         << "#$ NumberOfData " << ndat << '\n'
         << "#$ Spacing " << dr << '\n'
         << "#$ Offset " << dr << '\n';

    double edmax=-1e6, ed;
    for (double r=dr; r<=R_CUT; r+=dr) {
        ed=CalcElecDens(r);
        if (edmax<ed) edmax=ed;
        cout << ed << '\n';
    }
    cout<< "# END-PART EDENS\n";

    //------------------- embedding function ------------------

    dr=edmax/(double)ndat;

    cout<< "# BEGIN-PART EMBED\n"
        << "#$ NumberOfData " << ndat << '\n'
        << "#$ Spacing " << dr << '\n'
        << "#$ Offset " << dr << '\n';

    for (double r=dr; r<=edmax; r+=dr) {
        double emb, demb;
        Embedd_energy(r, &emb, &demb);
        cout << emb << '\n';
    }
    cout<< "# END-PART EMBED\n";
    return 0;
}
