/* 
 *****************************************************
 * Object-MD
 * Object Oriented Molecular Dynamics Program
 * 
 * (c)Yudi Rosandi, 2005
 *
 * EAM potential generator
 * 
 * High energy: ZBL
 * Low energy: tight binding
 * Cut radius: splined to zero with Tersoff function
 *
 *****************************************************
*/

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <omd/conversion.hpp>
#include <omd/paramhandler.hpp>

using std::ifstream;
using std::string;
using std::ios;
using std::cout;
using std::cerr;

double R_CUT;
double PHI0;
double ALPHA;
double G0;
double BETA;
double ZC0, ZC1, ZC2, ZC3;
double ZEXP0, ZEXP1, ZEXP2, ZEXP3;
double SPLINE0, SPLINE1;
double P0, P1, P2, P3, P4, P5;
double TERSOFF_DELTA;

string ELEMENT;
bool setfl=false;

#define Tersoff(r)  (0.5 - 0.5 * sin(TersoffFreq * r + TersoffPhase))

double ElecDens, TersoffRadius, TersoffFreq, TersoffAmpl, TersoffPhase;
double ZZBLFactor;
double zze0, zze1, zze2, zze3;

void read_parameter(const char* pf) {
try {
	param_handler p; p.read(pf);
//	p.dump(cerr);
	ELEMENT=p.string_value("ELEMENT");
	R_CUT=p.double_value("R_CUT");
	PHI0=p.double_value("PHI0");
	ALPHA=p.double_value("ALPHA");
	G0=p.double_value("G0");     
	BETA=p.double_value("BETA");
	ZC0=p.double_value("ZBL_C0");
	ZC1=p.double_value("ZBL_C1");    
	ZC2=p.double_value("ZBL_C2");
	ZC3=p.double_value("ZBL_C3");
	ZEXP0=p.double_value("ZBL_EXP0");
	ZEXP1=p.double_value("ZBL_EXP1");
	ZEXP2=p.double_value("ZBL_EXP2");
	ZEXP3=p.double_value("ZBL_EXP3");
	SPLINE0=p.double_value("SPLINE0");
	SPLINE1=p.double_value("SPLINE1");
	P0=p.double_value("P0");
	P1=p.double_value("P1");
	P2=p.double_value("P2");
	P3=p.double_value("P3");
	P4=p.double_value("P4");
	P5=p.double_value("P5");

	TERSOFF_DELTA=p.double_value("TERSOFF_DELTA");
	TersoffRadius =   R_CUT - 2.* TERSOFF_DELTA;
	TersoffFreq   =   M_PI/2./TERSOFF_DELTA;
	TersoffAmpl   = - 0.5 * TersoffFreq;
	TersoffPhase  =   M_PI*(TERSOFF_DELTA-R_CUT)/2./TERSOFF_DELTA;
	
	if(p.exist("FORMAT"))
		if(p.lower_string_value("FORMAT")=="setfl"){
			setfl=true;
			cerr<<"using setfl format\n";
	}

} catch(const char* serr) {
	cerr<<"error: "<<serr<<std::endl;
	exit(1);
} catch(...){
	cerr<<"unknown error\n";
	exit(1);
}
}

void Init()
{
	double ZBLNorm;
	double Z1=78.0, Z2=78.0;	
	
	// This is in electron volts!! loss 1 e
	ZZBLFactor=1.0e10*Z1*Z2*E_CHARGE/(4.*M_PI*F_CONST);
	ZBLNorm=0.88534*0.5292/(pow(Z1,0.23)+pow(Z2, 0.23));

	zze0=ZEXP0/ZBLNorm;
	zze1=ZEXP1/ZBLNorm;
	zze2=ZEXP2/ZBLNorm;
	zze3=ZEXP3/ZBLNorm;
	
	TersoffRadius=R_CUT-2.*TERSOFF_DELTA;
	TersoffFreq=M_PI/(2.*TERSOFF_DELTA);
	TersoffAmpl=-0.5*TersoffFreq;
	TersoffPhase= M_PI*(TERSOFF_DELTA-R_CUT)/(2.*TERSOFF_DELTA);
}

double CalcZBL(double R)
{
    double pot= ZC0*exp(-zze0*R)+ZC1*exp(-zze1*R)+ZC2*exp(-zze2*R)+ZC3*exp(-zze3*R);
	return (pot*ZZBLFactor/R + 2.*sqrt(G0*exp(-BETA*R)));
}

double CalcTB(double R)
{
    double pot = PHI0*exp(-ALPHA*R);
   	return (R<TersoffRadius)?pot:pot*Tersoff(R);
}

double CalcSpline(double R)
{
   	double r=R-SPLINE0;
   	return P0+r*(P1+r*(P2+r*(P3+r*(P4+r*P5))));
}

double Density(double R)
{return (R<TersoffRadius)?(G0*exp(-BETA*R)):(G0*exp(-BETA*R)*Tersoff(R));}

double Embed(double rho)
{return -sqrt(rho);}


int main(int argc, char* argv[]) {

	if(argc==1){
		cerr << "syntax: "<<argv[0]<<" parameter_file\n";
		exit(1);
	}
	
	read_parameter(argv[1]);
		
	int ndat=5000;
	double dr=R_CUT/(double)ndat;
	cout.setf(ios::scientific);
	cout.precision(15);
	
	Init();	
	
	//------------------- pair potential -------------------
		
	cout<< "##################################\n"
        << "# Tight-binding eam potential\n"
        << "# Thomas J. Colla style\n#\n"
        << "# (c)Rosandi, 2007\n#\n"
        << "# Element = "<<ELEMENT<<"\n#\n"
        << "# BEGIN-PART PAIR\n"
	 	<< "#$ NumberOfData " << ndat << '\n'
		<< "#$ Spacing " << dr << '\n'
		<< "#$ Offset " << dr << '\n';
	if(setfl) cout<<"#$ Format SETFL\n";
	
	for (double r=dr; r<=R_CUT; r+=dr) {
		if (r<SPLINE0)         {cout << CalcZBL(r)*(setfl?r:1.0) << '\n';}
		else if (r<SPLINE1)    {cout << CalcSpline(r)*(setfl?r:1.0) << '\n';}
		else if (r<TersoffRadius){cout << CalcTB(r)*(setfl?r:1.0) << '\n';}
		else {cout<< CalcTB(r)*Tersoff(r)*(setfl?r:1.0) << '\n';}
	}

	cout << "# END-PART PAIR\n";
	
	//------------------ electron density --------------------
	
	cout << "# BEGIN-PART EDENS\n"
		 << "#$ NumberOfData " << ndat << '\n'
		 << "#$ Spacing " << dr << '\n'
		 << "#$ Offset " << dr << '\n';
	if(setfl) cout<<"#$ Format SETFL\n";
	
	double edmax=-1000, ed;
	for (double r=dr; r<=R_CUT; r+=dr) {
		ed=Density(r);
		if (edmax<ed) edmax=ed;
		cout << ed*(setfl?r:1.0) << '\n';
	}
	cout<< "# END-PART EDENS\n";
	
	//------------------- embedding function ------------------
	
	dr=edmax/(double)ndat;

	cout<< "# BEGIN-PART EMBED\n"
		<< "#$ NumberOfData " << ndat << '\n'
		<< "#$ Spacing " << dr << '\n'
		<< "#$ Offset " << dr << '\n';
	if(setfl) cout<<"#$ Format SETFL\n";

	for (double r=dr; r<=edmax; r+=dr) {cout << Embed(r)*(setfl?r:1.0) << '\n';}
	cout<< "# END-PART EMBED\n";
	return 0;
}
