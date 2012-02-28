/* omdlib
 ********************************************************
 *name
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
 * Table reader:
 * reads potential file and do table splining
 *
*/

#include <fstream>
#include <sstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <cmath>
#include "omdtool.h"
#include "param.h"
#include "treader.h"

using namespace omd;
using std::ifstream;
using std::ofstream;
using std::string;
using std::istringstream;

// Calculate spline coeficients
void TableReader::calculate_coef() {
	double m;
	int    n=n_data-2;
	
	double *b=new double[n];
	double *d=new double[n];
	double *S=new double[n_data];
	
	for (int i=1;i<n+1;i++) d[i-1]=coef[i+1].d-2.0*coef[i].d + coef[i-1].d; 
	for (int i=1;i<n-1;i++) b[i]=4.0;
	b[0]=b[n-1]=5.0;	
	
	// Calculate matrix with Thomas algorithm
	for (int k=1;k<n;k++) {m=1.0/b[k-1]; b[k]-=m; d[k]-=m*d[k-1];}
	S[n+1]=S[n]=d[n-1]/b[n-1];
	for (int k=n-2; k>=0; k--) S[k+1]=(d[k]-1.0*S[k+2])/b[k];
	S[0]=S[1];
	
	for (int i=0; i<n_data-1;i++) {
		coef[i].a=(S[i+1]-S[i])/dr/dr/dr;
		coef[i].b= 3.0*S[i]/dr/dr;
		coef[i].c=(coef[i+1].d-coef[i].d-2.0*S[i]-S[i+1])/dr;
	}
	
	delete[] b;
	delete[] d;
	delete[] S;
}

TableReader::TableReader(string table_filename, string table_name) {
	roffset = 0.0;
	coef=NULL;
	name="";
	if (table_filename!="") open(table_filename, table_name);
	format=plain;
	ready=allow_outrange_low=allow_outrange_hi=false;
	outrange_hi=outrange_low=0.0;
}

void TableReader::allocate() {
	coef = new CoefStruct[n_data];
	for (int i=0; i<n_data; i++)
    coef[i].a=coef[i].b=coef[i].c=coef[i].d=0.0;
}

// open file and read table data
bool TableReader::open_omd(string table_filename, string table_name) {
	
	//---------initialization-------------			
	string sformat("PLAIN");

	filename.assign(table_filename);
	name.assign(table_name);
  
	if(!param.read_pseudo(table_filename, table_name)) return false;
  
	n_data=param.int_value("NumberOfData");
	allocate();
	dr=param.double_value("Spacing");
	param.peek("Offset", roffset, 0.0);
	param.peek("Format", sformat);	
	
	if(param.exist("OutOfRangeLower")){
		outrange_low=param.double_value("OutOfRangeLower");
		allow_outrange_low=true;
	}

	if(param.exist("OutOfRangeUpper")){
		outrange_hi=param.double_value("OutOfRangeUpper");
		allow_outrange_hi=true;
	}
	
	if(sformat=="PLAIN") format=plain;
	else if(sformat=="RMULT") format=rmult;
	else throw string("unsuported format ")+sformat;

	ifstream fl(filename.c_str());
	if(!fl.good()) throw "error reading data file";
	
	int lnum=0;
	char line[1024];
	bool partfound=false;
	string stag="#$"+name;
	
	while(fl.good()) {
		fl.getline(line, 1023); lnum++;
		istringstream sstr(line); string s, ss; sstr>>s>>ss;
		if(s==stag&&ss=="--") {partfound=true;break;}
	}
	if(!partfound) 
    throw string("can not find table ")+name+"@"+filename;
	
	for(int i=0;i<n_data;i++){
		fl>>coef[i].d;
	}

	fl.close();	
	calculate_coef();
  return true;
}

// raw table containes x y tabulated value pair
// with equispaced x.

bool TableReader::open_raw(string table_filename) {
  char line[2048];
  int lnum=0;
  double aold,xmin,delta=0.0;
  vector<double> vco;
  
  ifstream fl(filename.c_str());
  if(!fl.good()) throw "error reading data file";
  while(fl.good()) {
    double a,b;
    fl.getline(line,2047);
    istringstream ss(line);
    
    if(ss>>a>>b) {
      if(vco.empty()) xmin=a;
      else delta+=a-aold;
      vco.push_back(b);
      aold=a;
    }
    lnum++;
  }
  fl.close();
  if(vco.empty()) return false;
  
  // transfer
  name=table_filename.substr(0,table_filename.find('.'));
  n_data=vco.size();
  dr=delta/double(n_data-1);
  roffset=xmin;
  format=plain;
  allow_outrange_low=false;
  allow_outrange_hi=false;
  
  allocate();
  for(int i=0;i<n_data;i++) coef[i].d=vco[i];
  calculate_coef();
  
  return true;
}

void TableReader::open(string table_filename, string tablename) {
  
  if(!open_omd(table_filename, tablename))
    if(!open_raw(table_filename)) throw "failed loading table";
  
	ready=true;
}

TableReader::~TableReader() {
  delete[] coef;
}

void TableReader::read(double r, double& val, double& dval) {
	
	double rval=r-roffset;
	int    x = (int)(rval/dr);
	double dx = rval-(double)x*dr;
	
	if (x>n_data-1) {
		if(allow_outrange_hi) {val=dval=outrange_hi;return;}
		else {
			die("high limit exeeded in "+name+"@"+filename+" (read) "
				"index="+as_string(x)+" r="+as_string(r)+" delta="+as_string(dx)+
				" low_limit="+as_string(min_range())+" hi_limit="+as_string(max_range()));
		}
	}
	
	if (x<-1) {
		if(allow_outrange_low) {
			val=dval=outrange_low;
			return;
		} else {
			die("low limit exeeded in "+name+"@"+filename+" (read) "
				"index="+as_string(x)+" r="+as_string(r)+" delta="+as_string(dx)+
				" low_limit="+as_string(min_range())+" hi_limit="+as_string(max_range()));
				
		}
	}
	
	if (x==n_data-1) {dx+=dr; x--;}
	if (x==-1) {x=0;}
	val=(coef[x].d + dx*(coef[x].c + dx*(coef[x].b + dx*coef[x].a)));
	dval=(coef[x].c + dx*(2.0*coef[x].b + 3.0*dx*coef[x].a));
	
	if(format==rmult){
		val/=r;
		dval=(dval-val)/r;
	}
}

double TableReader::read(double r) {
	
	double rval=r-roffset;
	int x = (int)(rval/dr);
	double dx = rval-(double)x*dr;

	if (x>n_data-1) {
		if(allow_outrange_hi) {
			return outrange_hi;
		} else {
			die("high limit exeeded in "+name+"@"+filename+" (read) "
				"index="+as_string(x)+" r="+as_string(r)+" delta="+as_string(dx)+
 				" low_limit="+as_string(min_range())+" hi_limit="+as_string(max_range()));

		}
	}
	
	if (x<(-1)) {
		if(allow_outrange_low) {
			return outrange_low;
		} else {
			die("low limit exeeded in "+name+"@"+filename+" (read) "
				"index="+as_string(x)+" r="+as_string(r)+" delta="+as_string(dx)+
				" low_limit="+as_string(min_range())+" hi_limit="+as_string(max_range()));
		}
	}

	if (x==n_data-1) {dx+=dr; x--;}
	if (x==-1) {x=0;}
	double val=(coef[x].d + dx*(coef[x].c + dx*(coef[x].b + dx*coef[x].a)));
	if(format==rmult) val/=r;
	
	return val;
}

double TableReader::dread(double r) {
	double val, dval;
	read(r, val, dval);
	return dval;
}

double TableReader::dread2(double r) {
	double rval=r-roffset;
	int x = (int)(rval/dr);
	double dx = rval-(double)x*dr;
	

	if (x>n_data-1) {
		if(allow_outrange_hi) {
			return outrange_hi;
		} else {
			die("high limit exeeded in "+name+"@"+filename+" (dread2)"
				"index="+as_string(x)+" r="+as_string(r)+" delta="+as_string(dx)+
				" low_limit="+as_string(min_range())+" hi_limit="+as_string(max_range()));

		}
	}
	
	if (x<-1) {
		if(allow_outrange_low) {
			return outrange_low;
		} else {
			die("low limit exeeded in "+name+"@"+filename+" "
				"index="+as_string(x)+" r="+as_string(r)+" delta="+as_string(dx)+
				" low_limit="+as_string(min_range())+" hi_limit="+as_string(max_range()));
		}
	}

	if(x==n_data-1){dx+=dr; x--;}
	if(x==-1)x=0;

	double ddval=2.0*coef[x].b + 6.0*dx*coef[x].a;
	
	if(format==rmult) {
		double dval=dread(r);
		ddval=(ddval-2.0*dval)/r;
	}

	return ddval;
}

void TableReader::dump(std::ostream& ofl, int resolution) { 
	double rmax=roffset+(double)(n_data-1)*dr;
	double dx=rmax/(double)resolution;
  std::cerr<<" max="<<rmax<<" dr="<<dr<<" dx="<<dx<<std::endl;
	ofl << "# Table of value first and second derivative\n";
	for (double x=roffset; x<=rmax; x+=dx) {
		ofl << x << ' ' << read(x) << ' ' << dread(x) << ' ' << dread2(x) << '\n';
	}
}

void TableReader::dump(string filename, int resolution) { 
	ofstream ofl(filename.c_str());
  if(!ofl.good())
    throw string("can not open file for writing ")+filename;
  dump(ofl,resolution);
  ofl.close();
}

void TableReader::dump_var() {
	std::cout << "offset = " << roffset << ", " << "dr = " << dr << "\n";
}
