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
#include <omd/tablereader.hpp>
#include <omd/paramhandler.hpp>

using std::ifstream;
using std::ofstream;
using std::string;
using std::istringstream;

// Calculate spline coeficients
void TableReader::calculate_coef() {
	OMD_FLOAT m;
	int    n=n_data-2;
	
	OMD_FLOAT *b=new OMD_FLOAT[n];
	OMD_FLOAT *d=new OMD_FLOAT[n];
	OMD_FLOAT *S=new OMD_FLOAT[n_data];
	
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
	set_name("TABLE READER");
	register_class(get_name());
	roffset = 0.0;
	coef=NULL;
	tablename="";
	if (table_filename!="") open(table_filename, table_name);
	format=plain;
	ready=allow_outrange_low=allow_outrange_hi=false;
	outrange_hi=outrange_low=0.0;
}

void TableReader::allocate() {
	coef = new CoefStruct[n_data];
	for (int i=0; i<n_data; i++) coef[i].a=coef[i].b=coef[i].c=coef[i].d=0.0;
}

// open file and read table data
void TableReader::open_omd(string table_filename, string table_name) {
	
	//---------initialization-------------			
	string sformat("PLAIN");

	filename.assign(table_filename);
	tablename.assign(table_name);
	
	param.read_pseudo(table_filename, table_name);
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
	else die("unsuported format "+sformat);

	ifstream fl(filename.c_str());
	assert(fl.good(), "error in reading data file", tablename);
	
	int lnum=0;
	char line[1024];
	bool partfound=false;
	string stag="#$"+tablename;
	
	while(fl.good()) {
		fl.getline(line, 1023); lnum++;
		istringstream sstr(line); string s, ss; sstr>>s>>ss;
		if(s==stag&&ss=="--") {partfound=true;break;}
	}
	assert(partfound, "can not find table "+tablename+"@"+filename);
	blog("table '"+tablename+"' found at line "+as_string(lnum), LOGCREATE);

	for(int i=0;i<n_data;i++){
		fl>>coef[i].d;
	}

	fl.close();	
	calculate_coef();
}

void TableReader::open(string table_filename, string tablename) {
	try{
		open_omd(table_filename, tablename);
	} catch(...){
		die("error reading potential table file: "+tablename+"@"+table_filename);
	}
	ready=true;
}

TableReader::~TableReader() {delete[] coef;}

void TableReader::read(OMD_FLOAT r, OMD_FLOAT& val, OMD_FLOAT& dval) {
	
	OMD_FLOAT rval=r-roffset;
	int    x = (int)(rval/dr);
	OMD_FLOAT dx = rval-(OMD_FLOAT)x*dr;
	
	if (x>n_data-1) {
		if(allow_outrange_hi) {val=dval=outrange_hi;return;}
		else {
			die("high limit exeeded in "+tablename+"@"+filename+" (read) "
				"index="+as_string(x)+" r="+as_string(r)+" delta="+as_string(dx)+
				" low_limit="+as_string(min_range())+" hi_limit="+as_string(max_range()));
		}
	}
	
	if (x<-1) {
		if(allow_outrange_low) {
			val=dval=outrange_low;
			return;
		} else {
			die("low limit exeeded in "+tablename+"@"+filename+" (read) "
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

OMD_FLOAT TableReader::read(OMD_FLOAT r) {
	
	OMD_FLOAT rval=r-roffset;
	int x = (int)(rval/dr);
	OMD_FLOAT dx = rval-(OMD_FLOAT)x*dr;
/*
	if(debug_flag) {
		std::cerr << "table="<<tablename<<"@"<<filename<<";"
		          << "r="<<r<<";"
		          << "offset="<<roffset<<";"
		          << "dr="<<dr<<";"
		          << "rval="<<rval<<";"
		          << "x="<<x<<";"
		          << "dx="<<dx<<";\n";
	}
*/
	if (x>n_data-1) {
		if(allow_outrange_hi) {
			return outrange_hi;
		} else {
			die("high limit exeeded in "+tablename+"@"+filename+" (read) "
				"index="+as_string(x)+" r="+as_string(r)+" delta="+as_string(dx)+
 				" low_limit="+as_string(min_range())+" hi_limit="+as_string(max_range()));

		}
	}
	
	if (x<(-1)) {
		if(allow_outrange_low) {
			return outrange_low;
		} else {
			std::cerr << "die here...!!!!!!!!!!!!! "<<tablename<<"@"<<filename<<"\n";
			die("low limit exeeded in "+tablename+"@"+filename+" (read) "
				"index="+as_string(x)+" r="+as_string(r)+" delta="+as_string(dx)+
				" low_limit="+as_string(min_range())+" hi_limit="+as_string(max_range()));

		}
	}

	if (x==n_data-1) {dx+=dr; x--;}
	if (x==-1) {x=0;}
	OMD_FLOAT val=(coef[x].d + dx*(coef[x].c + dx*(coef[x].b + dx*coef[x].a)));
	if(format==rmult) val/=r;

/*	
	if(debug_flag) {
		std::cerr<<tablename<<"@"<<filename<<" returned="<<val<<"\n";
	}
*/
	
	return val;
}

OMD_FLOAT TableReader::dread(OMD_FLOAT r) {
	OMD_FLOAT val, dval;
	read(r, val, dval);
	return dval;
}

OMD_FLOAT TableReader::dread2(OMD_FLOAT r) {
	OMD_FLOAT rval=r-roffset;
	int x = (int)(rval/dr);
	OMD_FLOAT dx = rval-(OMD_FLOAT)x*dr;
	

	if (x>n_data-1) {
		if(allow_outrange_hi) {
			return outrange_hi;
		} else {
			die("high limit exeeded in "+tablename+"@"+filename+" (dread2)"
				"index="+as_string(x)+" r="+as_string(r)+" delta="+as_string(dx)+
				" low_limit="+as_string(min_range())+" hi_limit="+as_string(max_range()));

		}
	}
	
	if (x<-1) {
		if(allow_outrange_low) {
			return outrange_low;
		} else {
			die("low limit exeeded in "+tablename+"@"+filename+" "
				"index="+as_string(x)+" r="+as_string(r)+" delta="+as_string(dx)+
				" low_limit="+as_string(min_range())+" hi_limit="+as_string(max_range()));
		}
	}

	if(x==n_data-1){dx+=dr; x--;}
	if(x==-1)x=0;

	OMD_FLOAT ddval=2.0*coef[x].b + 6.0*dx*coef[x].a;
	
	if(format==rmult) {
		OMD_FLOAT dval=dread(r);
		ddval=(ddval-2.0*dval)/r;
	}

	return ddval;
}

void TableReader::dump(string filename, int resolution) { 
	OMD_FLOAT rmax=roffset+(OMD_FLOAT)(n_data-1)*dr;
	OMD_FLOAT dx=rmax/(OMD_FLOAT)resolution;
	ofstream ofl(filename.c_str());
	assert(ofl.good(), "can not open file for writing", filename);
	ofl << "# Table of value and its first derivative\n";
	for (OMD_FLOAT x=roffset; x<=rmax; x+=dx) {
		ofl << x << ' ' << read(x) << ' ' << dread(x) << ' ' << dread2(x) << '\n';
	}
}

void TableReader::dump_var() {
	std::cout << "offset = " << roffset << ", " << "dr = " << dr << "\n";
}
