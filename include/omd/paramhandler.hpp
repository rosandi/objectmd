/* omdlib
 ********************************************************
 *
 * ObjectMD
 * Molecular Dynamics Class Library
 * Version 3.0 (2009,2011)
 *
 * Released under:
 *       GNU GENERAL PUBLIC LICENSE Version 3
 *       Copyright (C) 2005  Yudi Rosandi
 *
 * Read license.txt file in the source's root directory.
 *
 ********************************************************
 * ObjectMD header file
*/

#ifndef _PARAM_HANDLER_HPP_
#define _PARAM_HANDLER_HPP_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctype.h>

#ifndef OMD_FLOAT
#define OMD_FLOAT double
#endif

using std::string;
using std::vector;

/**
 * @ingroup tool
 * @brief Parameter Handler
 * 
 * This class handles parameters, are stored in a string array.
 * An individual string is accessed by the index operator "[]". Parameters
 * can be threat as tag-value pair. A tag considers the next element as its
 * value. ParamHandler provides some conversion functions, returning the
 * value in different data types (int, double, string).
 * 
 * The class is able to take the parameters in the following ways:
 * 
 *   - from string array, like the parameter of C main function,
 *     main(argc, argv).
 *   - a single string contains one or more words separated by space.
 * 
 * Parameters may be assigned (using assign() function) or appended 
 * (append() function) into the existing parameters.
 * 
 * A list of parameters can be stored in a plain text file, and loaded using
 * read() function. Every word in the file is one parameter element, and indexed
 * sequentialy according to the reading position. A line begining with hash (#)
 * is cosidered as a comment and ignored. read() appends the parameters, not
 * assigning.
 * 
 * Using read_pseudo(), a set of special character mark (pseudo mark) can be put
 * in a text file, by default "#$". This kind of line is called the pseudo comment. 
 * ParamHandler will scan for lines begining with this mark, and load the text 
 * as parameter (excluding the mark). To stop scanning before eof, a double 
 * dash "--" must be written after the pseudo mark, separated by a space.
 * 
 * See an application, for example, in @link functable omd function tables. 
 * @endlink
 * 
*/

class ParamHandler {
	vector<string> spar;
	string prefix;
	string pseu_code;

public:
	// prefix: this is simply the parameter prefix like "-" or "--"
	// or anything

	ParamHandler() {
		prefix.assign("");
		pseu_code.assign("#$");
	}
	ParamHandler(string fpar){
		prefix.assign("");
		pseu_code.assign("#$");
		read(fpar);
	}
	ParamHandler(int argc, char* argv[]) {
		prefix.assign("");
		pseu_code.assign("#$");
		for(int i=0;i<argc;i++) {
			spar.push_back(argv[i]);
		}
	}
	
	int  size(){return spar.size();}
	void push_val(string p){spar.push_back(p);}
	void push_par(string p){spar.push_back(prefix+p);}
	void push_pair(string p, string v){push_par(p);push_val(v);}
	
	void clear(){spar.clear();}
	void set_prefix(const char* pp){prefix.assign(pp);}

	void copy(ParamHandler &p){
		prefix.assign(p.prefix);
		spar.clear();
		for(int i=0;i<(int)p.spar.size();i++){
			spar.push_back(p.spar[i]);
		}
	}

	void append(int argc, char* argv[]){
		for(int i=0;i<argc;i++) {
			spar.push_back(argv[i]);
		}
	}
	
	void append(string parstr) {
		std::istringstream s(parstr);
		while(s.good()){
			string ts;
			s >> ts;
			if(s.bad())break;
			if(ts!="")spar.push_back(ts);
		}
	}
	
	void assign(int argc, char* argv[]){clear(); append(argc, argv);}
	void assign(string parstr) {clear(); append(parstr);}

	// appending a file
	void read(string filename) {
		std::ifstream fl(filename.c_str());
		
		if(fl.bad()) {
			std::cerr<<"error reading parameter from file: "<<filename<<"\n";
			throw "error in param handler";
		}

		while(fl.good()) {
			char sln[1024];
			fl.getline(sln,1024);
			if(string(sln)=="--") break;
			std::istringstream s(sln);
			while(s.good()){
				string ts;s>>ts;
				if(s.bad())break;
				if(ts[0]!='#'&&ts!="")spar.push_back(ts);
				else break;
			}
		}
		fl.close();
	}
	
	/**
	 * Scans pseudo parameters in a file. Scanning stops at eof of when the end
	 * mark found. The end mark is a double dash "--" exactly after the pseudo 
	 * code "#$"
	 */
	
	void read_pseudo(string filename, string psmark="") {
		
		std::ifstream fl(filename.c_str());

		if(fl.bad()) {
			std::cerr<<"error reading parameter from file: "<<filename<<"\n";
			throw "error in param handler";
		}
		
		// pseudo format: #$psmark
		psmark=pseu_code+psmark;

		while(fl.good()) {
			char sln[1024];
			string pseu,ts;
			
			fl.getline(sln,1024);
			std::istringstream s(sln); s>>pseu>>ts;
			if(pseu!=psmark) continue;

			if(ts=="--")break;
			if(ts!="")spar.push_back(ts); // push it if its a valid value
			
			while(s>>ts){ // read the rest values
				if(s.bad())break;
				if(ts!="")spar.push_back(ts); else break;
			}
		}
		fl.close();
	}

	
	string& operator[](int idx){
		if(idx<(int)spar.size()) return spar[idx];
		throw "index out of range";	return spar[0]; // avoids warning
	}
	
	void shift(){
		if(spar.empty()) return;
		spar.erase(spar.begin());
	}
	
	void shift(int n){
		if(n>(int)spar.size())n=spar.size();
		for(int i=0;i<n;i++)shift();
	}
	
	bool exist(string p){
		for(int i=0;i<(int)spar.size();i++) {
			if(spar[i]==(prefix+p)) return true;
		}
		return false;
	}
	
	// returns the index of first occurrence
	int index_of(string p){
		for(int i=0;i<(int)spar.size();i++) {
			if(spar[i]==(prefix+p))return i;
		}
		throw (string("parameter ")+p+" does not exist").c_str();
	}

	vector<string>::iterator iter_of(string p){
		vector<string>::iterator pv;
		for(pv=spar.begin();pv!=spar.end();pv++)
			if(*pv==p) return pv;
		return spar.end();
	}

	void remove(string p) {
		vector<string>::iterator idx=iter_of(p);
		if(idx==spar.end()) return;
		spar.erase(idx);
	}

	void remove_pair(string p) {
		vector<string>::iterator idx=iter_of(p);
		if(idx==spar.end()) return;
		spar.erase(idx);
		if(idx!=spar.end())spar.erase(idx);
	}

	double double_value(int idx) {
		double retd=0.0;
		if(idx<(int)spar.size()){
			char st[256],*pt;
			memset(st,0,256);
			spar[idx].copy(st,256);
			retd=strtod(st,&pt);
			if(pt==st) throw "paramhandler: conversion to double failed";
		} else throw "paramhandler: conversion out of range";
		return retd;
	}
	
	int int_value(int idx){
		int reti=0;
		if(idx<(int)spar.size()){
			char st[256],*pt;
			memset(st,0,256);
			spar[idx+1].copy(st,256);
			reti=strtol(st,&pt,10);
			if(pt==st)throw "conversion to integer failed";
		} else throw "out of range";
		return reti;
	}
	
	string string_value(int idx){
		if(idx>=(int)spar.size()) throw "conversion: out of range";
		return spar[idx];
	}

    //--------------------------------------------------------------

	string string_value(string p, int index=0){
		if(!exist(p)) throw (string("parameter ")+p+" doesn't exist").c_str();
		int idx=index_of(p)+index;
		string rets("");
		if((idx+1)<(int)spar.size()) rets=spar[idx+1];
		if(rets=="") throw "no string value found";
		return rets;
	}
	
	// lower case
	string lower_string_value(string p){
		string rets=string_value(p);
		for(int i=0;i<(int)rets.size();i++) {
			rets[i]=tolower(rets[i]);
		}
		return rets;
	}
	    
	double double_value(string p, int index=0) {
		if(!exist(p)) throw (string("parameter ")+p+" doesn't exist").c_str();
		int idx=index_of(p)+index;
		
		if(idx<(int)spar.size()-1){
			char st[256],*pt;
			memset(st,0,256);
			spar[idx+1].copy(st,256);
			double retd=strtod(st,&pt);
			if(pt==st)
				throw (string("invalid double value for parameter ")+
				       spar[idx]+": "+st).c_str();
			return retd;
		}

		string msg=string("double value required for parameter ");

		if (index==0) msg.append(spar[idx]);
		else msg.append(spar[idx]+"(array)");
		throw (msg.c_str());

	}
	
	int int_value(string p, int index=0){
		if(!exist(p)) throw (string("parameter ")+p+" doesn't exist").c_str();
		int idx=index_of(p)+index;
		
		if(idx<(int)spar.size()-1){
			char st[256],*pt;
			memset(st,0,256);
			spar[idx+1].copy(st,256);
			int reti=strtol(st,&pt,10);
			if(pt==st)
				throw (string("invalid integer value for parameter ")+
				       spar[idx]+": "+st).c_str();
			return reti;
		}

		string msg=string("integer value required for parameter ");

		if (index==0) msg.append(spar[idx]);
		else msg.append(spar[idx]+"(array)");
		throw (msg.c_str());
	}

	// -- peek if exist, otherwise don't touch the variable
	bool peek(string p, string& var){
		if(exist(p)){var=string_value(p);return true;}
		return false;
	}

	bool peek(string p, double &var){
		if(exist(p)){var=double_value(p);return true;} 
		return false;
	}
	
	bool peek(string p, int &var){
		if(exist(p)){var=int_value(p);return true;}
		return false;
	}
	
	bool peek(string p, bool &var){
		if(exist(p)){
			string st=string_value(p);
			if(st=="true"||st=="yes") var=true;
			else var=false;
			return true;
		}
		return false;
	}
	
	// --- peek if exist, otherwise give default value
	
	bool peek(string p, string& var, string def){
		if(exist(p)){var=string_value(p);return true;}
		var=def; return false;
	}
	
	bool peek(string p, double &var, double def){
		if(exist(p)){var=double_value(p);return true;}
		var=def;return false;
	}
	
	bool peek(string p, int &var, int def){
		if(exist(p)){var=int_value(p);return true;}
		var=def;return false;
	}
	
	bool peek(string p, bool &var, bool def){
		if(exist(p)){
			string st=string_value(p);
			if(st=="true"||st=="yes") var=true;
			else var=false;
			return true;
		}
		var=def;return false;
	}

//----------------------------
	
	void dump(std::ostream& out){
		out << "number of parameters = " << spar.size() << "\n";
		for(int i=0;i<(int)spar.size();i++)
			out <<"["<<i<<"] "<<spar[i] << " ";
		out << "\n#total "<<spar.size()<<" parameters\n";
	}
	
	string raw_string() {
		string ss="";
		for(int i=0;i<(int)spar.size();i++) ss.append(spar[i]+" ");
		return ss;
	}
	
	// the parameter name and the end string is excluded 
	string raw_string(string p, string ends="--") {
		string ss="";
		for(int i=index_of(p)+1;i<(int)spar.size();i++) {
			if(spar[i]==ends) break;
			ss.append(spar[i]+" ");
		}
		return ss;
	}
	
	void set_pair(string p, string val) {
		if(exist(p)) {
			int idx=index_of(p);
			if(idx==(int)spar.size()-1) push_val(val);
			else spar[idx]=val;
		} else push_pair(p, val);			
	}
	
};

#endif
