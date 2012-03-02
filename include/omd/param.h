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

#ifndef _PARAM_HPP_
#define _PARAM_HPP_

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

#ifdef OMD_WITH_MUPARSER
#include <muParser.h>
#endif

using std::string;
using std::vector;

namespace omd {

class ParamHandler {
	vector<string> spar;
	string prefix;
	string pseu_code;
  
#ifdef OMD_WITH_MUPARSER
  mu::Parser* parser;
#else
  void* parser;
#endif
  
public:
	ParamHandler();
	ParamHandler(string fpar);
	ParamHandler(int argc, char* argv[]);
	int  size(){return spar.size();}
	void push_val(string p){spar.push_back(p);}
	void push_par(string p){spar.push_back(prefix+p);}
	void push_pair(string p, string v){push_par(p);push_val(v);}
	void clear(){spar.clear();}
	void set_prefix(const char* pp){prefix.assign(pp);}
	void copy(ParamHandler &p);
	void append(int argc, char* argv[]);
	void append(string parstr);
	void assign(int argc, char* argv[]) {clear(); append(argc, argv);}
	void assign(string parstr) {clear(); append(parstr);}
	void read(string filename);
	bool read_pseudo(string filename, string psmark="");
	string& operator[](int idx);
	void shift();
	void shift(int n);
	bool exist(string p);
	int index_of(string p);
	vector<string>::iterator iter_of(string p);
	void remove(string p);
	void remove_pair(string p);
	double double_value(int idx);
	int int_value(int idx);
	string string_value(int idx);
	string string_value(string p, int index=0);
	string lower_string_value(string p);
	double double_value(string p, int index=0);
	int int_value(string p, int index=0);
	bool peek(string p, string& var);
	bool peek(string p, double &var);
	bool peek(string p, int &var);
	bool peek(string p, bool &var);
	bool peek(string p, string& var, string def);
	bool peek(string p, double &var, double def);
	bool peek(string p, int &var, int def);
	bool peek(string p, bool &var, bool def);
	void dump(std::ostream& out);
	string raw_string();
	string raw_string(string p, string ends="--");
	void set_pair(string p, string val);
  
#ifdef OMD_WITH_MUPARSER
  void SetParser(mu::Parser* pars){parser=pars;}
#endif
  
};

}

#endif
