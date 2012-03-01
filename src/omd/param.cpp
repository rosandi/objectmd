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
#include <omd/omdtool.h>
#include <omd/param.h>

using namespace omd;

ParamHandler::ParamHandler() {
  prefix.assign("");
  pseu_code.assign("#$");
  parser==NULL;
}

ParamHandler::ParamHandler(string fpar){
  prefix.assign("");
  pseu_code.assign("#$");
  read(fpar);
}

ParamHandler::ParamHandler(int argc, char* argv[]) {
  prefix.assign("");
  pseu_code.assign("#$");
  for(int i=0;i<argc;i++) {
    spar.push_back(argv[i]);
  }
}

void ParamHandler::copy(ParamHandler &p){
  prefix.assign(p.prefix);
  spar.clear();
  for(int i=0;i<(int)p.spar.size();i++){
    spar.push_back(p.spar[i]);
  }
}

void ParamHandler::append(int argc, char* argv[]){
  for(int i=0;i<argc;i++) {
    spar.push_back(argv[i]);
  }
}

void ParamHandler::append(string parstr) {
  std::istringstream s(parstr);
  while(s.good()){
    string ts;
    s >> ts;
    if(s.bad())break;
    if(ts!="")spar.push_back(ts);
  }
}

void ParamHandler::read(string filename) {
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

bool ParamHandler::read_pseudo(string filename, string psmark) {
  
  std::ifstream fl(filename.c_str());
  
  if(fl.bad()) {
    std::cerr<<"error reading parameter from file: "<<filename<<"\n";
    throw "error in param handler";
  }
  
  // pseudo format: #$psmark
  psmark=pseu_code+psmark;
  int npse=0;
  
  while(fl.good()) {
    char sln[1024];
    string pseu,ts;
    
    fl.getline(sln,1024);
    std::istringstream s(sln); s>>pseu>>ts;
    if(pseu!=psmark) continue;
    npse++;
    
    if(ts=="--")break;
    if(ts!="")spar.push_back(ts); // push it if its a valid value
    
    while(s>>ts){ // read the rest values
      if(s.bad())break;
      if(ts!="")spar.push_back(ts); else break;
    }
  }
  fl.close();
  if(!npse) return false; // no pseudo command found
  return true;
}


string& ParamHandler::operator[](int idx){
  if(idx<(int)spar.size()) return spar[idx];
  throw "index out of range";	return spar[0]; // avoids warning
}

void ParamHandler::shift(){
  if(spar.empty()) return;
  spar.erase(spar.begin());
}

void ParamHandler::shift(int n){
  if(n>(int)spar.size())n=spar.size();
  for(int i=0;i<n;i++)shift();
}

bool ParamHandler::exist(string p){
  for(int i=0;i<(int)spar.size();i++) {
    if(spar[i]==(prefix+p)) return true;
  }
  return false;
}

// returns the index of first occurrence
int ParamHandler::index_of(string p){
  for(int i=0;i<(int)spar.size();i++) {
    if(spar[i]==(prefix+p))return i;
  }
  throw (string("parameter ")+p+" does not exist").c_str();
}

vector<string>::iterator ParamHandler::iter_of(string p){
  vector<string>::iterator pv;
  for(pv=spar.begin();pv!=spar.end();pv++)
    if(*pv==p) return pv;
  return spar.end();
}

void ParamHandler::remove(string p) {
  vector<string>::iterator idx=iter_of(p);
  if(idx==spar.end()) return;
  spar.erase(idx);
}

void ParamHandler::remove_pair(string p) {
  vector<string>::iterator idx=iter_of(p);
  if(idx==spar.end()) return;
  spar.erase(idx);
  if(idx!=spar.end())spar.erase(idx);
}

double ParamHandler::double_value(int idx) {
  if(idx>(int)spar.size()) throw "param: out of range";  
  return as_double(spar[idx],parser);
}

int ParamHandler::int_value(int idx){
  if(idx>=(int)spar.size()) throw "param: out of range";
  return as_int(spar[idx],parser);
}

string ParamHandler::string_value(int idx){
  if(idx>=(int)spar.size()) throw "conversion: out of range";
  return spar[idx];
}

//--------------------------------------------------------------

string ParamHandler::string_value(string p, int index){
  if(!exist(p)) throw (string("parameter ")+p+" doesn't exist").c_str();
  int idx=index_of(p)+index;
  string rets("");
  if((idx+1)<(int)spar.size()) rets=spar[idx+1];
  if(rets=="") throw "no string value found";
  return rets;
}

// lower case
string ParamHandler::lower_string_value(string p){
  string rets=string_value(p);
  for(int i=0;i<(int)rets.size();i++) {
    rets[i]=tolower(rets[i]);
  }
  return rets;
}

double ParamHandler::double_value(string p, int index) {
  if(!exist(p)) throw (string("parameter ")+p+" doesn't exist").c_str();
  int idx=index_of(p)+index+1;
  if(idx>=(int)spar.size()) throw "no parameter value found";
  return(as_double(spar[idx],parser));
}

int ParamHandler::int_value(string p, int index){
  if(!exist(p)) throw (string("parameter ")+p+" doesn't exist").c_str();
  int idx=index_of(p)+index+1;
  if(idx>=(int)spar.size()) throw "no parameter value found";
  return(as_int(spar[idx],parser));
}

// -- peek if exist, otherwise don't touch the variable
bool ParamHandler::peek(string p, string& var){
  if(exist(p)){var=string_value(p);return true;}
  return false;
}

bool ParamHandler::peek(string p, double &var){
  if(exist(p)){var=double_value(p);return true;} 
  return false;
}

bool ParamHandler::peek(string p, int &var){
  if(exist(p)){var=int_value(p);return true;}
  return false;
}

bool ParamHandler::peek(string p, bool &var){
  if(exist(p)){
    string st=string_value(p);
    if(st=="true"||st=="yes") var=true;
    else var=false;
    return true;
  }
  return false;
}

// --- peek if exist, otherwise give default value

bool ParamHandler::peek(string p, string& var, string def){
  if(exist(p)){var=string_value(p);return true;}
  var=def; return false;
}

bool ParamHandler::peek(string p, double &var, double def){
  if(exist(p)){var=double_value(p);return true;}
  var=def;return false;
}

bool ParamHandler::peek(string p, int &var, int def){
  if(exist(p)){var=int_value(p);return true;}
  var=def;return false;
}

bool ParamHandler::peek(string p, bool &var, bool def){
  if(exist(p)){
    string st=string_value(p);
    if(st=="true"||st=="yes") var=true;
    else var=false;
    return true;
  }
  var=def;return false;
}

//----------------------------

void ParamHandler::dump(std::ostream& out){
  out << "number of parameters = " << spar.size() << "\n";
  for(int i=0;i<(int)spar.size();i++)
    out <<"["<<i<<"] "<<spar[i] << " ";
  out << "\n#total "<<spar.size()<<" parameters\n";
}

string ParamHandler::raw_string() {
  string ss="";
  for(int i=0;i<(int)spar.size();i++) ss.append(spar[i]+" ");
  return ss;
}

// the parameter name and the end string is excluded 
string ParamHandler::raw_string(string p, string ends) {
  string ss="";
  for(int i=index_of(p)+1;i<(int)spar.size();i++) {
    if(spar[i]==ends) break;
    ss.append(spar[i]+" ");
  }
  return ss;
}

void ParamHandler::set_pair(string p, string val) {
  if(exist(p)) {
    int idx=index_of(p);
    if(idx==(int)spar.size()-1) push_val(val);
    else spar[idx]=val;
  } else push_pair(p, val);			
}
