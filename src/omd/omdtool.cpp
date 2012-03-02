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
 * 
*/

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <cfloat>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <sstream>
#ifdef OMD_WITH_MUPARSER
#include <muParser.h>
#endif
#include <omd/omdtool.h>

using std::string;
using std::vector;
using std::ofstream;
using std::ostream;
using std::istringstream;

namespace omd {
    
  string msg;
  
  void die(string message) {
    msg=string("pid[")+as_string(getpid())+"] "+message;
    throw msg.c_str();
  }
  
  void die(string message, string val) {
    msg=string("pid[")+as_string(getpid())+"] "+message+": "+val;
    throw msg.c_str();
  }
  
  void mdassert(bool cond, string message) {
    if(!cond) die(message);
  }
  
  void mdassert(bool cond, string message, string val){
    if(!cond) die(message,val);
  }
  
  bool file_exist(string fname){
    struct stat filestat;
    if(stat(fname.c_str(), &filestat)) return false;
    return true;
  }
  
  string as_string(int val, const char* format){
    char st[256];
    if(format)sprintf(st,format,val);
    else sprintf(st,"%d",val);
    return string(st);
  }
  
  string as_string(double val, const char* format){
    char st[256];
    if(format)sprintf(st,format,val);
    else sprintf(st,"%0.5E",val);
    return string(st);
  }
  
  string as_string(void* val){
    char st[256];
    sprintf(st,"%lX",(long int)val);
    return string(st);
  }

  int as_int(string val) {
    int ii;
    istringstream ss(val);
    if(!(ss>>ii)) throw (string("conversion error: ")+val).c_str();
    return ii;
  }
  
  double as_double(string val) {
    double dd;
    istringstream ss(val);
    if(val=="inf") return DBL_MAX;
    if(val=="-inf") return -DBL_MAX;
    if(!(ss>>dd)) throw (string("conversion error: ")+val).c_str();
    return dd;
  }
  
#ifdef OMD_WITH_MUPARSER
  int as_int(string val, mu::Parser* parser) {
    int ii;
    if(val[0]=='{') {
      val=remove_char(val,"{ }");
      try {
        parser->SetExpr(val);
        ii=(int)parser->Eval();
      } catch (mu::Parser::exception_type &e){
        throw (string("parser error: "+e.GetMsg()+" eq:"+val).c_str());
      }
    } else ii=as_int(val);
    return ii;
  }

  double as_double(string val, mu::Parser* parser) {
    double dd;
    if(val[0]=='{') {
      val=remove_char(val,"{ }");
      try {
        parser->SetExpr(val);
        dd=parser->Eval();
      } catch (mu::Parser::exception_type &e){
        throw (string("parser error: "+e.GetMsg()+" eq:"+val).c_str());
      }
    } else dd=as_double(val);
    return dd;
  }
#endif
  
  string lower_case(string str){
    for(int i=0;i<(int)str.size();i++)str[i]=tolower(str[i]);
    return str;
  }
  
  string replace_char(string str, char cold, char cnew) {
    for(int i=0;i<(int)str.size();i++){
      if(str[i]==cold) str[i]=cnew;
    }
    return str;		
  }
  
  string replace_char(string str, string cold, char cnew) {
    for(int i=0;i<(int)str.size();i++){
      for(int c=0;c<(int)cold.size();c++) {
        if(str[i]==cold[c]) str[i]=cnew;
      }
    }
    return str;		
  } 
  
  string remove_char(string str, char delc) {
    string rets;
    for(int i=0;i<(int)str.size();i++) {
      if(str[i]!=delc)
        rets.append(1,str[i]);
    }
    return rets;    
  }
  
  string remove_char(string str, string delc) {
    for(int i=0;i<(int)delc.size();i++)
      str=remove_char(str,delc[i]);
    return str;
  }
  
  string trim(string str) {
    int i,j;
    for(i=0;i<(int)str.size();i++){
      if(str[i]==' ') continue;
      if(str[i]=='\t') continue;
      if(str[i]=='\r') continue;
      if(str[i]=='\n') continue;
      break;
    }
    
    for(j=(int)str.size()-1;j>=0;j--){
      if(str[j]==' ') continue;
      if(str[j]=='\t') continue;
      if(str[j]=='\r') continue;
      if(str[j]=='\n') continue;
      break;
    }
    return str.substr(i,j+1);
    
  }
  
  bool char_exist(string& str, char ch) {
    return (str.find(ch)!=str.npos);
  }
  
}
