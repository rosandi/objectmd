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

using std::string;
using std::vector;
using std::ofstream;
using std::ostream;

#ifndef _MD_TOOLKIT_H_
#define _MD_TOOLKIT_H_

#define LOGMEMORY  1
#define LOGCREATE  2
#define LOGDESTROY 4
#define LOGINFO    8
#define LOGWARNING 16

#ifndef _OMD_TYPES_
#define _OMD_TYPES_

#define OMD_FLOAT double

#endif

// these two files must be defined and opened else where....
// OMD does it in toolkit.cpp

extern ofstream blogstream;
extern ofstream memlogstream;

/**
 * @ingroup tool
 * @brief The collection of tool functions
 * 
 * This class is the root of all Object-MD classes. It provides usefull
 * functionalities, used commonly by all objects.
 * 
 * The variables and functions can be classified into:
 * 
 *   - class type and identity management: type registerring, checking, naming.
 *   - logging: manage the run-time log and memory allocation log files.
 *   - tools: checking file, path. 
 *   - data conversion and string manipulation.
 *   - proccess management: warning, assertion, and terminating the running
 *     program.
 * 
*/

class MDToolkit {

	string object_name;
	int logflags;
	vector<string> class_type;
	
	
protected:
	MDToolkit* logger;
	bool debug_flag;

    // class name will be registered in lower case
    void register_class(string type) {
    	assert(!type_of(type), "not allowed: registering the same type as parent", type);
    	type=lower_case(type);
    	type=replace_char(type, ' ', '_');
    	class_type.push_back(type);
	}

public:

	MDToolkit(){
		// default log flag is defind in config.hpp
		logflags=LOGFLAG;
		object_name.assign("md_toolkit");
		register_class("md_toolkit");
		debug_flag=false;
		logger=this;
	}
	
	virtual ~MDToolkit(){}

	string get_name(){return object_name;}
	void set_name(string newname){object_name.assign(newname);}	
	void enable_log(int flags=0xFF){logflags|=flags;}
	void disable_log(int flags=0xFF){logflags&=(~flags);}
	void log_flagset(int flags){logflags=flags;}
	int log_enabled(int flags=0xFF){return logflags&flags;}

	void debug_enable() {debug_flag=true;}
	void debug_disable() {debug_flag=false;}

	MDToolkit* get_logger(){return logger;}
	virtual void set_logger(MDToolkit* mdlog){logger=mdlog;}

	virtual void sendlog(string slog){blogstream<<slog<<std::endl;}
	virtual void sendmemlog(string slog){memlogstream<<slog<<std::endl;}

	virtual void blog(string logtext, int flags=LOGINFO){
		if(logflags&flags){
			logger->sendlog("("+get_name()+") "+logtext);
		}
	}

	virtual void blog(string logprefix, string logtext, int flags=LOGINFO){
		if(logflags&flags){
			logger->sendlog("("+get_name()+") "+logprefix+" "+logtext);
		}
	}
	
	virtual void logmem(void* ptr){ // free
		if(logflags&LOGMEMORY)
			logger->sendmemlog("free ("+get_name()+") "+as_string(ptr));
	}

	virtual void logmem(void* ptr, int size){ // allocate
		if(logflags&LOGMEMORY)
			logger->sendmemlog("malloc ("+get_name()+") "+
			        as_string(ptr)+" "+as_string(size));
	}
	
	virtual void logmem(void* oldptr, void* newptr, int size){ // reallocate
		if(logflags&LOGMEMORY)
			logger->sendmemlog("realloc ("+get_name()+") "+
					as_string(oldptr)+" "+as_string(newptr)+" "+as_string(size));
	}

	void die(string message) {
		std::cerr<<"pid:"<<getpid()<<" ("<<get_name()<<")=> "<<message<<std::endl;
		throw "die...";
	}

	void die(string message, string data) {
		std::cerr<<"pid:"<<getpid()<<" ("<<get_name()<<")=> "
		         <<message<<": "<<data<<std::endl;
		throw "die...";
	}

	void assert(bool cond, string message) {
		if(!cond) die(message);
	}

	void assert(bool cond, string message, string val){
		if(!cond) die(message+val);
	}
	
	void warn(string message) {
		if(logflags&LOGWARNING) sendlog("warning ("+get_name()+"): "+message);
	}

	bool file_exist(string fname){
		struct stat filestat;
		if(stat(fname.c_str(), &filestat)) return false;
		return true;
	}

	string as_string(int val, const char* format=NULL){
		char st[256];
		if(format)sprintf(st,format,val);
		else sprintf(st,"%d",val);
		return string(st);
	}

	string as_string(OMD_FLOAT val, const char* format=NULL){
		char st[256];
		if(format)sprintf(st,format,val);
		else sprintf(st,"%0.2f",val);
		return string(st);
	}

	string as_string(void* val){
		char st[256];
		sprintf(st,"%lX",(long int)val);
		return string(st);
	}
  
  double as_double(string str) {
    char st[256],*pt;
    memset(st,0,256);
    str.copy(st,256);
    double retd=strtod(st,&pt);
    if(pt==st) die("invalid double value: "+str);
    return retd;
  }
  
  int as_int(string str) {
    char st[256],*pt;
    memset(st,0,256);
    str.copy(st,256);
    int reti=strtol(st,&pt,10);
    if(pt==st) die("invalid integer value: "+str);
    return reti;    
  }

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
  
  bool char_exist(string str, char ch) {
    for(int i=0;i<str.length();i++)
      if(str[i]==ch) return true;
    return false;
  }

	// class polling
	bool type_of(string type) {
    	type=lower_case(type);
    	type=replace_char(type, ' ', '_');		
		for(int i=0;i<(int)class_type.size();i++) {
			if(class_type[i]==type) return true;
		}
		return false;
	}

	string get_type() {
		return class_type.back();
	}

	string descent_type() {
		string ty;
		for(int i=0;i<(int)class_type.size();i++){
			ty.append(class_type[i]+" ");
		}
		return ty;
	}

};


#endif
