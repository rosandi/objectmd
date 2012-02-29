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

#include <string>

using std::string;

#ifndef _OMD_TOOL_H_
#define _OMD_TOOL_H_

namespace omd {

  void die(string message);
  void die(string message, string val);
  void mdassert(bool cond, string message);
  void mdassert(bool cond, string message, string val);
  void warn(string message);
  bool file_exist(string fname);
  string as_string(int val, const char* format=NULL);
  string as_string(double val, const char* format=NULL);
  string as_string(void* val);
  int    as_int(string);
  double as_double(string);
  string lower_case(string str);
  string trim(string str);
  string replace_char(string str, char cold, char cnew);
  string replace_char(string, string, char);
  string remove_char(string, char);
  string remove_char(string, string);
  bool char_exist(string&, char);

}

#endif
