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
 *      Table reader
 *      The spline data is stored in memory
 *
*/

#ifndef _TABLE_READER_
#define _TABLE_READER_

#include <omd/base.hpp>

class TableReader:public MDClass {
	// Coefficient a is the value it self
	struct CoefStruct {OMD_FLOAT a, b, c, d;} *coef;
	enum tabel_format {plain, rmult} format;
	int	n_data;
	bool allow_outrange_low;
	bool allow_outrange_hi;
	OMD_FLOAT outrange_low;
	OMD_FLOAT outrange_hi;
	OMD_FLOAT  dr, roffset;
	void calculate_coef();
	void allocate();
	void open_omd(string, string);
	bool ready;

	public:
		string  filename;
		string  tablename;
		TableReader(string table_filename="", string table_name="");
		~TableReader();
		void   open(string table_filename, string table_name="");
		OMD_FLOAT max_range() {return (roffset+(OMD_FLOAT)(n_data-1)*dr);}
		OMD_FLOAT min_range() {return roffset;}
		bool   is_ready() {return ready;}
		OMD_FLOAT read(OMD_FLOAT r);
		void   read(OMD_FLOAT r,OMD_FLOAT& val,OMD_FLOAT& dval);
		OMD_FLOAT dread(OMD_FLOAT r);
		OMD_FLOAT dread2(OMD_FLOAT r);
		void   dump(string filename, int resolution=1000);
		void   dump_var();
};

#endif
