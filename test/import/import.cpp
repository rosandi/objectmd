/* app
 ************************************************
 * Object-MD
 * Object Oriented Molecular Dynamics Program
 * 
 * (c)Yudi Rosandi, 2005
 *
 * Version 2.0 (10.11.08)
 *
 * Project started on July 2005.
 *
 ************************************************
 *
*/

#include <iostream>
#include <unistd.h>
#include <fstream>

#include <omd/container.hpp>
#include <omd/paramhandler.hpp>

int main(int argc, char* argv[]) {
		
		ParamHandler p(argc, argv);
		
		if(p[1]==""){
			std::cerr<<"file name needed\n";
			return 1;
		}
		
		AtomContainer A,B;
		A.Import(p[1], 0)->DumpAtoms("type0");
		B.Import(p[1], 1)->DumpAtoms("type1");
	
	    AtomContainer C("aluminum");
	    C.Import(p[1])
	    	->SetName("imported")
	    	->SetWriteMode(WM_XID)
	    	->DumpAtoms();	    	
}
