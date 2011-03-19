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

int main() {
    AtomContainer A,B;
    A.Load("save.bin", "hot")->DumpAtoms("hot.load");
    B.Load("save.bin", "cool")->DumpAtoms("cool.load");

/*	
	    AtomContainer C("aluminum");
	    C.Import(p[1])
	    	->SetName("imported")
	    	->SetWriteMode(WM_XID)
	    	->DumpAtoms();	    	
*/

}
