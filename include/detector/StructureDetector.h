//-------------------------Ackland Structure Detector -----------------------//

#ifndef _STRUCTURE_DETECTOR_HPP_
#define _STRUCTURE_DETECTOR_HPP_

#include <vector>
#include <detector/DataDumper.h>

struct struct_neig {int idx; double rsq;};

/**
  @ingroup detector
  @brief Ackland's near order structure detector
  
  The class detects the local structure of an atom in a crystal.
  The algorithm used is the ackland near order method.
  
  ref:
    - Ackland, G. J. & Jones, A. P. "Applications of local crystal structure 
      measures in experiment and simulation",  Phys. Rep. B, 2006, 73, 054104
    - Gerolf Ziegenhain, the same detector \@lammps.

**/

class StructureDetector: public DataDumper {
	
	vector<struct_neig> *neig;
	int nalloc,stidx;
	double sqrcut;
	
public:
	StructureDetector(double tm=0.0, string fn="Data"):
	DataDumper("",tm,fn,3) {
		set_name("cstr");
		register_class("STRUCTURE_DETECTOR");
		neig=NULL;
	}
	
	virtual ~StructureDetector(){MemDeleteArray(neig);}
	void Init(MDSystem* WorkSys);
	void Measure(){FindStructure();DataDumper::Measure();}
	void IterationNode(int at, int to);
	void FindStructure();
};

#endif
