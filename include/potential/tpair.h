/* omdlib
 ********************************************************
 *
 * ObjectMD
 * Molecular Dynamics Class Library
 * Version 2.0 (2009)
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
 *       Pair potential using splined table
 *
*/

#ifndef _T_PAIR_H_
#define _T_PAIR_H_

#include <vector>
#include <omd/gadget.h>
#include <omd/forcekernel.h>
#include <omd/treader.h>

/**
 * @brief Table splining pair potential
 * 
 * This is the force kernel for pair interaction potential. The potential
 * is loaded from a table. See @link functable the Object-MD table @endlink
 * chapter.
 * 
 */

class TForcePair: public ForceKernel {
	TableReader phi;
	string tablefile;
	
	public:	
		
		TForcePair(string TableFile) {
			char sname[128];
			set_name("TABLE PAIR FORCE");
			register_class(get_name());
			phi.name.assign("READER@"+get_name());
			tablefile=TableFile;
		}
		
		void   PrintInfo(ostream& ost) 
		{ost << "ID." << id << " " << get_name() << "(cut radius = " << CutRadius << ", Table=" << phi.filename << ")\n";}
		
		void Init(MDSystem* WorkSys) {
			ForceKernel::Init(WorkSys);
			tablefile=WorkSys->search_path("$OMD_TABLE","pair."+tablefile);
			phi.open(tablefile);
			CutRadius=phi.max_range();
		}

		void ComputeHalf(Atom &at, Atom &to) {
			double ff, dx, dy, dz, r;
			double pot,drpot;
			r = sqrt(CalcSqrDistance(at, to, dx, dy, dz));
			if (r<CutRadius) {
				phi.read(r, pot, drpot);
				ff=-drpot/r;
				ReturnForce(at,to,dx,dy,dz,ff,pot);
			}
		}
};

#endif
