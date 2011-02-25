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
#include <fstream>
#include <unistd.h>
#include <omd/systemgrid.hpp>
#include <omd/team.hpp>
#include <omd/tpair.hpp>
#include <class/NeighborCell.hpp>
#include <class/FCC100.hpp>
#include <class/TempPressDetector.hpp>
#include <class/StructureDetector.hpp>
#include <class/SysMonitor.hpp>

class MySim:public MDSystemGrid {
	
	
	void CreateSystem() {
		param.read(search_path("$OMD_TABLE", "def.platinum"));
		
		CrystalFCC100 CC(8,8,8, "argon");
		CC.Create();
		CC.Shift(5.,5.,4*5.3);

		CrystalFCC100 CD(16,16,8, "platinum");
		CD.Create();	    

	    CC.DumpAtoms("ar");
	    CD.DumpAtoms("pt");
	    
	    AddAtom(new AtomContainer)
	    ->Combine(CC)
	    ->SetName("argon");

	    AddAtom(new AtomContainer)
	    ->Combine(CD)
	    ->SetName("platinum");

	}
	
	void CreateGadget() {
		SetIntegrator(new MDIntegrator);

		AddForce(new TForceEAM("platinum"), "platinum");
		AddForce(new TForcePair("argon.lj"), "argon");
		AddForce(new TForcePair("ar-pt.lj"), "argon", "platinum");
	
		AddConditioner(new NeighborCell);		
		AddDetector(new SysMonitor("md.out"));
		Detector* main_detector=new TempPressDetector(0.1,"Data");
		AddDetector(main_detector);
		AddDetector(new StructureDetector)->Join(main_detector);
		
	}

	void SystemSetting() {
		PBoundary=PERIODIC_X+PERIODIC_Y;
		MaxTime=5.0;
		BorderOffset(0.25*param.double_value("lattice_constant"));
		SetWriteMode(WM_ID);
	}

	void BeforeRun() {
		PrintInfo("info.out");
		DumpAtoms("init.data");
	}	

public:
	MySim(int& argc, char** &argv):MDSystemGrid(argc,argv){}
};

int main(int argc, char* argv[]) {
	MySim TheSim(argc,argv);
	TheSim.SetClusterArch(2,2,1);
	return TheSim.Run();
}
