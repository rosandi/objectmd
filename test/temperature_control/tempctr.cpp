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
#include <omd/systemgrid.hpp>
#include <omd/team.hpp>
#include <class/NeighborCell.hpp>
#include <class/FCC100.hpp>
#include <class/TempController.hpp>
#include <class/TempPressDetector.hpp>
#include <class/SysMonitor.hpp>

class MySim:public MDSystemGrid {
	
	
	void CreateSystem() {
	    AddAtom(new CrystalFCC100(16, 16, 16, "aluminum"))
	    ->Create()
	    ->SetTemperature(100);
	}
	
	void CreateGadget() {
		SetIntegrator(new MDIntegrator);

		AddForce(new TForceEAM("aluminum"));
	
		AddConditioner(new NeighborCell);
		AddConditioner(new TempController(300.0,0.1));
		AddDetector(new SysMonitor("md.out"));
		AddDetector(new TempPressDetector(0.1,"Data"));
		
	}

	void SystemSetting() {
		PBoundary=PERIODIC_X|PERIODIC_Y|PERIODIC_Z;
		MaxTime=10.0;
		BorderOffset(0.25*param.double_value("lattice_constant"));
	}

	void BeforeRun() {
		PrintInfo("info.out");
		DumpAtoms("init.dat");
	}	

public:

	MySim(int& argc, char** &argv):MDSystemGrid(argc,argv){
		param.read(search_path("$OMD_TABLE", "def.aluminum"));
	}
};

int main(int argc, char* argv[]) {
	MySim TheSim(argc,argv);
	TheSim.SetClusterArch(2,2,1);
	return TheSim.Run();
}
