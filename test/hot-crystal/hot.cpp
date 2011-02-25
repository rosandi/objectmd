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
#include <class/TempPressDetector.hpp>
#include <class/StructureDetector.hpp>
#include <class/SysMonitor.hpp>

class MySim:public MDSystemGrid {
	
	
	void CreateSystem() {
		
	    AddAtom(new CrystalFCC100(32, 32, 32, "aluminum"))
	    ->Create()
	    ->SetTemperature(300.0);

	}
	
	void CreateGadget() {
		SetIntegrator(new MDIntegrator);

		AddForce(new TForceEAM("aluminum"));
	
		AddConditioner(new NeighborCell);		

		AddDetector(new SysMonitor("md.out"));
		AddDetector(new TempPressDetector(0.1,"Data"));
		
	}

	void SystemSetting() {
		PBoundary=PERIODIC_X|PERIODIC_Y|PERIODIC_Z;
		MaxTime=5.0;
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
	TheSim.SetClusterArch(2,1,1);
	return TheSim.Run();
}
