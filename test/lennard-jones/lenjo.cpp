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
#include <omd/simgrid.h>
#include <omd/tpair.h>
#include <class/NeighborCell.hpp>
#include <class/CrystalFCC100.hpp>
#include <class/TempPressDetector.hpp>
#include <class/StructureDetector.hpp>
#include <class/SysMonitorGrid.hpp>
#include <fstream>

class MySim:public SimSystemGrid {
	
	
	void CreateSystem() {
		param.read(search_path("$OMD_MATERIAL", "argon"));
		
		AddAtom(new CrystalFCC100(16,16,16, "argon"))
		->Create()
		->SetTemperature(30)
	    ->SetName("argon");

	}
	
	void CreateGadget() {
		SetIntegrator(new MDIntegrator);
		AddForce(new TForcePair("lj-ar.tab"));
		AddConditioner(new NeighborCell);		
		AddDetector(new SysMonitorGrid("md.out"));
		Detector* main_detector=new TempPressDetector(0.05,"Data");
		AddDetector(main_detector);
		AddDetector(new StructureDetector)->Join(main_detector);
		
	}

	void SystemSetting() {
		PBoundary=PERIODIC_X+PERIODIC_Y;
		MaxTime=5.0;
		BorderOffset(0.25*param.double_value("lattice_constant"));
	}

	void BeforeRun() {
		std::ofstream finfo("info.out");
		PrintInfo(finfo);
		finfo.close();
		DumpAtoms("init_data");
	}	

public:
	MySim(int& argc, char** &argv):SimSystemGrid(argc,argv){}
};

MD_MAIN_BEGIN
	MySim TheSim(argc,argv);
	TheSim.SetClusterArch(2,1,1);
	TheSim.Run();
MD_MAIN_END
