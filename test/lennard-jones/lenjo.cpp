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
#include <potential/tpair.hpp>
#include <conditioner/VerletList.hpp>
#include <crystal/FCC100.hpp>
#include <detector/ThermoDetector.hpp>
#include <detector/StructureDetector.hpp>
#include <detector/SysMonitor.hpp>
#include <fstream>

class MySim:public MDSystemGrid {
	
	void CreateSystem() {
		param.read(search_path("$OMD_TABLE", "def.argon"));
		
		AddAtom(new CrystalFCC100(16,16,16, "argon"))
		->Create()
		->SetTemperature(30)
	    ->SetName("argon");

	}
	
	void CreateGadget() {
		SetIntegrator(new MDIntegrator);
		AddForce(new TForcePair("argon.lj"));
		AddConditioner(new VerletList);		
		AddDetector(new SysMonitor);
		Detector* main_detector=new ThermoDetector(0.05);
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
	MySim(int& argc, char** &argv):MDSystemGrid(argc,argv){}
};

int main(int argc, char* argv[]) {
	MySim TheSim(argc,argv);
	TheSim.SetClusterArch(2,1,1);
	TheSim.Run();
}
