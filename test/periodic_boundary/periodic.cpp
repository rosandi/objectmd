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
#include <omd/systemgrid.hpp>
#include <potential/team.hpp>
#include <conditioner/NeighborCell.hpp>
#include <crystal/FCC100.hpp>
#include <detector/ThermoDetector.hpp>
#include <detector/StructureDetector.hpp>
#include <detector/SysMonitor.hpp>

#define LC 3.924

class MySim:public MDSystemGrid {	
	
	void CreateSystem() {

		param.read(search_path("$OMD_TABLE", "def.platinum"));
		
		CrystalFCC100 CC(16,16,32, "platinum");
		CC.Create()
			->SetTemperature(800.0)
			->SetName("hot")
			->DumpAtoms();

		CrystalFCC100 CD(16,16,32, "platinum");
		CD.Create()
			->Shift(0.0,0.0,-16.0*LC)
			->SetName("cool")->DumpAtoms();
	    
	    AddAtom(new AtomContainer(CC,CD))
// or like this:	    ->Combine(CC)
//              	    ->Combine(CD)
	    ->SetName("TargetCrystal");
//	    ->param.dump(std::cout);
	}
	
	void CreateGadget() {
		SetIntegrator(new MDIntegrator);
		AddForce(new TForceEAM("platinum"));
		AddConditioner(new NeighborCell);		
		AddDetector(new SysMonitor("md.out"));
//		Detector* main_detector=new ThermoDetector(0.01);
//		AddDetector(main_detector);
		AddDetector(new ThermoDetector(0.01));
//		AddDetector(new StructureDetector)->Join(main_detector);
		
	}

	void SystemSetting() {
//		BorderOffset(0.25*LC);
		AcceptSignal(SIGINT);
		AcceptSignal(SIGUSR1);
	}

	void BeforeRun() {
		PrintInfo("info.out");
		DumpAtoms("init_data");
	}
	
	void AfterRun() {
		SaveSimulation();
	}

};

int main(int argc, char* argv[]) {
	MySim TheSim;
	TheSim.SetArgument(argc,argv);
	return TheSim.Run();
}
