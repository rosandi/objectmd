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
		AddAtom(new CrystalFCC100(16,16,64, "platinum"))
			->SetName("target")
			->Create();
	}
	
	void PostCreation() {
/*		AddAtomGroup("hot slab")
			->SelectGT(WEST_EDGE, SOUTH_EDGE, -16*LC)
			->SetTemperature(600);
*/
	}
	
	void CreateGadget() {
		SetIntegrator(new MDIntegrator);
		AddForce(new TForceEAM("platinum"));
		AddConditioner(new NeighborCell);		
		AddDetector(new SysMonitor("md.out"));
		AddDetector(new ThermoDetector(0.01));
	}

	void SystemSetting() {
		BorderOffset(0.25*LC);
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
