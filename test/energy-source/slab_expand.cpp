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

#include <omd/systemgrid.hpp>
#include <potential/team.hpp>
#include <conditioner/VerletList.hpp>
#include <conditioner/EnergySource.hpp>
#include <crystal/FCC100.hpp>
#include <detector/ThermoDetector.hpp>
#include <detector/SysMonitor.hpp>

class MySim:public MDSystemGrid {	
	
	void CreateSystem() {
	    AddAtom(new CrystalFCC100(16,16,64,"aluminum"))
	    	->Create()
	    	->SetName("Target");
	}
	
	void CreateGadget() {
		SetIntegrator(new MDIntegrator);
		AddForce(new TForceEAM("aluminum"));
		AddConditioner(new VerletList);
		AddConditioner(new EnergySource);		
		AddDetector(new SysMonitor("md.out"));
		AddDetector(new ThermoDetector(0.01));
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
