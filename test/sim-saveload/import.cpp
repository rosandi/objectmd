#include <omd/systemgrid.hpp>
#include <omd/integrator.hpp>
#include <potential/team.hpp>
#include <crystal/FCC100.hpp>
#include <conditioner/VerletList.hpp>
#include <detector/ThermoDetector.hpp>
#include <detector/SysMonitor.hpp>
#include <omd/paramhandler.hpp>

class creator:public MDSystemGrid {

	void CreateGadget() {
		SetIntegrator(new MDIntegrator);
		AddForce(new TForceEAM("aluminum"));
		AddConditioner(new VerletList);
		AddDetector(new SysMonitor("import.out"));
		AddDetector(new ThermoDetector(0.1));
	}
	
	void BeforeRun() {
		PrintInfo("import.info");
		DumpAtoms("import-init.cry");
	}
	
	void AfterRun() {
		SaveSimulation("import.bin");
	}

};

int main(int argc, char* argv[]) {
    creator c;
    c.SetArgument(argc,argv);
    return c.Run();
}
