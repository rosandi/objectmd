#include <omd/systemgrid.hpp>
#include <omd/integrator.hpp>
#include <potential/team.hpp>
#include <omd/paramhandler.hpp>
#include <crystal/FCC100.hpp>
#include <conditioner/VerletList.hpp>
#include <detector/ThermoDetector.hpp>
#include <detector/SysMonitor.hpp>
#include <conditioner/TempController.hpp>

class creator:public MDSystemGrid {
	void CreateSystem() {
		AddAtom(new CrystalFCC100(16, 16, 16, "aluminum"))
			->Create()
			->SetTemperature(param.double_value("temperature"))
			->SetName("target");
	}
	
	void CreateGadget() {
		SetIntegrator(new MDIntegrator);
		AddForce(new TForceEAM("aluminum"));
		AddConditioner(new VerletList);
		AddConditioner(new TempController(param.double_value("temperature"), 0.5));
		AddDetector(new SysMonitor("export.out"));
		AddDetector(new ThermoDetector(0.1));
	}
	
	void BeforeRun() {
		PrintInfo("pbc_export.info");
		DumpAtoms("pbc.init");
	}
	
	void AfterRun() {
		SaveSimulation("export.bin");
		DumpAtoms("pbc_export.cry");
	}

};

int main(int argc, char* argv[]) {
	if(argc==1) {
		std::cerr<<"syntax:\n ./pbc_export --param pbc.par\n";
		exit(1);
	}
	creator c;
	c.SetArgument(argc,argv);
	return c.Run();
}
