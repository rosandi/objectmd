#include <omd/system.hpp>
#include <omd/integrator.hpp>
#include <omd/team.hpp>
#include <omd/paramhandler.hpp>
#include <class/FCC100.hpp>
#include <class/NeighborCell.hpp>
#include <class/TempPressDetector.hpp>
#include <class/SysMonitor.hpp>
#include <class/TempController.hpp>

double tempe=55.0;

class creator:public MDSystem {
	void CreateSystem() {
		AddAtom(new CrystalFCC100(16, 16, 16, "aluminum"))
			->Create()
			->SetTemperature(tempe)
			->SetName("target");
	}
	
	void CreateGadget() {
		SetIntegrator(new MDIntegrator);
		AddForce(new TForceEAM("aluminum"));
		AddConditioner(new NeighborCell);
		AddConditioner(new TempController(tempe, 0.5));
		AddDetector(new SysMonitor("md-export-s.out"));
		AddDetector(new TempPressDetector(0.1));
	}
	
	void SystemSetting() {
		PBoundary=PERIODIC_X|PERIODIC_Y;
		SetWriteMode(WM_ID|WM_XID);
		MaxTime=1.0;
		SetOutputDirectory("output");
		BorderOffset(0.25*param.double_value("lattice_constant"));
	}

	void BeforeRun() {
		PrintInfo("pbc_info.out");
		DumpAtoms("pbc_init.dat");
	}
	
	void AfterRun() {
		SaveSimulation();
		DumpAtoms("pbc_export.dat");
	}

public:
    creator(int& argc, char** &argv):MDSystem(argc,argv) {
		param.read(search_path("$OMD_TABLE", "def.aluminum"));
	}

};

int main(int argc, char* argv[]) {
	creator c(argc, argv);
	ParamHandler p(argc,argv);
	p.peek("-t", tempe);
	c.SetName("export serial "+c.as_string(tempe));
	return c.Run();
}
