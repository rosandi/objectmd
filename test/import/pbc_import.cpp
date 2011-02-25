#include <omd/system.hpp>
#include <omd/integrator.hpp>
#include <omd/team.hpp>
#include <class/FCC100.hpp>
#include <class/NeighborCell.hpp>
#include <class/TempPressDetector.hpp>
#include <class/SysMonitor.hpp>
#include <omd/paramhandler.hpp>

string datafile;

class creator:public MDSystem {
	void CreateSystem() {
		Import(datafile);
	}
	
	void CreateGadget() {
		SetIntegrator(new MDIntegrator);
		AddForce(new TForceEAM("aluminum"));
		AddConditioner(new NeighborCell);
		AddDetector(new SysMonitor("md-import-s.out"));
		AddDetector(new TempPressDetector(0.1));
	}
	
	void SystemSetting() {
		SetWriteMode(WM_ID|WM_XID);
		MaxTime=2.0;
		SetOutputDirectory("output");
	}

	void BeforeRun() {
		PrintInfo("pbc_info.out");
		DumpAtoms("pbc_init.dat");
	}
	
	void AfterRun() {
		SaveSimulation();
	}

public:
    creator(int& argc, char** &argv):MDSystem(argc,argv){
		param.read(search_path("$OMD_TABLE", "def.aluminum"));
	}

};

int main(int argc, char* argv[]) {
OMD_BEGIN
	ParamHandler p(argc,argv);
	datafile=p.string_value("-i");
	creator c(argc, argv);
	return c.Run();
OMD_END
}
