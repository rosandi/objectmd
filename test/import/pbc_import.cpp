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
		PrintInfo("pbc.info");
		DumpAtoms("pbc.init");
	}
	
	void AfterRun() {
		SaveSimulation("pbc_import.bin");
	}

};

int main(int argc, char* argv[]) {
    if(argc==1){
        std::cerr << "syntax:\n ./pbc_import --param pbc.par import filename\n";
        exit(1);
    } 
    creator c;
    c.SetArgument(argc,argv);
    return c.Run();
}
