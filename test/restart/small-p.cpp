#include <iostream>
#include <unistd.h>
#include <omd/systemgrid.hpp>
#include <omd/gadget.hpp>
#include <potential/team.hpp>
#include <crystal/FCC111.hpp>
#include <crystal/FCC100.hpp>
#include <conditioner/NeighborCell.hpp>
#include <conditioner/DynamicTimeStep.hpp>
#include <detector/ThermoDetector.hpp>
#include <detector/StructureDetector.hpp>
#include <detector/SysMonitor.hpp>
#include <detector/RestartSaver.hpp>

class MySim:public MDSystemGrid {

	void CreateSystem() {
		CrystalFCC100 a(32,32,10, param.string_value("material"));
		a.Create()->Shift(-5.,-5.,50.0)->AddVelocity(0.0, 0.0, -10.0);

		CrystalFCC111 b(32,32,10, param.string_value("material"));
		b.Create();
		   
		AddAtom(new AtomContainer(a,b))->SetName("system");
	}
	
	void CreateGadget() {
		SetIntegrator(new MDIntegrator);
		AddForce(new TForceEAM(param.string_value("material")));
		AddConditioner(new NeighborCell);
		AddConditioner(new DynamicTimeStep);
		AddDetector(new RestartSaver);
		AddDetector(new SysMonitor);
		AddDetector(new ThermoDetector);
		AddDetector(new StructureDetector)
			->Join(SearchGadget("THERMO DETECTOR"));
	}
	
	void BeforeRun() {
		PrintInfo("info.out");
	}
	
	void AfterRun() {
		PrintTime(std::cout);
		SaveSimulation("end.bin");
	}
};

int main(int argc, char* argv[]) {
   	MySim TheSim;
   	TheSim.SetArgument(argc, argv);
   	TheSim.set_name("OMD PARALLEL TEST");
	return TheSim.Run();
}
