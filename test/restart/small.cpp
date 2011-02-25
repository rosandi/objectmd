#include <iostream>
#include <unistd.h>
#include <omd/system.hpp>
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

int mode=NORMAL_MODE;
double maxt=1.0;

class MySim:public MDSystem {

	void CreateSystem() {
		CrystalFCC100 a(32,32,10, "platinum");
		a.Create()->Shift(-5.,-5.,50.0)->AddVelocity(0.0, 0.0, -10.0);

		CrystalFCC111 b(32,32,10, "platinum");
		b.Create();
		   
		AddAtom(new AtomContainer(a,b))->SetName("system");
	}
	
	void CreateGadget() {
		SetIntegrator(new MDIntegrator);
		AddForce(new TForceEAM("platinum"));
		AddConditioner(new NeighborCell);
		AddConditioner(new DynamicTimeStep(0.1));
		AddDetector(new RestartSaver(1.0, "minute"));
		AddDetector(new SysMonitor);
		AddDetector(new ThermoDetector);
		AddDetector(new StructureDetector)->Join(SearchGadget("temperature pressure detector"));
	}
	
	void SystemSetting() {
		PBoundary=NONPERIODIC;
		MaxTime=maxt;
	}
	
	void BeforeRun() {
		PrintInfo("info.out");
	}
	
	void AfterRun() {
		PrintTime(std::cout);
	}
};

int main(int argc, char* argv[]) {
	ParamHandler p(argc, argv);
	if(p.exist("-restart")) mode=RESTART_MODE;
	p.peek("-t", maxt);
   	MySim TheSim;
   	TheSim.set_name("OMD TEST");
	return TheSim.Run(mode);
}
