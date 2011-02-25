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
#include <omd/systemgrid.hpp>
#include <omd/paramhandler.hpp>
#include <omd/team.hpp>
#include <class/NeighborCell.hpp>
#include <class/FCC100.hpp>
#include <class/TempPressDetector.hpp>
#include <class/StructureDetector.hpp>
#include <class/SysMonitor.hpp>

#define LC 4.032

int mode=NORMAL_MODE;

class MySim:public MDSystemGrid {
	void CreateSystem() {

		AddAtom(new CrystalFCC100(16,16,8, "aluminum"))
		   ->Create()
		   ->SetKineticEnergy(0.8);

	}
	
	void CreateGadget() {
		SetIntegrator(new MDIntegrator);
		AddForce(new TForceEAM("aluminum"));
		AddConditioner(new NeighborCell);

		AddDetector(new SysMonitor("md.out"));
		Detector* main_detector=new TempPressDetector(0.01,"Data");
		AddDetector(main_detector);
		AddDetector(new StructureDetector)->Join(main_detector);
		
	}
	
	void SystemSetting() {
		PBoundary=NONPERIODIC;
		MaxTime=5.0;
		BorderOffset(0.25*LC);
		AcceptSignal(SIGINT); // ctrl-c works only on serial/single proc
		AcceptSignal(SIGUSR1); 
		SetOutputDirectory("output");
	}

	void BeforeRun() {
		PrintInfo("info.out");
		DumpAtoms("data.init");
	}

};

int main(int argc, char* argv[]) {
	ParamHandler p(argc,argv);
	if(p.exist("-restart")) mode=RESTART_MODE;
	
	MySim TheSim;
	TheSim.SetArgument(argc,argv);
	TheSim.SetName("signal test");
	TheSim.SetClusterArch(GRID_AUTOZ);
	return TheSim.Run(mode);
}
