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
#include <omd/team.hpp>
#include <class/NeighborCell.hpp>
#include <class/FCC100.hpp>
#include <class/TempPressDetector.hpp>
#include <class/StructureDetector.hpp>
#include <class/SysMonitor.hpp>

class MySim:public MDSystemGrid {
	double lc;
	void CreateSystem() {

		CrystalFCC100 CC(16,16,16, "platinum");
		CC.Create();
		CC.SetName("hot");
		CC.SetTemperature(Unit->InKelvin(1000.0));

		CrystalFCC100 CD(16,16,16, "platinum");
		CD.Create();
		CD.SetName("cool");
	    CD.Shift(0.0,0.0,-16.0*lc);
	    CC.DumpAtoms("hot");
	    CD.DumpAtoms("cool");
	    
	    AddAtom(new AtomContainer(CC, CD));

	}
	
	void CreateGadget() {
		SetIntegrator(new MDIntegrator);
		AddForce(new TForceEAM("platinum"));
		AddConditioner(new NeighborCell);
		
		AddDetector(new SysMonitor("md.out"));
		Detector* main_detector=new TempPressDetector(0.01,"Data");
		AddDetector(main_detector);
		AddDetector(new StructureDetector)->Join(main_detector);
		
	}

	void SystemSetting() {
		PBoundary=PERIODIC_X+PERIODIC_Y+PERIODIC_Z;
		MaxTime=1.0;
		BorderOffset(0.25*lc);
	}
	
	void BeforeRun() {
		PrintInfo("info.out");
		DumpAtoms("init.data");
	}	
	
public:
	MySim() {
		ParamHandler p(search_path("$OMD_TABLE", "def.platinum"));
		lc=p.double_value("lattice_constant");
	}
};

int main(int argc, char* argv[]) {
	MySim TheSim;
	TheSim.SetArgument(argc,argv);
	TheSim.SetClusterArch(2,2,1);
	return TheSim.Run();
}

