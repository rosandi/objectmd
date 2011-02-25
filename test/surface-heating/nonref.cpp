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
#include <omd/simgrid.h>
#include <mat/Aluminum>
#include <omd/fmmea.h>
#include <class/NeighborCell.hpp>
#include <class/CrystalFCC100.hpp>
#include <class/TempPressDetector.hpp>
#include <class/StructureDetector.hpp>
#include <class/SysMonitorGrid.hpp>

class fmm_detect:public fmmea {
	int auxi;
public:

	void Init(SimSystem* WorkSys) {
		fmmea::Init(WorkSys);
		auxi=ClaimAuxVariable(this,true,"fbottom");
	}

	void ClearAccumulators(){
		int na=GetNAtom();
		for(int i=0; i<na; i++) Atoms(i).aux[auxi]=0.0;
	}

	void ReturnForce(int at, int to, double dx, double dy, double dz, double fr, double pot){
		fmmea::ReturnForce(at,to,dx,dy,dz,fr,pot);
		if(Atoms(at).z>Atoms(to).z)
			Atoms(at).aux[auxi]+=fr*dz;
		else
			Atoms(to).aux[auxi]-=fr*dz;
	}

	fmm_detect(string tab):fmmea(tab){}
};

class MySim:public SimSystemGrid {

	void CreateSystem() {
		CrystalFCC100 h(16,16,10, LC_Al, MASS_Al, Z_Al);
		h.Create();
		h.SetName("Al-heated");
		h.SetTemperature(3000.0);
		
		CrystalFCC100 c(16,16,118, LC_Al, MASS_Al, Z_Al);
		c.Create();
		c.SetName("Al-cool");
		c.Shift(0.0,0.0,-5.0*LC_Al);

		AddAtom(new AtomContainer(h,c));
	}
	
	void CreateGadget() {
		SetIntegrator(new Integrator);
		AddForce(new fmm_detect("fmm-aluminum"));
		AddConditioner(new NeighborCell);		
		AddDetector(new SysMonitorGrid("md.out"));
		Detector* main_detector=new TempPressDetector(0.05,"Data");
		AddDetector(main_detector);
	}

	void SystemSetting() {
		TimeStep=1e-3;
		DynamicTimeStep=false;
		PBoundary=PERIODIC_X+PERIODIC_Y;
		MaxTime=20.0;
		BorderOffset(0.25*LC_Al);
		SetWriteMode("fv");
	}
	
	void BeforeRun() {
		PrintInfo(std::cout);
		DumpAtoms("init_data");
//		die("panic mode");
	}	
};

MD_MAIN_BEGIN
	MySim TheSim;
	TheSim.SetArgument(argc,argv);
	TheSim.SetClusterArch(2,2,1);
	TheSim.SetOutputDirectory("output");
	TheSim.Run();
MD_MAIN_END
