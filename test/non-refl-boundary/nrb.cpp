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

#include <omd.hpp>
#include <conditioner/NonReflecting.h>

#define ALPHA  -0.0818
#define FSTAT   0.1439
#define ZML 2.0

class MySim:public MDSystemGrid {
	double lc;

	void CreateSystem() {
		AddAtom(new FCC("100",16,16,132, "aluminum"))->SetName("Target");    
	}
  
  void CreateGroup() {
    AddAtomGroup("heated")
      ->SelectGT(0.0,0.0,Box.z1-10.0*ZML)
      ->SetTemperature(1000.0);
    AddAtomGroup("nrb")->SelectLT(DBL_MAX,DBL_MAX,Box.z0+4.0*ZML);
  }

	void CreateGadget() {
		SetIntegrator(new MDIntegrator);

		NonReflecting* noref=new NonReflecting;
		AddForce(new TForceEAM("aluminum"))->SetEvaluator(noref);

		AddConditioner(new VerletList);
		AddConditioner(noref);

		AddDetector(new SysMonitor("md.out"));
		AddDetector(new ThermoDetector(0.05));
	}

	void BeforeRun() {
		PrintInfo("info.out");
		DumpAtoms("init_data");
	}

};

int main(int argc, char* argv[]) {
	MySim TheSim;
	TheSim.SetArgument(argc,argv);
	return TheSim.Run();
}
