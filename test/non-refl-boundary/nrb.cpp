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

#define TERMID 13

class MySim:public MDSystemGrid {
	double lc;

	void CreateSystem() {
		FCC h("100",16,16,10, "aluminum");
		h.Create()
			->SetName("Al-heated")
			->SetTemperature(1000.0);
		
		FCC c("110",16,16,118, "aluminum");
		c.Create()
			->SetName("Al-cool")
			->Shift(0.0,0.0,-5.0*lc);

		FCC b("110",16,16,4, "aluminum");
		b.Create()
			->SetName("terzone")
			->Shift(0.0,0.0,-64.0*lc)
			->SetXID(TERMID);
		
		AddAtom(new AtomContainer)
			->Combine(h)
			->Combine(c)
			->Combine(b)
			->SetName("Target");
			
	}

	void CreateGadget() {
		SetIntegrator(new MDIntegrator);

		NonReflecting* noref=new NonReflecting(ALPHA, FSTAT, TERMID, lc/2.0, true);
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
