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

#include <omd.h>

class MySim:public MDSystemGrid {
	void CreateSystem() {
		AddAtom(new FCC("111",30,20,10, "platinum"))
		   ->Create()
		   ->SetName("Crystal");

/*		   
		AddAtom(new FreeAtom(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, MASS_Pt, Z_Pt));
		AddAtom(new FreeAtom(2.0, 0.0, 0.0, 0.0, 0.0, 0.0, MASS_Pt, Z_Pt));
		AddAtom(new FreeAtom(0.0, 2.0, 0.0, 0.0, 0.0, 0.0, MASS_Pt, Z_Pt));
		AddAtom(new FreeAtom(2.0, 2.0, 0.0, 0.0, 0.0, 0.0, MASS_Pt, Z_Pt));
*/
	}
	
	void CreateGadget() {
		SetIntegrator(new MDIntegrator);
		AddForce(new TForceEAM("platinum"));
//		AddConditioner(new NeighborCell);
//		AddConditioner(new VerletListFull);
		AddConditioner(new VerletList);
		AddDetector(new SysMonitor("md.out"));
		AddDetector(new ThermoDetector(0.1));
	}
	
	void BeforeRun() {
		DumpAtoms("init.dat");
		PrintInfo("info.out");
	}

};


int main(int argc, char* argv[]) {
	MySim TheSim;
	TheSim.SetArgument(argc,argv);
	return TheSim.Run();
}
