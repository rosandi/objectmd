#include <omd/system.hpp>
#include <omd/integrator.hpp>
#include <crystal/FCC100.hpp>
#include <potential/DummyForce.hpp>
#include <conditioner/NeighborCell.hpp>

class creator:public MDSystem {
	
	void CreateGadget() {
		SetIntegrator(new MDIntegrator);
		AddForce(new DummyForce, "hot");
		AddForce(new DummyForce, "cool");
		AddForce(new DummyForce, "hot", "cool");
		AddConditioner(new NeighborCell);
	}
	
	void BeforeRun() {
		PrintInfo("simtest.info");
		DumpAtoms("output.cry");
		die("normal termination");
	}
	
};

int main(int argc, char* argv[]) {
	if(argc<3) {
		std::cerr<<"syntax:\n./simtest paramfile simtest.par\n";
		exit(1);
	}
	creator c;
	c.SetArgument(argc,argv);
	return c.Run();
}
