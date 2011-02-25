#include <omd/system.hpp>
#include <omd/integrator.hpp>
#include <class/FCC100.hpp>
#include <class/DummyForce.hpp>
#include <class/NeighborCell.hpp>

class creator:public MDSystem {
	void CreateSystem() {
		Import("system.cry");
	}
	
	void CreateGadget() {
		SetIntegrator(new MDIntegrator);
		AddForce(new DummyForce, "hot");
		AddForce(new DummyForce, "cool");
		AddForce(new DummyForce, "hot", "cool");
		AddConditioner(new NeighborCell);
	}
	
	void SystemSetting() {
		SetWriteMode(WM_ID|WM_XID);
	}

	void BeforeRun() {
		DumpAtoms("output.cry");
		die("normal termination");
	}
	
};

int main(int argc, char* argv[]) {
	creator c;
	return c.Run();
}
