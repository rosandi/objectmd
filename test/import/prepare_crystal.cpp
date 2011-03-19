#include <omd/system.hpp>
#include <omd/integrator.hpp>
#include <crystal/FCC100.hpp>
#include <potential/DummyForce.hpp>
#include <conditioner/NeighborCell.hpp>

class creator:public MDSystem {
	void CreateSystem() {
		AtomContainer* A=AddAtom(new CrystalFCC100(24,24,32, "argon"))
			->Create()
			->SetTemperature(800.0)
			->SetID(0)->SetXID(1)
			->SetName("hot");
		
		double lc=A->param.double_value("lattice_constant");

		AddAtom(new CrystalFCC100(32,32,32, "platinum"))
			->Create()
			->Shift(0.0,0.0,-16.0*lc)
			->SetID(1)->SetXID(2)
			->SetName("cool");

	}
	
	void CreateGadget() {
		SetIntegrator(new MDIntegrator);
		AddForce(new DummyForce, "hot");
		AddForce(new DummyForce, "cool");
		AddForce(new DummyForce, "hot", "cool");
		AddConditioner(new NeighborCell);
	}
	
	void SystemSetting() {
		BorderOffset(10.0);
		SetWriteMode(WM_ID|WM_XID);
	}

	void BeforeRun() {
		DumpAtoms("system.cry");
		die("normal termination");
	}
	
};

int main(int argc, char* argv[]) {
	creator c;
	return c.Run();
}
