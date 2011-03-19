#include <crystal/FCC100.hpp>

int main(int argc, char* argv[]) {
	CrystalFCC100 CC(32,32,32, "platinum");
	CC.Create()
		->SetTemperature(800.0)
		->SetID(0)->SetXID(1)
		->SetName("hot")
		->Save("save.bin", "w");

	CrystalFCC100 CD(32,32,32, "platinum");
	CD.Create()
		->Shift(0.0,0.0,-16.0*CD.param.double_value("lattice_constant"))
		->SetID(1)->SetXID(2)
		->SetName("cool")
		->Save("save.bin","a");
}
