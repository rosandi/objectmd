#include <omd/system.hpp>   
#include <omd/paramhandler.hpp>
#include <conditioner/NeighborCell.hpp>
#include <detector/ThermoDetector.hpp>

string cryfile;
string outfile("out");
string material;
double rcut;
double lconst;
int pbc=3;

class MyMDClass:public MDSystem {

    void CreateSystem() {
      ReadMaterial(material);
      Import(cryfile);
    }

    void CreateGadget() {
        AddConditioner(new NeighborCell(0, rcut));
        AddDetector(new ThermoDetector(0.0, rcut, outfile, 0));
    }

    void SystemSetting() {
        PBoundary=pbc;
    }
    
    void BeforeRun() {
      PrintInfo("info.out");
    }

};

int main(int argc, char* argv[]) {
try{
  ParamHandler p(argc, argv);
  rcut=p.double_value("-rcut");
  lconst=p.double_value("-lc");
  cryfile=p.string_value("-i");
  material=p.string_value("-mat");
  p.peek("-pbc", pbc);
  p.peek("-o", outfile);
} catch (const char* st) {
  std::cerr << st << "\n";
  return EXIT_FAILURE;
}
  
  std::cerr <<"r_cut=" << rcut 
  << " material=" << material
  << " latt_const=" << lconst 
  <<" pbc="  << pbc<<"\n";

  MyMDClass TheSim;
  return TheSim.Run(STATIC_MODE);

}
