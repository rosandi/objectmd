#include <omd/systemgrid.hpp>   
#include <omd/paramhandler.hpp>
#include <conditioner/VerletList.hpp>
#include <detector/LocalOrderParameter.hpp>

string cryfile;
string outfile("out");
double rcut;
double lconst;
int pbc=3;

class MyMDClass:public MDSystemGrid {

    void CreateSystem() {
      ReadMaterial("aluminum");
      Import(cryfile);
    }

    void CreateGadget() {
        param.append("--silent verlet.radius "+as_string(rcut));
        AddConditioner(new VerletList);
        AddDetector(new LocalOrderParameter(outfile, lconst, rcut, false));
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
  p.peek("-pbc", pbc);
  p.peek("-o", outfile);
} catch (const char* st) {
  std::cerr << st << "\n";
  return EXIT_FAILURE;
}
  
  std::cerr <<"r_cut=" << rcut << " latt_const=" << lconst 
            <<" pbc="<<pbc<<"\n";

  MyMDClass TheSim;
  TheSim.SetArgument(argc,argv);
  return TheSim.Run(STATIC_MODE);

}
