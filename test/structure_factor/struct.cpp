#include <omd/system.hpp>   
#include <omd/paramhandler.hpp>
#include <class/NeighborCell.hpp>
#include <class/StructureFactor.hpp>

string cryfile;
string outfile("out");
int nbin=100;
int pbc=3;
double rmax=-1.0;
double dsurf=0.0;

class MyMDClass:public MDSystem {

    void CreateSystem() {
      ReadMaterial("aluminum");
      Import(cryfile);
    }

    void CreateGadget() {
        AddDetector(new StructureFactor(outfile, false, nbin, 10, rmax))
          ->param.push_pair("skip_surface", as_string(dsurf));
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
  nbin=p.int_value("-bin");
  cryfile=p.string_value("-i");
  p.peek("-o", outfile);
  p.peek("-rmax", rmax);
  p.peek("-pbc", pbc);
  p.peek("-d", dsurf);
  
  std::cerr << "bin="<<nbin
      <<" R_max="<<rmax
      <<" boundary_condition="<<pbc
      <<" skip distance="<<dsurf
      <<"\n";
  
} catch (const char* st) {
  std::cerr << st << "\n";
  return EXIT_FAILURE;
}

  MyMDClass TheSim;
  TheSim.SetArgument(argc, argv);
  return TheSim.Run(STATIC_MODE);

}
