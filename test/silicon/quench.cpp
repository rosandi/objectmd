#include <omd/systemgrid.hpp>
#include <omd/container.hpp>
#include <potential/sw.hpp>
#include <crystal/Diamond.hpp>
#include <conditioner/VerletListFull.hpp>
#include <conditioner/CoordClamp.hpp>
#include <conditioner/Quencher.hpp>
#include <detector/ThermoDetector.hpp>
#include <detector/SysMonitor.hpp>

class MyMDClass:public MDSystemGrid {

    void CreateGroup() {

        AddAtomGroup("border")
          ->SelectInverseBox(Box.x0+2.4, Box.x1-2.4, 
                             Box.y0+2.4, Box.y1-2.4,
                             Box.z0+2.4, INFINITE);
                             
    }
    
    void CreateGadget() {
        SetIntegrator(new MDIntegrator);
        AddForce(new StillingerWeber("silicon"));
        AddConditioner(new VerletListFull);
        AddConditioner(new Quencher);
        AddConditioner(new CoordClamp("border"));
        AddDetector(new SysMonitor);  
        AddDetector(new ThermoDetector);
    }
    void BeforeRun() {
        DumpAtoms("quench-init.dat", WM_VELOCITY|WM_NID);
        DumpAtoms(SearchContainer("border")->GetAtomStorage(), "border.dat", WM_VELOCITY);
        PrintInfo("quench-info.out");
    }
                            
    void AfterRun() {
      DumpAtoms("quenched.dat", WM_VELOCITY|WM_NID);
      SaveSimulation("quenched.bin");
    }

};

int main(int argc, char* argv[]) {
    MyMDClass TheSim;
    TheSim.SetParameterFile("omd-quench");
    TheSim.SetArgument(argc,argv);
    return TheSim.Run();
}
