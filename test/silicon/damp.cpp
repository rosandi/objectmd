#include <omd/systemgrid.hpp>
#include <omd/container.hpp>
#include <potential/sw.hpp>
#include <crystal/Diamond.hpp>
#include <conditioner/VerletListFull.hpp>
#include <conditioner/CoordClamp.hpp>
#include <conditioner/ForceDamper.hpp>
#include <detector/ThermoDetector.hpp>
#include <detector/SysMonitor.hpp>

class MyMDClass:public MDSystemGrid {

    void CreateGroup() {

        AddAtomGroup("border")
          ->SelectInverseBox(Box.x0+2.4, Box.x1-2.4, 
                             Box.y0+2.4, Box.y1-2.4,
                             Box.z0+2.4, INFINITE);
                             
       AddAtomGroup("damp-area")
          ->SelectInverseBox(Box.x0+5.43, Box.x1-5.43, 
                             Box.y0+5.43, Box.y1-5.43,
                             -INFINITE, INFINITE);
//          ->Commit()
//          ->DumpAtoms("damp-area-"+as_string(GetRank()), WM_NID);
    }
    
    void CreateGadget() {
        SetIntegrator(new MDIntegrator);
        AddForce(new StillingerWeber("silicon"));
        AddConditioner(new VerletListFull);
        AddConditioner(new ForceDamper("damp-area"));
        AddConditioner(new CoordClamp("border"));
        AddDetector(new SysMonitor);  
        AddDetector(new ThermoDetector);
    }
    void BeforeRun() {
        PrintInfo("damp-info.out");
    }
                            
    void AfterRun() {
      DumpAtoms("damped.dat", WM_VELOCITY|WM_NID);
      SaveSimulation("damped.bin");
    }

};

int main(int argc, char* argv[]) {
    MyMDClass TheSim;
    TheSim.SetParameterFile("omd-damp");
    TheSim.SetArgument(argc,argv);
    return TheSim.Run();
}
