#include <omd/systemgrid.hpp>
#include <omd/container.hpp>
#include <potential/sw.hpp>
#include <crystal/Diamond.hpp>
#include <conditioner/VerletListFull.hpp>
#include <detector/ThermoDetector.hpp>
#include <detector/SysMonitor.hpp>

class MyMDClass:public MDSystemGrid {

    void CreateSystem() {
        if(param.exist("load.crystal")) {
            AddAtom(new AtomContainer("silicon"))
              ->Import(param.string_value("load_crystal"))
              ->SetName("Crystal");
        } else {
            int xml=param.int_value("monolayer",0);
            int yml=param.int_value("monolayer",1);
            int zml=param.int_value("monolayer",2);
            string ori=param.string_value("orientation");
            AddAtom(new Diamond(ori,xml,yml,zml,"silicon"))
              ->SetName("Crystal");
        }
    }

    void CreateGadget() {
        SetIntegrator(new MDIntegrator);
        AddForce(new StillingerWeber("silicon"));
        AddConditioner(new VerletListFull);
        AddDetector(new SysMonitor);  
        AddDetector(new ThermoDetector);
    }
    
    void SystemSetting() {
      if(param.exist("temperature"))
        SetTemperature(param.double_value("temperature"));
    }
    
    void BeforeRun() {
        DumpAtoms("init.dat", WM_VELOCITY);
        PrintInfo("info.out");
    }
    
    void AfterRun() {
      DumpAtoms("final.dat", WM_VELOCITY|WM_NID);
      SaveSimulation("saved.bin");
    }
    
};

int main(int argc, char* argv[]) {
    MyMDClass TheSim;
    TheSim.SetArgument(argc,argv);
    return TheSim.Run();
}
