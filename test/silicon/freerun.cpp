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
            AddAtom(new Diamond("100",
                    param.double_value("monolayer",0),
                    param.double_value("monolayer",1),
                    param.double_value("monolayer",2),
                    "silicon"))
                    
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
    
};

int main(int argc, char* argv[]) {
    MyMDClass TheSim;
    TheSim.SetArgument(argc,argv);
    return TheSim.Run();
}
