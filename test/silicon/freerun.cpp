#include <omd/systemgrid.hpp>
#include <omd/container.hpp>
#include <potential/sw.hpp>
#include <conditioner/VerletListFull.hpp>
#include <detector/ThermoDetector.hpp>
#include <detector/SysMonitor.hpp>

class MyMDClass:public MDSystem {

    void CreateSystem() {
        AtomContainer *A=new AtomContainer("silicon");
        A->Import(param.string_value("load_crystal"));
//        if(param.exist("temperature"))
//            A->SetTemperature(param.double_value("temperature"));
        AddAtom(A)->SetName("Crystal");
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
