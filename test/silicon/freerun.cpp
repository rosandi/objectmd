#include <omd/system.hpp>
#include <omd/container.hpp>
#include <potential/sw.hpp>
#include <potential/DummyForce.hpp>
#include <conditioner/VerletList.hpp>
#include <detector/ThermoDetector.hpp>
#include <detector/SysMonitor.hpp>

class MyMDClass:public MDSystem {

    void CreateSystem() {
        AtomContainer *A=new AtomContainer("silicon");
        A->Import(param.string_value("load_crystal"));
        if(param.exist("temperature"))
            A->SetTemperature(param.double_value("temperature"));
        AddAtom(A)->SetName("Crystal");
    }

    void CreateGadget() {
        SetIntegrator(new MDIntegrator);
        AddForce(new StillingerWeber("silicon"));
        AddConditioner(new VerletList);
        AddDetector(new SysMonitor);  
        AddDetector(new ThermoDetector);
    }

    void BeforeRun() {
        DumpAtoms("init.dat");
        PrintInfo("info.out");
    }
    
};

int main(int argc, char* argv[]) {
    MyMDClass TheSim;
    TheSim.SetArgument(argc,argv);
    return TheSim.Run();
}
