#include <omd/system.hpp>
#include <omd/container.hpp>
#include <potential/tersoff.hpp>
#include <potential/DummyForce.hpp>
#include <conditioner/NeighborCell.hpp>
#include <detector/ThermoDetector.hpp>
#include <detector/SysMonitor.hpp>

class MyMDClass:public MDSystem {

    void CreateSystem() {
        AtomContainer *A=new AtomContainer(param.string_value("tersoff.material"));
        A->Import(param.string_value("load_crystal"));
        if(param.exist("temperature"))
            A->SetTemperature(param.double_value("temperature"));
        AddAtom(A)->SetName("Crystal");
    }

    void CreateGadget() {
        SetIntegrator(new MDIntegrator);
        AddForce(new Tersoff);
        AddConditioner(new NeighborCell);
        AddDetector(new SysMonitor);  
        AddDetector(new ThermoDetector);
    }

    void SystemSetting() {
        PBoundary=NONPERIODIC;
        param.peek("maxtime", MaxTime, 1.0);
        SetOutputDirectory("output");
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
