#include <omd/system.hpp>
#include <crystal/FCC111.hpp>
#include <potential/team.hpp>
#include <conditioner/VerletList.hpp>
#include <detector/ThermoDetector.hpp>
#include <detector/SysMonitor.hpp>

class MyMDClass:public MDSystem {

    void CreateSystem() {
        AddAtom(new CrystalFCC111(30,20,16, "platinum"))
           ->Create()
           ->SetTemperature(100.0)
           ->SetName("Crystal");
    }

    void CreateGadget() {
        SetIntegrator(new MDIntegrator);
        AddForce(new TForceEAM("platinum"));
        AddConditioner(new VerletList);
        AddDetector(new SysMonitor("md.out"));  
        AddDetector(new ThermoDetector(0.02));
    }

    void SystemSetting() {
        PBoundary=NONPERIODIC;
        MaxTime=0.2;
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
