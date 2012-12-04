#include <omd.h>

ofstream trjfl;

class MyMDClass:public MDSystem {

    void CreateSystem() {
        AddAtom(new FCC("111",30,20,10, "platinum"))
           ->Create()
           ->SetTemperature(100.0)
           ->SetName("Crystal");
    }

    void CreateGadget() {
        SetIntegrator(new MDIntegrator);
        AddForce(new TForceEAM("platinum"));
        AddModify(new VerletList);
        AddDetector(new SysMonitor("md.out"));  
        AddDetector(new ThermoDetector(0.02));
    }

    void SystemSetting() {
        PBoundary=NONPERIODIC;
        MaxTime=0.2;
        SetOutputDirectory("output");
    }
// ------------     
    void InlineFunction() {
        trjfl << Atoms(0).x << " " << Atoms(0).y << " " << Atoms(0).z << " "
              << Atoms(10).x << " " << Atoms(10).y << " " << Atoms(10).z
              << std::endl;
    }
//-----------------------

    void BeforeRun() {
        DumpAtoms("init.dat");
        PrintInfo("info.out");
    }

};

int main(int argc, char* argv[]) {
    MyMDClass TheSim;
    
    // buka file
    trjfl.open("traj.dat");
    TheSim.SetArgument(argc,argv);
    TheSim.Run();
    trjfl.close();
    return 0;
}
