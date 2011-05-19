#include <omd/systemgrid.hpp>
#include <omd/container.hpp>
#include <potential/sw.hpp>
#include <crystal/Diamond.hpp>
#include <conditioner/VerletListFull.hpp>
#include <detector/ThermoDetector.hpp>
#include <detector/SysMonitor.hpp>

class MyMDClass:public MDSystemGrid {

    void CreateSystem() {
        AddAtom(new Diamond("100",
                param.double_value("monolayer",0),
                param.double_value("monolayer",1),
                param.double_value("monolayer",2),
                "silicon"))                
                ->SetName("TARGET");
    }
    
    void CreateGroup() {
        int lc=param.double_value("silicon.lattice_constant");
        
        AddAtomGroup("fixed")
        ->SelectLT(Box.x0+lc/8.0, "all")
        ->SelectGT(Box.x1-lc/8.0, "all")
        ->SelectLT(Box.y0+lc/8.0, "all")
        ->SelectGT(Box.y1-lc/8.0, "all")
        ->SelectLT(Box.z0+lc/8.0, "all");
        
        AddAtomGroup("damped")
        ->SelectLT(Box.x0+lc, "all")
        ->SelectGT(Box.x1-lc, "all")
        ->SelectLT(Box.y0+lc, "all")
        ->SelectGT(Box.y1-lc, "all")
        ->SelectLT(Box.z0+lc, "all");        

    }

    void CreateGadget() {
        SetIntegrator(new MDIntegrator);
        AddForce(new StillingerWeber("silicon"));
        AddConditioner(new VerletListFull);
        AddConditioner(new CoordClamp("fixed", "xyz"));
        AddConditioner(new ForceDamper("damped"));
        AddDetector(new SysMonitor);  
        AddDetector(new ThermoDetector);
    }
    
    void SystemSetting() {
        SetTemperature(param.double_value("temperature"));
        SetBoundaryCondition(PERIODIC_X|PERIODICY);
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
