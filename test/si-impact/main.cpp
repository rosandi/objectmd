#include <omd/systemgrid.hpp>
#include <omd/container.hpp>
#include <potential/sw.hpp>
#include <potential/tpair.hpp>
#include <crystal/Diamond.hpp>
#include <crystal/Projectile.hpp>
#include <conditioner/VerletListFull.hpp>
#include <conditioner/CoordClamp.hpp>
#include <conditioner/ForceDamper.hpp>
#include <conditioner/DynamicTimeStep.hpp>
#include <conditioner/RangeLimiter.hpp>
#include <detector/ThermoDetector.hpp>
#include <detector/SysMonitor.hpp>
#include <detector/TrajectoryWatcher.hpp>

class MyMDClass:public MDSystemGrid {

    void CreateSystem() {
        MDSystemGrid::CreateSystem();
        
        Projectile pa("argon",
          param.double_value("projectile",0)+
          param.double_value("projectile.offset",0),
          param.double_value("projectile",1)+
          param.double_value("projectile.offset",1),
          param.double_value("projectile",2),
          param.double_value("projectile.theta"),
          param.double_value("projectile.phi"),
          param.double_value("projectile.energy"));
          
        Atom* a=pa.AtomPtr(0);
        Atom* p=SearchContainer("projectile")->AtomPtr(0);
        
        p->x=a->x;
        p->y=a->y;
        p->z=a->z;
        p->vx=a->vx;
        p->vy=a->vy;
        p->vz=a->vz;
        
        blog("projectile: x("+
              as_string(p->x)+","+
              as_string(p->y)+","+
              as_string(p->z)+") v("+
              as_string(p->vx)+","+
              as_string(p->vy)+","+
              as_string(p->vz)+")");
        
    }

    void CreateGroup() {
    
        double bt=param.double_value("border.fix");
        double bd=param.double_value("border.damp");

        AddAtomGroup("border")
          ->SelectType(0)
          ->SelectInverseBox(Box.x0+bt, Box.x1-bt, 
                             Box.y0+bt, Box.y1-bt,
                             Box.z0+bt, INFINITE);
                             
       AddAtomGroup("damp-area")
          ->SelectType(0)
          ->SelectInverseBox(Box.x0+bd, Box.x1-bd, 
                             Box.y0+bd, Box.y1-bd,
                             Box.z0+bd, INFINITE);
    }
    
    void CreateGadget() {
        SetIntegrator(new MDIntegrator);
        AddForce(new StillingerWeber("silicon"), "target", "target");
        AddForce(new TForcePair("ar-si.zbl"), "projectile", "target");
        AddConditioner(new VerletListFull);
        AddConditioner(new ForceDamper("damp-area"));
        AddConditioner(new CoordClamp("border"));
        AddConditioner(new DynamicTimeStep);
        AddConditioner(new RangeLimiter("projectile", 1000.0));
        AddDetector(new SysMonitor);
        AddDetector(new ThermoDetector);
        AddDetector(new TrajectoryWatcher("projectile", "proj.trj"));
    }

    void BeforeRun() {
        PrintInfo("info.out");
        DumpAtoms("init.dat");
    }
                            
    void AfterRun() {
      DumpAtoms("final.dat", WM_VELOCITY|WM_NID|WM_TID);
      SaveSimulation("final.bin");
    }

};

int main(int argc, char* argv[]) {
    MyMDClass TheSim;
    TheSim.SetArgument(argc,argv);
    return TheSim.Run();
}
