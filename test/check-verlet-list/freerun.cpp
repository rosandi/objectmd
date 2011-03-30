#include <omd/systemgrid.hpp>
#include <omd/container.hpp>
//#include <potential/sw.hpp>
//#include <conditioner/VerletList.hpp>
#include <conditioner/VerletListFull.hpp>
#include <detector/ThermoDetector.hpp>
#include <detector/SysMonitor.hpp>


class CheckForce: public ForceKernel {
    VerletListFull* Verlet;
	public:
	CheckForce() {CutRadius=5.0;}
	void Init(MDSystem* ws) {
	  ForceKernel::Init(ws);
	    Verlet=dynamic_cast<VerletListFull*>(System->GetIterator());
        assert(Verlet->type_of("verlet list full neighbor"), "die");
	             
	}
	
	void ComputeHalf(Atom& at, Atom& to) {
      std::cerr << "pair("<<at.nid<<","<<to.nid<<")\n";
	}
	
	void ComputeFull(Atom& at, Atom& to) {
      std::cerr << "full("<<at.nid<<","<<to.nid<<"):\n";
      int iat,ito,nidx,lstart,lend;
      NeighborList* nlist;
      Verlet->GetIterationVariables(iat,ito,nidx,nlist);
      for(int i=nidx+1;i<nlist->end;i++) {
        int ka=nlist->list[i];
        std::cerr  << " ang("<<at.nid<<","<<to.nid<<","<<Atoms(ka).nid<<")\n";
      }
	}
};

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
        AddForce(new CheckForce);
//        AddForce(new StillingerWeber("silicon"));
        AddConditioner(new VerletListFull);
        AddDetector(new SysMonitor);  
        AddDetector(new ThermoDetector);
    }
    
    void SystemSetting() {
//      SetTemperature(param.double_value("temperature"));
    }
    
    void BeforeRun() {
        DumpAtoms("init.dat", WM_VELOCITY);
        PrintInfo("info.out");
        die("check...");
    }
    
};

int main(int argc, char* argv[]) {
    MyMDClass TheSim;
    TheSim.SetArgument(argc,argv);
    return TheSim.Run();
}
