/* app
 ************************************************
 * Object-MD
 * Object Oriented Molecular Dynamics Program
 * 
 * (c)Yudi Rosandi, 2005
 *
 * Version 2.0 (10.11.08)
 *
 * Project started on July 2005.
 *
 ************************************************
 *
 */

#include <omd.hpp>

class MySim:public MDSystemGrid {
  
	void CreateSystem() {
    if(param.exist("atom:")) {
      istringstream ss("atom:", "end");
      
      while(ss.good()) {
        string snm, crys;
        mdassert(ss>>snm>>crys,"invalid atom definition: "+ss.str());
        
        if(crys=="fcc") {
          string ori,mat;
          double xm,ym,zm;
          mdassert(ss>>ori>>xm>>ym>>zm>>mat, "invalid creation options: "+ss.str());
          AddAtom(new FCC(ori,xm,ym,zm,mat))->SetName(snm);
        }
        
        else die("not implemented yet");
      }
    }
	}
  
  void CreateGroup() {
    
    if(param.exist("group:")) {
      istringstream ss("group:","end");
      
      while(ss.good()) {
        string snm, alg;
        mdassert(ss>>snm>>alg,"invalid group options: "+ss.str());
        AtomGroup* AG=AddAtomGroup(snm);
        
        if(alg=="box") {
          double x0,x1,y0,y1,z0,z1;
          assert(ss>>x0>>y0>>z0>>x1>>y1>>z1, "invalid group box region");
          AG->SelectBox(x0,x1,y0,y1,z0,z1);
        }
        
        else if(alg=="ibox") {
          double x0,x1,y0,y1,z0,z1;
          assert(ss>>x0>>y0>>z0>>x1>>y1>>z1, "invalid group box region");
          AG->SelectBoxInverse(x0,x1,y0,y1,z0,z1);
        }
        
        else if(alg=="gt") {
          double x0,y0,z0;
          assert(ss>>x0>>y0>>z0, "invalid coordinate");
          AG->SelectGT(x0,y0,z0);
        }

        else if(alg=="ge") {
          double x0,y0,z0;
          assert(ss>>x0>>y0>>z0, "invalid coordinate");
          AG->SelectGE(x0,y0,z0);
        }        
        
        else if(alg=="lt") {
          double x0,y0,z0;
          assert(ss>>x0>>y0>>z0, "invalid coordinate");
          AG->SelectLT(x0,y0,z0);
        }  
        
        else if(alg=="le") {
          double x0,y0,z0;
          assert(ss>>x0>>y0>>z0, "invalid coordinate");
          AG->SelectLE(x0,y0,z0);
        }

      }
    }
  }
  
	void CreateGadget() {
    if(!param.exist("integrator")) SetIntegrator(new MDIntegrator);
    
		NonReflecting* noref=new NonReflecting;
		AddForce(new TForceEAM("aluminum"))->SetEvaluator(noref);
    
		AddConditioner(new VerletList);
		AddConditioner(noref);
    
		AddDetector(new SysMonitor("md.out"));
		AddDetector(new ThermoDetector(0.05));
	}
  
	void BeforeRun() {
		PrintInfo("info.out");
		DumpAtoms("init_data");
	}
  
};

int main(int argc, char* argv[]) {
	MySim TheSim;
	TheSim.SetArgument(argc,argv);
	return TheSim.Run();
}
