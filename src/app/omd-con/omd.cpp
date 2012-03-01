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
      istringstream ss(param.raw_string("atom:", "end"));
      
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
      istringstream ss(param.raw_string("group:","end"));
      
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
          AG->SelectInverseBox(x0,x1,y0,y1,z0,z1);
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
  
  void PostCreation() {
    if(param.exist("modify:")) {
      istringstream ss(param.raw_string("modify:","end"));
      
      while(ss.good()) {
        string snm,op;
        mdassert(ss>>snm>>op,"invalid modify options: "+ss.str());
        
        if(op=="temperature") {
          double tempe;
          mdassert(ss>>tempe, "temperature value required");
          SearchContainer(snm)->SetTemperature(tempe);
        }
        
        else if(op=="velocity") {
          double vx,xy,xz;
          mdassert(ss>>vx>>vy>>vz, "velocity vector required");
          SearchContainer(snm)->SetVelocity(vx,xy,xz);
        }
        
        else if(op=="kinetic") {
          double kine;
          mdassert(ss>>tempe, "kinetic energy required");
          SearchContainer(snm)->SetKinetic(kine);
          
        }
        
        else if(op=="shift") {
          double sx,sy,sz;
          mdassert(ss>>sx>>sy>>sz, "shift distances required");
          SearchContainer(snm)->Shift(sx,sy,sz);          
        }
        
        else die("operation not implemented: "+op);
      }
      
    }
  }
  
	void CreateGadget() {
    vector evaluator;
    vector evalforce;
    
    SetIntegrator(new MDIntegrator);
    
    if(param.exist("interaction:")) {
      istringstream ss(param.raw_string("interaction:"));
      
      while(ss.good()) {
        string a,b,sty,fty,eval;
        ForceKernel* FF;
        
        mdassert(ss>>a>>b>>sty, "invalid interaction options: "+ss.str());
        int dot=sty.find('.');
        mdassert(dot!=sty.npos, "can not define interaction type: "+sty);
        
        fty=sty.substr(0,dot);
        sty.erase(0,dot+1);
        
        if(fty=="eam") FF==new TForceEAM(sty);
        else if(fty=="pair") FF=new TForcePair(sty);
        else if(fty=="sw") FF=new StillingerWeber(sty);
        
        else die("interaction type not implemented: "+fty);
        
        AddForce(FF,a,b);
        
        // check if this force wants to be evaluated
        if(ss>>eval) {
          if(eval=="eval") {
            ss>>eval;
            evaluator.push_back(eval);
            evalforce.push_back(FF);
          }
        }
        
      }
    }
    
    if(param.exist("cond:")) {
      istringstream ss(param.raw_string("cond:","end"));
      
      while(ss.good()) {
        string cnm;
        mdassert(ss>>cnm, "invalid conditioner option");
        
        if(cnm=="verlet") AddConditioner(new VerletList);
        else if(cnm=="verlet.full") AddConditioner(new VerletListFull);
        else if(cnm=="clamp") AddConditioner(new CoordClamp);
        else if(cnm=="dyndt") AddConditioner(new DynamicTimeStep);
        else if(cnm=="source") AddConditioner(new EnergySource);
        else if(cnm=="damp") AddConditioner(new ForceDamper);
        else if(cnm=="boundary.nrb") AddConditioner(new NonReflecting);
        else die("conditioner not yet enabled: "+cnm);
        
      }
    }
        
    if(param.exist("detect:")) {
      istringstream ss(param.raw_string("detect:","end"));
      
      while(ss.good()) {
        string dnm;
        mdassert(ss>>dnm, "invalid detector options: "+ss.str());
        
        if(cnm=="monitor") AddDetector(new SysMonitor);
        else if(cnm=="thermo") AddDetector(new ThermoDetector);
        
        else die("detector not yet enabled: "+cnm);
        
      }

    }
    
    // set all evaluators
    for(int i=0;i<(int)evalforce.size();i++) {
      evalforce[i]->SetEvaluator(SearchGadget(evaluator[i]));
    }
    
    
	}
  
	void BeforeRun() {
    if(param.exist("prerun:")) {
      istringstream ss(param.raw_string("prerun:","end"));
      
      while(ss.good()) {
        string cmd;
        if(!(ss>>cmd)) break;
        
        if(cmd=="dump") {
          string cnm,fnm;
          ss>>cnm>>fnm;
          SearchContainer(cnm)->DumpAtoms(fnm);
        }
        
        else if(cmd=="info") {
          string fnm;
          assert(ss>>fnm, "filename required");
          PrintInfo(fnm);
        }
        
        else die("command not implemented: "+cmd);
      }
                       
    }
	}
  
};

int main(int argc, char* argv[]) {
	MySim TheSim;
	TheSim.SetArgument(argc,argv);
	return TheSim.Run();
}
