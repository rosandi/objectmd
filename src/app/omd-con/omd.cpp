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

#include <omd.h>
#include <muParser.h>

using namespace std;

class MySim:public MDSystemGrid {
  
  mu::Parser parser;
  
  // math variables
  double xml,yml,zml;
  vector<string> eqvarname;
  double eqvar[100]; // allocate for only 100 user defined variables
  
  double eval(string equ) {
    equ=trim(equ);
    if(equ.empty()) return 0.0;
    
    if(equ[0]=='{') { // write equation between {}
      equ=remove_char(equ,"{ }");
      try {
        parser.SetExpr(equ);
        return parser.Eval();
      } catch (mu::Parser::exception_type &e){
        std::cout << e.GetMsg() << std::endl;
        throw "equation parsing error";        
      }
    }
    
    return as_double(equ);
    
  }
    
	void CreateSystem() {
    if(param.exist("create:")) {
      istringstream ss(param.raw_string("create:", "end"));
      string snm, crys;
      
      while(ss>>snm>>crys) {
        
        if(crys=="fcc") {
          string ori,mat;
          string xm,ym,zm;
          mdassert(ss>>ori>>xm>>ym>>zm>>mat, "invalid creation options: "+ss.str());
          blog("create: "+snm+" "+crys);
          AddAtom(new FCC(ori,eval(xm),eval(ym),eval(zm),mat))->SetName(snm);
        }
        
        else die("not implemented yet: "+crys);
      }
    }
	}
  
  void CreateGroup() {
    
    if(param.exist("group:")) {
      istringstream ss(param.raw_string("group:","end"));
      string snm, alg;
      
      while(ss>>snm>>alg) {
        AtomGroup* AG=AddAtomGroup(snm);
        
        if(alg=="box") {
          string x0,x1,y0,y1,z0,z1;
          mdassert(ss>>x0>>y0>>z0>>x1>>y1>>z1, "invalid group box region");
          AG->SelectBox(eval(x0),eval(x1),eval(y0),eval(y1),eval(z0),eval(z1));
        }
        
        else if(alg=="ibox") {
          string x0,x1,y0,y1,z0,z1;
          mdassert(ss>>x0>>y0>>z0>>x1>>y1>>z1, "invalid group box region");
          AG->SelectInverseBox(eval(x0),eval(x1),eval(y0),eval(y1),eval(z0),eval(z1));
        }
        
        else if(alg=="gt") {
          string x0,y0,z0;
          mdassert(ss>>x0>>y0>>z0, "invalid coordinate");
          AG->SelectGT(eval(x0),eval(y0),eval(z0));
        }

        else if(alg=="ge") {
          string x0,y0,z0;
          mdassert(ss>>x0>>y0>>z0, "invalid coordinate");
          AG->SelectGE(eval(x0),eval(y0),eval(z0));
        }        
        
        else if(alg=="lt") {
          string x0,y0,z0;
          mdassert(ss>>x0>>y0>>z0, "invalid coordinate");
          AG->SelectLT(eval(x0),eval(y0),eval(z0));
        }  
        
        else if(alg=="le") {
          string x0,y0,z0;
          mdassert(ss>>x0>>y0>>z0, "invalid coordinate");
          AG->SelectLE(eval(x0),eval(y0),eval(z0));
        }

      }
    }
  }
  
  void PostCreation() {
    if(param.exist("modify:")) {
      istringstream ss(param.raw_string("modify:","end"));
      string snm,op;
      
      while(ss>>snm>>op) {
        
        if(op=="temperature") {
          double tempe;
          mdassert(ss>>tempe, "temperature value required");
          SearchContainer(snm)->SetTemperature(tempe);
        }
        
        else if(op=="velocity") {
          string vx,vy,vz;
          mdassert(ss>>vx>>vy>>vz, "velocity vector required");
          SearchContainer(snm)->SetVelocity(eval(vx),eval(vy),eval(vz));
        }
        
        else if(op=="kinetic") {
          string kine;
          mdassert(ss>>kine, "kinetic energy required");
          SearchContainer(snm)->SetKineticEnergy(eval(kine));
          
        }
        
        else if(op=="shift") {
          string sx,sy,sz;
          mdassert(ss>>sx>>sy>>sz, "shift distances required");
          SearchContainer(snm)->Shift(eval(sx),eval(sy),eval(sz));          
        }
        
        else die("operation not implemented: "+op);
      }
      
    }
  }
  
	void CreateGadget() {
    vector<string> evaluator;
    vector<ForceKernel*> evalforce;
    
    SetIntegrator(new MDIntegrator);
    
    if(param.exist("interaction:")) {
      istringstream ss(param.raw_string("interaction:","end"));
      string a,b,sty,fty,evalname;
      
      while(ss>>a>>b>>sty) {
        ForceKernel* FF;
        
        int dot=sty.find('.');
        mdassert(dot!=sty.npos, "can not define interaction type: "+sty);
        
        fty=sty.substr(0,dot);
        sty.erase(0,dot+1);
        
        if(fty=="eam") FF=new TForceEAM(sty);
        else if(fty=="pair") FF=new TForcePair(sty);
        else if(fty=="sw") FF=new StillingerWeber(sty);
        
        else die("interaction type not implemented: "+fty);
        
        AddForce(FF,a,b);
        
        // check if this force wants to be evaluated
        if(ss>>evalname) {
          if(evalname=="eval") {
            ss>>evalname;
            evaluator.push_back(evalname);
            evalforce.push_back(FF);
          }
        }
        
      }
    }
    
    if(param.exist("cond:")) {
      istringstream ss(param.raw_string("cond:","end"));
      string dnm,cnm;
      
      while(ss>>dnm) {
        
        int dot=dnm.find('.');
        if(dot!=dnm.npos) cnm=dnm.substr(0,dot);
        else cnm=dnm;        
        
        if(cnm=="verlet") {
          if(dnm=="verlet.full") AddConditioner(new VerletListFull);
          else AddConditioner(new VerletList);
        }
        
        else if(cnm=="clamp") AddConditioner(new CoordClamp)->set_name(dnm);
        else if(cnm=="dyndt") AddConditioner(new DynamicTimeStep)->set_name(dnm);
        else if(cnm=="source") AddConditioner(new EnergySource)->set_name(dnm);
        else if(cnm=="damp") AddConditioner(new ForceDamper)->set_name(dnm);
        else if(cnm=="nrb") AddConditioner(new NonReflecting)->set_name(dnm);
        else die("conditioner not yet enabled: "+cnm);
        
      }
    }
        
    if(param.exist("detect:")) {
      istringstream ss(param.raw_string("detect:","end"));
      string dnm,cnm;
      
      while(ss>>dnm) {
        
        int dot=dnm.find('.');
        if(dot!=dnm.npos) cnm=dnm.substr(0,dot);
        else cnm=dnm;
        
        // name is set to sync with parameter tag
        if(cnm=="monitor") AddDetector(new SysMonitor)->set_name(dnm);
        else if(cnm=="thermo") AddDetector(new ThermoDetector)->set_name(dnm);
        
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
      string cmd;
      
      while(ss>>cmd) {
        
        if(cmd=="dump") {
          string cnm,fnm;
          ss>>cnm>>fnm;
          SearchContainer(cnm)->DumpAtoms(fnm);
        }
        
        else if(cmd=="info") {
          string fnm;
          mdassert(ss>>fnm, "filename required");
          PrintInfo(fnm);
        }
        
        else if(cmd=="quit") die("quit");
        
        
        else die("command not implemented: "+cmd);
      }
                       
    }
	}
  
  void Init() {
    
    // variable initialization
    // if already exist->reassign
    
    if(param.exist("var:")) {
      istringstream ss(param.raw_string("var:","end"));
      string vname,vval;
      int fnd;
      while(ss>>vname>>vval) {
        int iv=(int)eqvarname.size();
        bool found=false;
        
        for(fnd=0;fnd<iv;fnd++) {
          if(vname==eqvarname[fnd]) {
            found=true;
            break;
          }
        }
        
        if(found) eqvar[fnd]=eval(vval);
        else {  
          int idx=(int)eqvarname.size();
          eqvar[idx]=eval(vval);
          eqvarname.push_back(vname);
          parser.DefineVar(vname,eqvar+idx);
        }
        
      }
    }
  }
  
public:
  MySim(int &argc, char** &argv) {
    SetArgument(argc,argv);
    parser.DefineConst("Pi",M_PI);
    parser.DefineConst("kb",8.6173324E-5); // eV/K
    parser.DefineVar("x0",&Box.x0);
    parser.DefineVar("y0",&Box.y0);
    parser.DefineVar("z0",&Box.z0);
    parser.DefineVar("x1",&Box.x1);
    parser.DefineVar("y1",&Box.y1);
    parser.DefineVar("z1",&Box.z1);
  }
    
};

int main(int argc, char* argv[]) {
  // analize command line arguments
  try {
    ParamHandler p(argc,argv);
  
    if(p.exist("--param")) {
      mdassert(file_exist(p.string_value("--param")), "parameter file not found");
    } else {
      mdassert(file_exist("omd-parameter"), "default parameter file not found");
    }
    
  } catch (const char* st) {
    std::cerr<<st<<std::endl;
    return -1;
  } catch (...) {
    std::cerr<<"undefined error\n";
    return -1;
  }
  
	MySim TheSim(argc,argv);
	return TheSim.Run();
}
