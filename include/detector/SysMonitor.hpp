
#ifndef _SYS_MONITOR_
#define _SYS_MONITOR_

#include <string>
#include <omd/detector.hpp>
#include <omd/dataslot.hpp>
#include <omd/paragadget.hpp>

using std::ios;
using std::string;

/**
  @ingroup detector
  @brief The System Monitor
  
  This class provides the simulation data to be printed to the
  screen and a file. The data are:
    
    - the step
    - the simulation elapsed time
    - the time step
    - the kinetic, potential, and total energy 

 Parameters:
    - @e monitor.filename (file name)
    - @e monitor.nofile : suppress creating output file
    - @e monitor.sample (step)  : write data every nstep steps (sample time). 
      Used also by MDSystem.

**/

class SysMonitor:public Detector, public ParallelGadget {

	ofstream fout;
	bool     usefile;
	bool     usescreen;
	bool     firsttime;
	int      every;
	
	OMD_FLOAT   RelPotential;
	OMD_FLOAT   MechanicEnergy;
	
public:

    SysMonitor(string fn="md.out", bool oneline=false) {
                  	
		set_name("SYSTEM_MONITOR_GRID");
		register_class("SYSTEM_MONITOR");
		usefile=true;
		every=1;
		SetFilename(fn);
		firsttime=true;
	}
	        
    virtual ~SysMonitor() {
    	if(GetRank()==0&&usefile) fout.close();
	}
	
	void ReadParameter() {
		SysParam->peek("monitor.filename", Filename);
	    if(SysParam->exist("monitor.nofile")) usefile=false;
		SysParam->peek("monitor.sample", every);
		if(every<1) every=1;
	}
    
    void Init(MDSystem* WorkSys) {
    	Detector::Init(WorkSys);
    	ParallelGadget::Init(WorkSys);

	    if(GetRank()==ROOT) {
	    	if(usefile)fout.open(GetFilename().c_str(), ios::trunc);

	    	RegisterMessageSlot(new DataSlot("@",0))
	    		->SetFormat("%d")
				->SetData(System->Step);

	    	RegisterMessageSlot(new DataSlot("t",0))
	    		->SetFormat("%0.3E")
				->SetData(System->ElapsedTime);
			
	    	RegisterMessageSlot(new DataSlot("dt",0))
	    		->SetFormat("%0.3E")
				->SetData(System->Integrator->TimeStep);
	    	
			RegisterMessageSlot(new DataSlot("Ek",0))
	    		->SetFormat("%0.5E")
	    		->SetData(System->Kinetic);

	    	RegisterMessageSlot(new DataSlot("Ep",0))
	    		->SetFormat("%0.5E")
	    		->SetData(RelPotential);

	    	RegisterMessageSlot(new DataSlot("E",0))
	    		->SetFormat("%0.5E")
	    		->SetData(MechanicEnergy);

		}

    }

	virtual void Measure() {
		if(GetRank()==ROOT){
			RelPotential=System->Potential-System->BasePotential;
			MechanicEnergy=System->Kinetic+System->Potential-System->BasePotential;
			
			if (usefile) {
				int nslot=System->MessageSlots.size();
				if(firsttime) {
					fout << "#";
					for (int i=0;i<nslot;i++)
						fout << System->MessageSlots[i]->GetLabel() << " ";
					fout << '\n';
					firsttime=false;
				}
	    		for(int i=0;i<nslot;i++)
	    			fout << System->MessageSlots[i]->AsString() << " ";
	    			
    			fout << '\n';
				fout.flush();
				
			}
		}
	}
	
	virtual void Detect() {
		if(!(System->Step%every)) Measure();
	}
	
};

#endif
