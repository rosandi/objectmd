/* omdlib
 ********************************************************
 *
 * ObjectMD
 * Molecular Dynamics Class Library
 * Version 3.0 (2009,2011) 
 *
 * Released under:
 *       GNU GENERAL PUBLIC LICENSE Version 3
 *       Copyright (C) 2005  Yudi Rosandi
 *
 * Read license.txt file in the source's root directory.
 *
 ********************************************************
 * ObjectMD header file
 *
 * Base class of all detectors
 *
*/

#ifndef _DETECTOR_H_
#define _DETECTOR_H_

#include <omd/gadget.h>
#include <omd/integrator.h>
#include <omd/system.h>

namespace omd {

/**
 * @ingroup detector
 * @brief Detector Base Class
 * 
 * This class is the parent of all Object-MD detectors. The descendant have to
 * implement the Detector::Measure() (abstract) function to perform measurement. 
 * The call to this function is triggered periodicaly using a sampling time,
 * passed to the constructor. The default value is zero, which
 * means that measurement will be executed every time step. 
 * The sampling time is defined in the simulation time (not simulation step).
 * Normally, it has unit of pico-second.
 *
 * Detector class preserves variables for filenames, the prefix, filename,
 * and the postfix. The prefix can be used to redirect the directory, where the
 * file is written. The postfix can be used to give file numbering, etc.
 * Files are not opened by this class. It is left to the descendant to 
 * manage the physical file.
 * 
*/

class Detector: public MDGadget {
protected:

    OMD_FLOAT  TSample;
    OMD_FLOAT  NextSample;
    
	string FilenamePrefix;
	string Filename;
	string FilenamePostfix;
	
public:
    Detector(OMD_FLOAT tm=0.0){
    	TSample=tm;
    	Target=NULL;
    	TargetName="";
    	Filename="";
    	FilenamePrefix="";
    	FilenamePostfix="";
    }
    
    virtual ~Detector(){}

	/**
	 * Reserved for joining/linking with other detector. To enable it the
	 * descendant must reimplement the function.
	 * 
	 **/

	virtual Detector* Join(MDGadget* det){
		die("unimplemented Join() in "+get_name()+". This class can't join"); 
		return this;
	}
	
	void Init(MDSystem* WorkSys){
		MDGadget::Init(WorkSys);
		NextSample=TSample;
		if(TSample<WorkSys->Integrator->TimeStep)
			TSample=WorkSys->Integrator->TimeStep;
			
		if(NextSample<WorkSys->ElapsedTime)
			NextSample=WorkSys->ElapsedTime;
	}

	virtual string GetFilename(){
		return (FilenamePrefix+Filename+FilenamePostfix);
	}
	
	virtual void SetFilename(string stname){Filename=stname;}
	virtual void SetFilenamePostfix(string stname){FilenamePostfix=stname;}
	virtual void SetFilenamePrefix(string stname){FilenamePrefix=stname;}

	void PrintInfo(ostream& ost){
		ost<<"id."<<id<<" "<<get_name()<<" -- target: "<<Target->get_name()<<std::endl;
	}
             
    
    void SetSampleTime(OMD_FLOAT tm) {TSample=tm;}
            
	// measurement
    virtual void Measure()=0;

    // detection timing     
    virtual void Detect() {
    	if(System->GetMode()==STATIC_MODE) {
    		Measure();
			NCalls++;
		} else {
			if(OnTime(NextSample)){
				Measure();
				NextSample+=TSample;
				NCalls++;
			}
		}
	}

    OMD_FLOAT GetSampleTime(){return TSample;}
	virtual Detector* SetName(string newname){set_name(newname); return this;}

};
}

#endif
