//---------------------RESTART FILE PERIODIC SAVER---------------------//

#ifndef _RESTART_SAVER_HPP_
#define _RESTART_SAVER_HPP_

#include <iostream>
#include <ctime>
#include <omd/detector.h>

namespace omd {


/**
 * @ingroup detector
 * @brief Save periodicaly binary restart file
 *
 * Parameters:
 * - save.filename (filename) : binary filename to save.
 * - save.every (time) (second|minute|hour|day) : saving period in real time.
 * - save.overwrite (yes|no) : overwrite the previously created file or enumerate binary files.
 *   The file sequence number will be appended to the filename.
 *
*/

class RestartSaver:public Detector, public ParallelGadget {
	double saving_period;
	double elapsed_time;
	double next_save;
	string time_unit;
	time_t current_time;
	time_t previous_time;
	bool overwrite;
	
public:
	RestartSaver(double speriod=-1.0, string tunit="second") {
		set_name("restart");
		saving_period=speriod;
		time_unit=tunit;
		overwrite=true;
	}
	
	void PrintInfo(ostream& ost) {
		double tp=saving_period;
		if(time_unit=="minute") tp/=60.0;
		if(time_unit=="hour") tp/=(60.0*60.0);
		if(time_unit=="day") tp/=(60.0*60.0*24.0);
		ost << "ID."<<id<<" "<<get_name()
		    <<" every="<<tp<<"("<<time_unit<<")"
		    <<" filename=" << GetFilename();
		if(overwrite) ost << "(overwrite)";
		else ost << "(no overwrite)";
		ost <<std::endl;
	}
  
  void ReadParameter() {
    SysParam->peek(mytag("every"),saving_period);
    SysParam->peek(mytag("unit"),time_unit);
  }
	
	void Init(MDSystem* WorkSys) {
		Detector::Init(WorkSys);
		ParallelGadget::Init(WorkSys);

		time(&current_time);
		previous_time=current_time;
		Filename="omd.restart";
		
		// omd-parameter has priority
		SysParam->peek("save.overwrite",overwrite);
		SysParam->peek("save.filename",Filename); // or use default restart filename

		if(SysParam->exist("save.every")) {
			saving_period=SysParam->double_value("save.every");
			time_unit=SysParam->string_value("save.every",1);
			if(!(time_unit=="minute"||time_unit=="hour")) time_unit="second";
			
		}

		mdassert(saving_period>0.0, "the saving period is undefined! (parameter save.every)");
		
		if(time_unit=="minute") saving_period*=60.0;
		if(time_unit=="hour") saving_period*=60.0*60.0;
		if(time_unit=="day") saving_period*=60.0*60.0*24.0;
		
		next_save=saving_period;
		
	}
	
	void Measure() {
		string binfile(GetFilename());
		if(!overwrite) binfile.append(string(".")+as_string(NCalls));
		System->SaveSimulation(binfile);
		blog("saved simulation restart file to "+GetFilename(), LOGINFO);
	}

	void Detect() {
		int a=0;
		if(GetRank()==MDROOT) {
			time(&current_time);
			double dt=difftime(current_time, previous_time);
			if(dt>saving_period) {
				a=1;
				previous_time=current_time;
				NCalls++;
			}
		}
		SyncProcesses();
		Broadcast(a);
		if(a==1) Measure();
	}
};

}

#endif
