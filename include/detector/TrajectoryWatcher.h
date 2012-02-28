#ifndef _TRAJECTORY_WATCHER_HPP_
#define _TRAJECTORY_WATCHER_HPP_

#include <cstring>
#include <vector>
#include <omd/detector.h>
#include <omd/paragadget.h>

using std::vector;

/**
 * @ingroup detector
 * @brief Watches one atom trajectory
 *
 * Parameters:
 *   - watch.file filename : sets trajectory filename
 *   - watch.index idx1+idx2+... : define atom indices to watch separated by + (plus)
 *
 * The tag "watch" may be altered using set_tag().
 */

class TrajectoryWatcher: public Detector, public ParallelGadget {

	ofstream fl;
	string filename;
	char st[4096];
	vector<int> vindex;
	
public:
	
	TrajectoryWatcher(string TargetAtom, string fname) {
		set_name("TRAJECTORY WATCHER");
		SetTag("watch");
		register_class(get_name());
		TargetName=TargetAtom;
		filename=fname;
	}
	
	~TrajectoryWatcher() {fl.close();}
	
	void ReadParameter() {
		SysParam->peek(mytag("file"),filename);
		if(SysParam->exist(mytag("index"))) {
			size_t na=GetNAtom();
			std::istringstream ist(replace_char(param.string_value(mytag("index")), '+', ' '));
			while(ist.good()) {
				int i;
				ist >> i;
				if(i<na) vindex.push_back(i);
			}
		}
	}

	void Init(MDSystem* WorkSys) {
		Detector::Init(WorkSys);
		ParallelGadget::Init(WorkSys);
		
		if(GetRank()==0)
			fl.open(filename.c_str(), std::ios::trunc);
		else
			fl.open(filename.c_str(), std::ios::app);
		
		mdassert(fl.good(), "can not open file for reading "+filename);
		
		if(vindex.size()==0) vindex.push_back(0);
		for(size_t i=0; i<vindex.size();i++) {
			if(vindex[i]>=Target->GetNAtom()) 
				die("index out of range: "+as_string(vindex[i])+
			        " "+Target->get_name()+" has "+as_string(Target->GetNAtom())+" atoms");	
		}
	}

	void Measure() {
		// may be the container is empty.....
		// consider parallel omd!
		if(GetNAtom()==0) return;
		
		for(size_t i=0;i<vindex.size();i++) {
			Atom* a=AtomPtr(vindex[i]);
			if (a->flag&FLAG_ACTIVE) {					
				sprintf(st, "%0.8f %0.8f %0.8f %0.8f %0.8f %0.8f %0.8f ",
						System->ElapsedTime, a->x, a->y, a->z, a->vx, a->vy, a->vz);
			}
		}
		
		fl << st << std::endl;
		fl.flush();
	}
	
	virtual void Detect() {Measure();}

};

#endif
