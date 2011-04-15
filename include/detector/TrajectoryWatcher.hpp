#ifndef _TRAJECTORY_WATCHER_HPP_
#define _TRAJECTORY_WATCHER_HPP_

#include <cstring>
#include <omd/detector.hpp>
#include <omd/paragadget.hpp>

/**
 * @ingroup detector
 * @brief Watches one atom trajectory
 */

class TrajectoryWatcher: public Detector, public ParallelGadget {

	ofstream fl;
	string filename;
	int gid;
	int watchflag;
	char st[256];
	
	public:
		
		TrajectoryWatcher(string TargetAtom, string fname) {
			set_name("TRAJECTORY WATCHER");
			register_class(get_name());
			TargetName=TargetAtom; gid=-1; filename=fname;
		}

		TrajectoryWatcher(int x_id, string fname) {
			gid=x_id; filename=fname;
		}
		
		~TrajectoryWatcher() {fl.close();}

		void Init(MDSystem* WorkSys) {
			Detector::Init(WorkSys);
			ParallelGadget::Init(WorkSys);
			
			if(GetRank()==0)
				fl.open(filename.c_str(), ios::trunc);
			else
				fl.open(filename.c_str(), ios::app);
			
			assert(fl.good(), "can not open file for reading "+filename);
	
			// only one atom!
			watchflag=ClaimFlagBit();
			
			if(GetNAtom()) {
				if(gid>=0) {
					for(int i=0; i<GetNAtom(); i++)
						if(Atoms(i).gid==gid){SetFlag(i, watchflag);break;}
				} else {SetFlag(0, watchflag);}
			}

		}

		void Measure() {
			// may be the container is empty.....
			int na=GetNAtom();
			int idx=0;

			if(na>1) {
				for(int i=0;i<na;i++) if(CheckFlag(i, watchflag)){idx=i;break;}
			}

			if(na){
				Atom* a=AtomPtr(idx);
				if((a->flag&watchflag)&&(a->flag&FLAG_ACTIVE)){					
					sprintf(st, "%0.8f %0.8f %0.8f %0.8f %0.8f %0.8f %0.8f\n",
							System->ElapsedTime, a->x, a->y, a->z, a->vx, a->vy, a->vz);
					fl << st; fl.flush();
				}
			}
		}
		
		virtual void Detect() {
			Measure();
		}

};

#endif
