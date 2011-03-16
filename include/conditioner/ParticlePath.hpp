
#include <vector>
#include <fstream>
#include <omd/conditioner.hpp>
#include <omd/tablereader.hpp>

using std::vector;
using std::ifstream;

/**
 * @ingroup conditioner
 * @brief Moves a particle in a defined path
 *
 * This class reads particle path from a file, contains x,y,z series of
 * coordinate. One coordinate entry is used in every single MD step.
 * Circulated.
 * 
*/

class ParticlePath: public Calc_Conditioner {	
	int Index, NPath, CurrentPath, Interrupt;
	OMD_FLOAT *x, *y, *z;
	bool Circulate;
	
	public:		
		void ReadFile(const char* tablefile) {
			vector<OMD_FLOAT> px;
			vector<OMD_FLOAT> py;
			vector<OMD_FLOAT> pz;
			ifstream fl(tablefile);
			if(fl.fail()) THROWSTR("Can not open file", tablefile);
			
			do {
				OMD_FLOAT a, b, c;
				fl >> a >> b >> c >> std::ws;
				px.push_back(a);
				py.push_back(b);
				pz.push_back(c);
			} while(!fl.eof());
			
			NPath=px.size();
			x=new OMD_FLOAT[NPath];
			y=new OMD_FLOAT[NPath];
			z=new OMD_FLOAT[NPath];
			for(int i=0;i<NPath;i++) {
				x[i]=px[i];
				y[i]=py[i];
				z[i]=pz[i];
			}
		}
		
		ParticlePath(const char* TargetAtom, const char* tablefile, bool circ=true, int idx=0) {
			TargetName=TargetAtom;
			Index=idx;
			set_name("Particle Path");
			ReadFile(tablefile);
			CurrentPath=-1;
			Circulate=circ;
		}

		void PreCalculation(){
			if(CurrentPath==-1){CurrentPath=0; return;}
			Target->Atoms[Index].x=x[CurrentPath];
			Target->Atoms[Index].y=y[CurrentPath];
			Target->Atoms[Index].z=z[CurrentPath];
			Target->Atoms[Index].vx=0.0;
			Target->Atoms[Index].vy=0.0;
			Target->Atoms[Index].vz=0.0;
			CurrentPath++;
			if(CurrentPath==NPath) {
				if(Circulate) CurrentPath=0;
				else System->ElapsedTime=System->MaxTime;
			}
		}
};
