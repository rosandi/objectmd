#ifndef _PARTICLE_PATH_HPP_
#define _PARTICLE_PATH_HPP_

#include <vector>
#include <fstream>
#include <omd/conditioner.h>
#include <omd/treader.h>
using namespace omd;


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
 * Parameters:
 *  - (name).file pathfilename : path file in x y z form
 *  - (name).index idx : index of atom whose cordinates defined by path (default=0)
 *  - (name).loop : loop the path on eof (default=false)
 *  - (name).terminate : terminate simulation on eof (default=true)
 *
 *  The name of the class by default is "Path". if no new name is given, this is the name used in parameter file. 
*/

class ParticlePath: public CalcConditioner {	
	int index, npath, current;
	vector<double> x, y, z;
	bool loop, finished, terminate;
	string pathfile;
	
public:		
	void readfile(string tablefile) {
		mdassert(file_exist(pathfile), "mising path data file. parameter: "+lower_case(get_name()));
		ifstream fl(tablefile.c_str());
		mdassert(!fl.fail(), "Can not open file"+tablefile);
		npath=0;
		do {
			double a, b, c;
			fl >> a >> b >> c >> std::ws;
			if(fl.bad()) die("bad data in "+tablefile);
			x.push_back(a);
			y.push_back(b);
			z.push_back(c);
			npath++;
		} while(!fl.eof());
	}
	
	ParticlePath(string TargetAtom) {
		set_name("Path");
		register_class("particle path");
		TargetName=TargetAtom;
		index=0;
		current=0;
		loop=false;
		finished=false;
		terminate=true;
		pathfile="-";
	}
	
	void ReadParameter() {
		string tag=replace_char(lower_case(get_name()),' ','_');
		SysParam->peek(mytag("file"),pathfile);
		SysParam->peek(mytag("loop"),loop);
		SysParam->peek(mytag("index"),index);
		SysParam->peek(mytag("terminate"),terminate);
		readfile(pathfile);
	}

	void PreCalculation(){
		if(finished) return;
		Target->Atoms(index).x=x[current];
		Target->Atoms(index).y=y[current];
		Target->Atoms(index).z=z[current];
		Target->Atoms(index).vx=0.0;
		Target->Atoms(index).vy=0.0;
		Target->Atoms(index).vz=0.0;
		current++;
		if(current==npath) {
			if(loop) current=0;
			else {
				if(terminate) System->Stop();
				finished=true;
			}
		}
	}
};

#endif
