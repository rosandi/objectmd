//-------Temperature-Pressure-Density--Detector-----------//

#ifndef _LOCAL_ORDER_PARAMETER_HPP_
#define _LOCAL_ORDER_PARAMETER_HPP_

#include <omd/dataslot.h>
#include <omd/iterator.h>
#include <detector/DataDumper.h>

using namespace omd;

/**
  @ingroup detector
  @brief Local Order Parameter Detector
  
  ref: Ivanov, D. S. & Zhigilei, L. V. Phys. Rev. B, 2003, 68, 064114
 
**/

class LocalOrderParameter: public DataDumper {
	int ilop;
	double Q0;
	double lc;
	double cutsqr;
	int *ncord;
	double *sinsum, *cossum;
	enum {do_calc, do_average} turn;
	bool enumerate;

public:
	LocalOrderParameter(string fn, double lconst, double rcut, bool enumerate_fn=true):
	DataDumper("", 0.0, fn, 3){
		set_name("local order parameter");
		register_class(get_name());		
		lc=lconst;
		cutsqr=rcut*rcut;
		Q0=4.0*M_PI/lc;
		ncord=NULL;
		sinsum=cossum=NULL;
		enumerate=enumerate_fn;
	}
	
	virtual void Init(MDSystem* WorkSys) {
		DataDumper::Init(WorkSys);
		ilop=ClaimAuxVariable(true, "lop");
	}

	virtual ~LocalOrderParameter() {
		MemFree(ncord);
		MemFree(sinsum);
		MemFree(cossum);
	}

	string GetFilename(){
		if(enumerate) return (FilenamePrefix+Filename+FilenamePostfix);
		return (FilenamePrefix+Filename);
	}

	void calculate_lop(int at, int to) {
		double dx, dy, dz;
		double d=CalcSqrDistance(at, to, dx, dy, dz);
		if(d<=cutsqr) {
			ncord[at]++; ncord[to]++;
			double qcos=(cos(Q0*dx)+cos(Q0*dy)+cos(Q0*dz));
			qcos+=cos(Q0*(dx+dy))+cos(Q0*(dy+dz))+cos(Q0*(dx+dz));
			double qsin=(sin(Q0*dx)+sin(Q0*dy)+sin(Q0*dz));
			qsin+=sin(Q0*(dx+dy))+sin(Q0*(dy+dz))+sin(Q0*(dx+dz));
			cossum[at]+=qcos;
			cossum[to]+=qcos;
			sinsum[at]+=qsin;
			sinsum[to]+=qsin;
		}
	}

	void average_lop(int at, int to) {
		double d=CalcSqrDistance(at, to);
		if(d<=cutsqr) {
			ncord[at]++;
			ncord[to]++;
			cossum[at]+=sinsum[to];
			cossum[to]+=sinsum[at];
		}
	}

	void IterationNode(int at, int to) {
		if(turn==do_calc)
			calculate_lop(at,to);
		else
			average_lop(at,to);
	}

	void find_local_order() {
		int na=GetNAtom();
		MemRealloc(ncord, na*sizeof(int));
		MemRealloc(sinsum, na*sizeof(double));
		MemRealloc(cossum, na*sizeof(double));
		
		for(int i=0;i<na;i++) {
			Atoms(i).aux[ilop]=0.0;
			ncord[i]=0;
			cossum[i]=0.0;
			sinsum[i]=0.0;
		}
		
		turn=do_calc;
		Iterator->Iterate(this);
		
		// use sinsum&cossum as buffer
		for(int i=0;i<na;i++) {
			double a=(cossum[i]*cossum[i]-sinsum[i]*sinsum[i]);
			double nn=(double)(36*ncord[i]*ncord[i]);
			sinsum[i]=cossum[i]=(nn>0.0)?a/nn:0.0;
		}
		
		turn=do_average;
		for(int i=0;i<na;i++) ncord[i]=1;
		Iterator->Iterate(this);
		
		for(int i=0;i<na;i++) {
			Atoms(i).aux[ilop]=cossum[i]/(double)(ncord[i]);
		}

	}

	void Measure() {
		find_local_order();
		DataDumper::Measure();
	}

};

#endif
