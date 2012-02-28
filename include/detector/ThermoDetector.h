//-------Temperature-Pressure-Density--Detector-----------//

#ifndef _THERMO_DETECTOR_HPP_
#define _THERMO_DETECTOR_HPP_

#include <omd/dataslot.h>
#include <omd/iterator.h>
#include <omd/paragadget.h>
#include <detector/DataDumper.h>

/**
  @ingroup detector
  @brief Temperature and Pressure Detector
  
  This detector calculates the temperature, pressure and density
  of the atoms. These are local quantities, measured in the neighbor
  volume of size \f$\frac{4}{3}\pi R_{cut}^3\f$, where \f$R_{cut}\f$
  is the cut radius, taken from the interaction potential.
  * 
  * The system temperature is measured with respect to the central of
  * mass velocity. If the system has a drift velocity, the measured
  * temperature will be lower than the one calculated using
  * the system kinetic energy.
  * 
  * Parameters: see DataDumper
  *
**/

class ThermoDetector: public DataDumper, public ParallelGadget {

	double *sumas,*sumek,*sumvir;
	double *sumvx,*sumvy,*sumvz;
	double ek,masat,masto,ercut,volcut;
	double system_temperature,system_pressure; // System's temperature and pressure
	double avg_temp; // average of atom temperature
	double avg_pres; // average of atom pressure
	int    *nneig;
	int    nalloc;
	int    tidx,pidx,nidx;
	bool   intensive_mode; // measure every timestep
	enum {tempe,press} process;

public:
	ThermoDetector(double tm=-1.0, double rcut=-1.0, string fn="Data", int extlen=3)
	:DataDumper("", tm, fn, extlen){
		set_name("THERMO DETECTOR");
		register_class(get_name());
		system_temperature=0.0;
		avg_temp=0.0;
		system_pressure=0.0;
		intensive_mode=false;
		sumas=sumek=sumvx=sumvy=sumvz=NULL;
		nneig=NULL;
		ercut=rcut;
	}
	
	void SetIntensive(bool inten=true){
		intensive_mode=inten;
		if(inten) blog("intensive mode: true");
		else blog("intensive mode: false");
	}
	
	void SetTemperature(int idx, double t) {Atoms(idx).aux[tidx]=t;}
	double GetTemperature(int idx) {return Atoms(idx).aux[tidx];}
	double GetTemperature(){return system_temperature;}
	
	/**
	 * read the average temperature of the atoms. monomers and dimers are
	 * excluded.
	 */

	double GetTemperatureAvg(){
		mdassert(intensive_mode, "attempt to invoke GetTemperatureAvg in non-intensive mode!");
		return avg_temp;
	}
	
	double GetPressure(int idx) {return Atoms(idx).aux[pidx];}
	double GetPressure(){return system_pressure;}
	double GetCutRadius(){return ercut;}
	void GetCMVelocities(double* &cvx, double* &cvy, double* &cvz){
		cvx=sumvx;cvy=sumvy;cvz=sumvz;
	}
	
	int GetAuxTempIndex(){return tidx;}
	int GetAuxPresIndex(){return pidx;}
	int GetAuxDensIndex(){return nidx;}
	double GetDetectVolume(){return volcut;}

	void Measure() {		
		Target->PushInfo("$ Temperature "+
		                 as_string(system_temperature, "%0.5e")+" K (kinetic)");
		Target->PushInfo("$ Pressure "+as_string(system_pressure,"%0.5e")+" GPa");
		DataDumper::Measure();
		Target->PopInfo(2);
	}

	void Init(MDSystem* WorkSys) {
		DataDumper::Init(WorkSys);
		ParallelGadget::Init(WorkSys);

		if(GetRank()==0) {
			RegisterMessageSlot(new DataSlot("temp"))
				->SetFormat("%0.3E")->SetData(system_temperature);
			RegisterMessageSlot(new DataSlot("avtemp"))
				->SetFormat("%0.3E")->SetData(avg_temp);
			RegisterMessageSlot(new DataSlot("pres"))
				->SetFormat("%0.3E")->SetData(system_pressure);
			RegisterMessageSlot(new DataSlot("avpres"))
				->SetFormat("%0.3E")->SetData(avg_pres);
		}

		nalloc=GetNAtom();
		tidx=ClaimAuxVariable(true,"temp");
		pidx=ClaimAuxVariable(true,"press");
		nidx=ClaimAuxVariable(true,"dens");
		MemAlloc(sumas,sizeof(double)*nalloc);
		MemAlloc(sumek,sizeof(double)*nalloc);
		MemAlloc(sumvx,sizeof(double)*nalloc);
		MemAlloc(sumvy,sizeof(double)*nalloc);
		MemAlloc(sumvz,sizeof(double)*nalloc);
		MemAlloc(sumvir,sizeof(double)*nalloc);
		MemAlloc(nneig,sizeof(int)*nalloc);
		if(ercut<0.0) ercut=System->GetMaxCutRadius();
		mdassert(ercut>0.0, "cut radius is zero... static mode?");
		volcut=4.0*M_PI*ercut*ercut*ercut/3.0;
		
		System->PushInfo("$ DetectorVolume "+as_string(volcut)+" Ang^3");
		if(System->GetMode()==STATIC_MODE) SetIntensive(true);
	}

	virtual ~ThermoDetector() {
		MemFree(sumas);
		MemFree(sumek);
		MemFree(sumvx);
		MemFree(sumvy);
		MemFree(sumvz);
		MemFree(nneig);
	}

	virtual bool CheckLoopTemp(int idx) {
		Atom* a=AtomPtr(idx);
		nneig[idx]++;
		masat=GetMass(a);
		sumas[idx]+=masat;
		sumvx[idx]+= masat*a->vx;
		sumvy[idx]+= masat*a->vy;
		sumvz[idx]+= masat*a->vz;			
		ek=0.5*masat*(a->vx*a->vx+a->vy*a->vy+a->vz*a->vz);
		sumek[idx]+=ek;
		return true;
	}

	virtual void IterationNodeTemp(int at, int to) {
		double d=CalcDistance(at,to);
		Atom* a=AtomPtr(at);
		Atom* b=AtomPtr(to);
		if (d<=ercut) {
			masto=GetMass(b);	
			nneig[to]++;
			nneig[at]++;
			sumas[to]+=masat;
			sumas[at]+=masto;
			sumvx[to]+=masat*a->vx;
			sumvy[to]+=masat*a->vy;
			sumvz[to]+=masat*a->vz;                
			sumvx[at]+=masto*b->vx;
			sumvy[at]+=masto*b->vy;
			sumvz[at]+=masto*b->vz;                
			sumek[to]+=ek;
			sumek[at]+=0.5*masto*(b->vx*b->vx+b->vy*b->vy+b->vz*b->vz);
		}
	}

	virtual void ExtractTemperature() {
		process=tempe;
		Iterator->IterateHalf(this);
		int na=GetNAtom();
		int navg=0;
		avg_temp=0.0;
		
		for (int i=0;i<na;i++) {
			Atom* a=AtomPtr(i);
			sumvx[i]/=sumas[i]; 
			sumvy[i]/=sumas[i];
			sumvz[i]/=sumas[i];
			sumek[i]-=0.5*sumas[i]*
				      (sumvx[i]*sumvx[i]+sumvy[i]*sumvy[i]+sumvz[i]*sumvz[i]);				      
			a->aux[tidx]=Unit->Temperature(sumek[i]/(double)nneig[i]);
			a->aux[nidx]=(double)nneig[i]/volcut;

			if(nneig[i]<2) continue;
			if(a->flag&FLAG_GHOST) continue;
			
			if(a->flag&FLAG_ACTIVE) {
				avg_temp+=Atoms(i).aux[tidx];
				navg++;
			}

		}

		avg_temp=TakeSUM(avg_temp);
		navg=TakeSUM(navg);
		if(navg>0) avg_temp/=(double)navg;
		else avg_temp=0.0;
		
	}

	virtual bool CheckLoopPressure(int idx) {
		sumvir[idx]+=Atoms(idx).virial;  // my self
		return true;
	}

	virtual void IterationNodePressure(int at, int to) {
		double d=CalcDistance(at,to);
		if(d<=ercut){
			sumvir[at]+=Atoms(to).virial;
			sumvir[to]+=Atoms(at).virial;
		}
    // from this we get: N.<sum f_ij.r_ij> = sum_k (f_ij.r_ij)_k 
	}
	
	virtual void ExtractPressure() {
		process=press;
		Iterator->IterateHalf(this);
		int na=GetNAtom();
		int navg=0;
		avg_pres=0.0;
		for (int i=0;i<na;i++) {
			Atom* a=AtomPtr(i);
			a->aux[pidx]=(nneig[i]==0)?0.0:
      // this is <sum f_ij.r_ij> /3.omega no double sum in half-loop iterator
      // for full-loop consider reimplement ReturnForce@ForceKernel
			Unit->Pressure((sumvir[i]/(6.0*volcut))+(sumek[i]/(1.5*volcut)));
			if(a->flag&FLAG_GHOST) continue;
			if(a->flag&FLAG_ACTIVE) {
				avg_pres+=a->aux[pidx];
				navg++;
			}
		}
		
		avg_pres=TakeSUM(avg_pres);
		navg=TakeSUM(navg);
		avg_pres/=navg;

	}

	virtual void ExtractQuantities() {
		int na=GetNAtom();
		
		if(nalloc<na) {
			nalloc=na;
			MemRealloc(sumas,sizeof(double)*nalloc);
			MemRealloc(sumek,sizeof(double)*nalloc);
			MemRealloc(sumvx,sizeof(double)*nalloc);
			MemRealloc(sumvy,sizeof(double)*nalloc);
			MemRealloc(sumvz,sizeof(double)*nalloc);
			MemRealloc(sumvir,sizeof(double)*nalloc);
			MemRealloc(nneig,sizeof(int)*nalloc);
		}

		int nnn=sizeof(double)*na;
		memset(sumas,0,nnn);
		memset(sumas,0,nnn);
		memset(sumvx,0,nnn);
		memset(sumvy,0,nnn);
		memset(sumvz,0,nnn);
		memset(sumek,0,nnn);
		memset(sumvir,0,nnn);
		memset(nneig,0,sizeof(int)*na);

		ExtractTemperature();
		ExtractPressure();
	}

	virtual void MeasureSystemTempPres() {
		
		// consider central of mass motion
		// Calculate system temperature and pressure

		// calculate center of mass velocity
		double cmvx=0.0;
		double cmvy=0.0;
		double cmvz=0.0;
		
		double mix=DBL_MAX;
		double miy=DBL_MAX;
		double miz=DBL_MAX;

		double max=-DBL_MAX;
		double may=-DBL_MAX;
		double maz=-DBL_MAX;
		
		double summass=0.0;
		
		int na=GetNAtom();

		for(int i=0;i<na;i++){
			Atom *a=AtomPtr(i);
			if(a->flag&FLAG_GHOST) continue;
			if(a->flag&FLAG_ACTIVE) {
				double ma=GetMass(a);
				cmvx+=ma*a->vx;
				cmvy+=ma*a->vy;
				cmvz+=ma*a->vz;
				summass+=ma;
				if(mix>a->x)mix=a->x;
				if(miy>a->y)miy=a->y;
				if(miz>a->z)miz=a->z;
				if(max<a->x)max=a->x;
				if(may<a->y)may=a->y;
				if(maz<a->z)maz=a->z;
			}
		}

		cmvx=TakeSUM(cmvx);
		cmvy=TakeSUM(cmvy);
		cmvz=TakeSUM(cmvz);
		mix=TakeMIN(mix);
		max=TakeMAX(max);
		miy=TakeMIN(miy);
		may=TakeMAX(may);
		miz=TakeMIN(miz);
		maz=TakeMAX(maz);
		summass=TakeSUM(summass);

		cmvx/=summass;
		cmvy/=summass;
		cmvz/=summass;
		
		int totatom=0;
		double temp=0.0;

		for(int i=0;i<na;i++){
			Atom* a=AtomPtr(i);
			if(a->flag&FLAG_GHOST) continue;
			if(a->flag&FLAG_ACTIVE) {
				double ax=a->vx-cmvx;
				double ay=a->vy-cmvy;
				double az=a->vz-cmvz;
				temp+=GetMass(a)*(ax*ax+ay*ay+az*az);
				totatom++;
			}
		}
		// temp=2*sum{E_kinetic}

		temp=TakeSUM(temp);
		totatom=TakeSUM(totatom);

		system_temperature=Unit->Temperature(0.5*temp/(double)totatom);

		// using total volume (L^d), d=dimension
		
		double lx=max-mix;
		double ly=may-miy;
		double lz=maz-miz;
		double totvol=(lx>0.0?lx:1.0)*(ly>0.0?ly:1.0)*(lz>0.0?lz:1.0);
		system_pressure=Unit->Pressure((System->Virial+temp)/3.0/totvol);
		
	}

	bool PreIterationNode(int idx) {
		if(process==tempe) return CheckLoopTemp(idx);
		else return CheckLoopPressure(idx);
	}
	
	void IterationNode(int at, int to) {
		if(process==tempe) IterationNodeTemp(at,to);
		else IterationNodePressure(at,to);		
	}
	
	void Detect(){
		MeasureSystemTempPres();
		if(intensive_mode) ExtractQuantities();
		else if(OnTime(NextSample)) ExtractQuantities();
		DataDumper::Detect();
	}

};

#endif
