//-------Temperature-Pressure-Density--Detector-----------//

#ifndef _THERMO_DETECTOR_HPP_
#define _THERMO_DETECTOR_HPP_

#include <omd/dataslot.hpp>
#include <omd/iterator.hpp>
#include <omd/paragadget.hpp>
#include <detector/DataDumper.hpp>

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

	OMD_FLOAT *sumas,*sumek,*sumvir;
	OMD_FLOAT *sumvx,*sumvy,*sumvz;
	OMD_FLOAT ek,masat,masto,ercut,volcut;
	OMD_FLOAT system_temperature,system_pressure; // System's temperature and pressure
	OMD_FLOAT avg_temp; // average of atom temperature
	OMD_FLOAT avg_pres; // average of atom pressure
	OMD_INT    *nneig;
	OMD_INT    nalloc;
	OMD_INT    tidx,pidx,nidx;
	bool   intensive_mode; // measure every timestep
	enum {tempe,press} process;

public:
	ThermoDetector(OMD_FLOAT tm=-1.0, OMD_FLOAT rcut=-1.0, string fn="Data", OMD_INT extlen=3)
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
	
	void SetTemperature(OMD_INT idx, OMD_FLOAT t) {Atoms(idx).aux[tidx]=t;}
	OMD_FLOAT GetTemperature(OMD_INT idx) {return Atoms(idx).aux[tidx];}
	OMD_FLOAT GetTemperature(){return system_temperature;}
	
	/**
	 * read the average temperature of the atoms. monomers and dimers are
	 * excluded.
	 */

	OMD_FLOAT GetTemperatureAvg(){
		assert(intensive_mode, "attempt to invoke GetTemperatureAvg in non-intensive mode!");
		return avg_temp;
	}
	
	OMD_FLOAT GetPressure(OMD_INT idx) {return Atoms(idx).aux[pidx];}
	OMD_FLOAT GetPressure(){return system_pressure;}
	OMD_FLOAT GetCutRadius(){return ercut;}
	void GetCMVelocities(OMD_FLOAT* &cvx, OMD_FLOAT* &cvy, OMD_FLOAT* &cvz){
		cvx=sumvx;cvy=sumvy;cvz=sumvz;
	}
	
	OMD_INT GetAuxTempIndex(){return tidx;}
	OMD_INT GetAuxPresIndex(){return pidx;}
	OMD_INT GetAuxDensIndex(){return nidx;}
	OMD_FLOAT GetDetectVolume(){return volcut;}

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
		MemAlloc(sumas,sizeof(OMD_FLOAT)*nalloc);
		MemAlloc(sumek,sizeof(OMD_FLOAT)*nalloc);
		MemAlloc(sumvx,sizeof(OMD_FLOAT)*nalloc);
		MemAlloc(sumvy,sizeof(OMD_FLOAT)*nalloc);
		MemAlloc(sumvz,sizeof(OMD_FLOAT)*nalloc);
		MemAlloc(sumvir,sizeof(OMD_FLOAT)*nalloc);
		MemAlloc(nneig,sizeof(OMD_INT)*nalloc);
		if(ercut<0.0) ercut=System->GetMaxCutRadius();
		assert(ercut>0.0, "cut radius is zero... static mode?");
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

	virtual bool CheckLoopTemp(OMD_SIZET idx) {
		Atom* a=AtomPtr(idx);
		nneig[idx]++;
		masat=GetMass(a->id);
		sumas[idx]+=masat;
		sumvx[idx]+= masat*a->vx;
		sumvy[idx]+= masat*a->vy;
		sumvz[idx]+= masat*a->vz;			
		ek=0.5*masat*(a->vx*a->vx+a->vy*a->vy+a->vz*a->vz);
		sumek[idx]+=ek;
		return true;
	}

	virtual void IterationNodeTemp(OMD_SIZET at, OMD_SIZET to) {
		OMD_FLOAT d=CalcDistance(at,to);
		Atom* a=AtomPtr(at);
		Atom* b=AtomPtr(to);
		if (d<=ercut) {
			masto=GetMass(b->id);	
			nneig[to]++;
			nneig[at]++;
			sumas[to]+=masat;
			sumas[at]+=masto;
			sumvx[to]+=masto*a->vx;
			sumvy[to]+=masto*a->vy;
			sumvz[to]+=masto*a->vz;                
			sumvx[at]+=masat*b->vx;
			sumvy[at]+=masat*b->vy;
			sumvz[at]+=masat*b->vz;                
			sumek[to]+=ek;
			sumek[at]+=0.5*masto*(b->vx*b->vx+b->vy*b->vy+b->vz*b->vz);
		}
	}

	virtual void ExtractTemperature() {
		process=tempe;
		Iterator->Iterate(this);
		OMD_INT na=GetNAtom();
		OMD_INT navg=0;
		avg_temp=0.0;
		
		for (OMD_INT i=0;i<na;i++) {
			Atom* a=AtomPtr(i);
			sumvx[i]/=sumas[i]; 
			sumvy[i]/=sumas[i];
			sumvz[i]/=sumas[i];
			sumek[i]-=0.5*sumas[i]*
				      (sumvx[i]*sumvx[i]+sumvy[i]*sumvy[i]+sumvz[i]*sumvz[i]);				      
			a->aux[tidx]=Unit->Temperature(sumek[i]/(OMD_FLOAT)nneig[i]);
			a->aux[nidx]=(OMD_FLOAT)nneig[i]/volcut;

			if(nneig[i]<2) continue;
			if(a->flag&FLAG_GHOST) continue;
			
			if(a->flag&FLAG_ACTIVE) {
				avg_temp+=Atoms(i).aux[tidx];
				navg++;
			}

		}

		avg_temp=TakeSUM(avg_temp);
		navg=TakeSUM(navg);
		if(navg>0) avg_temp/=(OMD_FLOAT)navg;
		else avg_temp=0.0;
		
	}

	virtual bool CheckLoopPressure(OMD_INT idx) {
		sumvir[idx]+=Atoms(idx).virial;  // my self
		return true;
	}

	virtual void IterationNodePressure(OMD_SIZET at, OMD_SIZET to) {
		OMD_FLOAT d=CalcDistance(at,to);
		if(d<=ercut){
			sumvir[at]+=Atoms(to).virial;
			sumvir[to]+=Atoms(at).virial;
		}
	}
	
	virtual void ExtractPressure() {
		process=press;
		Iterator->Iterate(this);
		OMD_INT na=GetNAtom();
		OMD_INT navg=0;
		avg_pres=0.0;
		for (OMD_INT i=0;i<na;i++) {
			Atom* a=AtomPtr(i);
			a->aux[pidx]=(nneig[i]==0)?0.0:
			Unit->Pressure((sumvir[i]/(6.0*volcut))+(sumek[i]/(1.5*volcut)));
			if(a->flag&FLAG_GHOST) continue;
			if(a->flag&FLAG_ACTIVE) {
				avg_pres+=a->aux[pidx];
				navg++;
				// 6=2*3 --> reduce double sum in virial done by force-kernel
			}
		}
		
		avg_pres=TakeSUM(avg_pres);
		navg=TakeSUM(navg);
		avg_pres/=navg;

	}

	virtual void ExtractQuantities() {
		OMD_INT na=GetNAtom();
		
		if(nalloc<na) {
			nalloc=na;
			MemRealloc(sumas,sizeof(OMD_FLOAT)*nalloc);
			MemRealloc(sumek,sizeof(OMD_FLOAT)*nalloc);
			MemRealloc(sumvx,sizeof(OMD_FLOAT)*nalloc);
			MemRealloc(sumvy,sizeof(OMD_FLOAT)*nalloc);
			MemRealloc(sumvz,sizeof(OMD_FLOAT)*nalloc);
			MemRealloc(sumvir,sizeof(OMD_FLOAT)*nalloc);
			MemRealloc(nneig,sizeof(OMD_INT)*nalloc);
		}

		OMD_INT nnn=sizeof(OMD_FLOAT)*na;
		memset(sumas,0,nnn);
		memset(sumas,0,nnn);
		memset(sumvx,0,nnn);
		memset(sumvy,0,nnn);
		memset(sumvz,0,nnn);
		memset(sumek,0,nnn);
		memset(sumvir,0,nnn);
		memset(nneig,0,sizeof(OMD_INT)*na);

		ExtractTemperature();
		ExtractPressure();
	}

	virtual void MeasureSystemTempPres() {
		
		// consider central of mass motion
		// Calculate system temperature and pressure

		// calculate center of mass velocity
		OMD_FLOAT cmvx=0.0;
		OMD_FLOAT cmvy=0.0;
		OMD_FLOAT cmvz=0.0;
		
		OMD_FLOAT mix=DBL_MAX;
		OMD_FLOAT miy=DBL_MAX;
		OMD_FLOAT miz=DBL_MAX;

		OMD_FLOAT max=-DBL_MAX;
		OMD_FLOAT may=-DBL_MAX;
		OMD_FLOAT maz=-DBL_MAX;
		
		OMD_FLOAT summass=0.0;
		
		OMD_INT na=GetNAtom();

		for(OMD_INT i=0;i<na;i++){
			Atom *a=AtomPtr(i);
			if(a->flag&FLAG_GHOST) continue;
			if(a->flag&FLAG_ACTIVE) {
				OMD_FLOAT ma=GetMass(a->id);
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
		
		OMD_INT totatom=0;
		OMD_FLOAT temp=0.0;

		for(OMD_INT i=0;i<na;i++){
			Atom* a=AtomPtr(i);
			if(a->flag&FLAG_GHOST) continue;
			if(a->flag&FLAG_ACTIVE) {
				OMD_FLOAT ax=a->vx-cmvx;
				OMD_FLOAT ay=a->vy-cmvy;
				OMD_FLOAT az=a->vz-cmvz;
				temp+=GetMass(a->id)*(ax*ax+ay*ay+az*az);
				totatom++;
			}
		}
		// temp=2*sum{E_kinetic}

		temp=TakeSUM(temp);
		totatom=TakeSUM(totatom);

		system_temperature=Unit->Temperature(0.5*temp/(OMD_FLOAT)totatom);

		// using total volume (L^d), d=dimension
		
		OMD_FLOAT lx=max-mix;
		OMD_FLOAT ly=may-miy;
		OMD_FLOAT lz=maz-miz;
		OMD_FLOAT totvol=(lx>0.0?lx:1.0)*(ly>0.0?ly:1.0)*(lz>0.0?lz:1.0);
		system_pressure=Unit->Pressure((System->Virial+temp)/3.0/totvol);
		
	}

	bool PreIterationNode(OMD_SIZET idx) {
		if(process==tempe) return CheckLoopTemp(idx);
		else return CheckLoopPressure(idx);
	}
	
	void IterationNode(OMD_SIZET at, OMD_SIZET to) {
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
