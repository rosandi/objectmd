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
 *
*/

#ifndef _ENERGY_INPUT_HPP_
#define _ENERGY_INPUT_HPP_

#include <omd/paragadget.h>
#include <omd/conditioner.h>

using namespace omd;


#define INPUT_PULSE    1
#define INPUT_FUNCTION 2

/**
 @brief Energy source
 
 This class inputs energy to the system as a pulse or any shape 
 defined by a normallized function table. The energy is pushed
 into the system using @e velocity @e scaling method, 
 \f$ v' \leftarrow \gamma v \f$. The scaling factor
 \f$\gamma \f$ is calculated by the following expression:
 
 \f[
 \gamma = \sqrt{1 + \frac{2\Delta E}{m v^2}}
 \f]
 
 Parameters:
 - @e source.energy (energy in eV): total energy given to the system
 - @e source.mode (pulse|function)  : input mode. 
 - @e source.start (time in ps) : time to start the energy input. If this
   parameter is zero, the atoms' velocity will be initiated by a corresponding
   amount of energy (\f$ \Delta E \f$) in the first step. The target crystal
   is assumed to have initially zero temperature.
 - @e source.width (width in ps) : the width of the pulse. Not applicable for
   @e mode=function.
 - @e source.function (filename) : the omd-table filename of a normallized
   function.
 
**/

class EnergySource: public PreConditioner, public ParallelGadget {
	int input_mode;
	double e_tally; // total energy/atom...
	double delta;
	double energy;
	double start_time;
	double stop_time;
	double width;
	double grad;
	string smode;
	string sfunc; // function must be normallized...
	bool firstcall;
	
public:
	EnergySource() {
		set_name("source");
		register_class(get_name());
		firstcall=true;
		e_tally=0.0;
	}
	
	void ReadParameter() {
		energy=SysParam->double_value(mytag("energy")); // energy per atom...
		SysParam->peek(mytag("start"), start_time, 0.0);
		SysParam->peek(mytag("mode"), smode, "pulse");
		SysParam->peek(mytag("width"), width, -1.0);
		SysParam->peek(mytag("function"), sfunc, "-");
	}
	
	bool CheckParameter() {
		
		double dt=GetTimeStep();
		
		if(smode=="pulse") {
			input_mode=INPUT_PULSE;
			if(width<dt) width=dt;
			delta=energy*dt/width;
			stop_time=start_time+width;
		}
		
		if(smode=="function") {
			input_mode=INPUT_FUNCTION;
			die("mode=function not implemented yet");
		}
		
		return true;
	}
	
	void Init(MDSystem* WorkSys) {
		Conditioner::Init(WorkSys);
		ParallelGadget::Init(WorkSys);

		if(SearchGadgetType("dynamic time step", false)&&smode!="sudden")
			die("incompatible with (DYNAMIC TIME STEP)");
	}
	
	void PreIntegration() {

		if(System->ElapsedTime<start_time) return;
		
		if(firstcall) {
			firstcall=false;
			if(start_time==0.0) System->SetKineticEnergy(delta);
			return;
		}
		
		if(System->ElapsedTime<stop_time) {
			int na=GetNAtom();
			double small_vel=System->GetMaxVelocity()*1e-6;
			small_vel*=small_vel;
			for(int i=0;i<na;i++){
				Atom* a=AtomPtr(i);
				double m=GetMass(a);						
				if(a->flag&FLAG_GHOST) continue;
				if(!(a->flag&FLAG_ACTIVE)) continue;
				double c2=a->vx*a->vx + a->vy*a->vy + a->vz*a->vz;
				if (c2<small_vel) continue;
				double chi=sqrt(1.0+2.0*delta/c2/m);
				a->vx*=chi;
				a->vy*=chi;
				a->vz*=chi;
			}
		}
		e_tally+=delta;

	}
	
	void PrintInfo(ostream& ost) {
		ost << "id."<<id<<" "<<get_name()<<": "<< "mode="<<smode;
		if(smode=="pulse") {
			ost <<" energy="<<energy
				<<" width="<<width
				<<" start="<<start_time<<"\n";
		}
		if(smode=="function") {
			ost <<" energy="<<energy<<" function="<<sfunc<<" start="<<start_time<<"\n";
		}
	}

};

#endif
