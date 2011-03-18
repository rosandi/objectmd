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

#include <omd/paragadget.hpp>
#include <omd/conditioner.hpp>

#define INPUT_PULSE    1
#define INPUT_RAMP     2
#define INPUT_FUNCTION 3

/**
 @brief Energy source
 
 This object input energy to the system as  pulse, ramp, or using
 any shape defined by a normallized function table.
 
 Parameters:
 - @e source.energy (energy in eV): total energy given to the system
 - @e source.mode (pulse|ramp|function)  : input mode. 
 - @e source.start (time in ps) : time to start the energy input. If this
   parameter is zero, the atoms' velocity will be initiated by a corresponding
   amount of energy (\f \Delta E \f) at the first step. The target crystal
   is assumed to have initially zero temperature.
 - @e source.width (width in ps) : the width of the pulse. Not applicable for
   @e mode=function.
 - @e source.function (filename) : the omd-table filename of a normallized
   function.
 
**/

class EnergySource: public Pre_Conditioner, public ParallelGadget {
	int input_mode;
	OMD_FLOAT e_tally; // total energy/atom...
	OMD_FLOAT delta;
	OMD_FLOAT energy;
	OMD_FLOAT start_time;
	OMD_FLOAT stop_time;
	OMD_FLOAT width;
	OMD_FLOAT grad;
	string smode;
	string sfunc; // function must be normallized...
	bool firstcall;
	
public:
	EnergySource() {
		set_name("ENERGY SOURCE");
		register_class(get_name());
		firstcall=true;
		e_tally=0.0;
	}
	
	void ReadParameter() {
		energy=SysParam->double_value("source.energy"); // energy per atom...
		SysParam->peek("source.start", start_time, 0.0);
		SysParam->peek("source.mode", smode, "pulse");
		SysParam->peek("source.width", width, -1.0);
		SysParam->peek("source.function", sfunc, "-");
	}
	
	bool CheckParameter() {
		
		double dt=GetTimeStep();
		
		if(smode=="pulse") {
			if(width<dt) width=dt;
			input_mode=INPUT_PULSE;
			delta=energy*dt/width;
			stop_time=start_time+width;
		}
		
		if(smode=="ramp") {
			die("mode=ramp not implemeted yet");
		}
		
		if(smode=="function") {
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
			OMD_FLOAT small_vel=System->GetMaxVelocity()*1e-6;
			small_vel*=small_vel;
			for(int i=0;i<na;i++){
				Atom* a=AtomPtr(i);
				OMD_FLOAT m=GetMass(a);						
				if(a->flag&FLAG_GHOST) continue;
				if(!(a->flag&FLAG_ACTIVE)) continue;
				OMD_FLOAT c2=a->vx*a->vx + a->vy*a->vy + a->vz*a->vz;
				if (c2<small_vel) continue;
				OMD_FLOAT chi=sqrt(1.0+2.0*delta/c2/m);
				a->vx*=chi;
				a->vy*=chi;
				a->vz*=chi;
			}
		}
		e_tally+=delta;

	}

};

#endif
