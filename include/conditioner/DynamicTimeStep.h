
#ifndef _DYNAMIC_TIME_STEP_HPP_
#define _DYNAMIC_TIME_STEP_HPP_

#include <omd/conditioner.h>
#include <omd/integrator.h>
#include <omd/paragadget.h>

using namespace omd;

/**
  @ingroup conditioner
  @brief Dynamic Time Step

  Parameters:
  - timestep.max (t) : maximum allowed time step (normally in ps)
  - timestep.maxpath (l) : maximum allowed distance to move in one step (normally in Angstrom)
  - timestep.update (n) : the updating period (in step)

  If the maxpath parameter is not defined in the constructor, @e timestep.maxpath @e must exist.

**/

class DynamicTimeStep: public PreConditioner, public ParallelGadget {
	double maxp;
	double maxdt;
	int upd;
	
public:

	// maxpath is typically 10% lattice constant
	DynamicTimeStep(double maxpath=-1.0, double max_time_step=0.001, int update_period=5) {
		set_name("DYNAMIC TIME STEP");
		register_class(get_name());
		maxp=maxpath;
		maxdt=max_time_step;
		upd=update_period;
	}
	
	void ReadParameter() {
		SysParam->peek("timestep.maxpath", maxp);
		SysParam->peek("timestep.max", maxdt);
		SysParam->peek("timestep.update", upd);
	}
	
	void Init(MDSystem *WorkSys) {
		PreConditioner::Init(WorkSys);
		ParallelGadget::Init(WorkSys);
		
		// calculate timestep first time...
		double ndt;
		if (System->SqrMaxVelocity>0.0) {
			ndt = maxp/(sqrt(System->SqrMaxVelocity)*1.10);
			if(ndt>maxdt) ndt=maxdt;
		} else ndt=maxdt;
		ndt=TakeMIN(ndt);
		System->Integrator->TimeStep=ndt;
		mdassert(maxp>0.0, "maxpath is undefined! (parameter timestep.maxpath)");
	}

	void PreIntegration(){
		if (!(System->Step%upd)) {
			double ndt;
			if (System->SqrMaxVelocity>0.0) {
    			ndt = maxp/(sqrt(System->SqrMaxVelocity)*1.10);
    			if(ndt>maxdt) ndt=maxdt;
			} else ndt=maxdt;
			ndt=TakeMIN(ndt);
			System->Integrator->TimeStep=ndt;			
		}
	}
};

#endif
