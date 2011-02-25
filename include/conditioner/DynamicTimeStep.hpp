
#ifndef _DYNAMIC_TIME_STEP_HPP_
#define _DYNAMIC_TIME_STEP_HPP_

#include <omd/conditioner.hpp>
#include <omd/integrator.hpp>
#include <omd/paragadget.hpp>

/**
  @ingroup conditioner
  @brief Dynamic Time Step

  Parameters:
  - timestep.max (t) : maximum allowed time step (normally in ps)
  - timestep.maxpath (l) : maximum allowed distance to move in one step (normally in Angstrom)
  - timestep.update (n) : the updating period (in step)

  If the maxpath parameter is not defined in the constructor, @e timestep.maxpath@e must exist.

**/

class DynamicTimeStep: public Pre_Conditioner, public ParallelGadget {
	OMD_FLOAT maxp;
	OMD_FLOAT maxdt;
	OMD_INT upd;
	
public:

	// maxpath is typically 10% lattice constant
	DynamicTimeStep(OMD_FLOAT maxpath=-1.0, OMD_FLOAT max_time_step=0.001, OMD_INT update_period=5) {
		set_name("DYNAMIC TIME STEP");
		register_class(get_name());
		maxp=maxpath;
		maxdt=max_time_step;
		upd=update_period;
	}
	
	void Init(MDSystem *WorkSys) {
		Pre_Conditioner::Init(WorkSys);
		ParallelGadget::Init(WorkSys);

		SysParam->peek("timestep.maxpath", maxp);
		SysParam->peek("timestep.max", maxdt);
		SysParam->peek("timestep.update", upd);
		assert(maxp>0, "maxpath is undefined! (parameter timestep.maxpath)");
	}

	void PreIntegration(){
		if(System->Step==1) return;
		if (!(System->Step%upd)) {
			OMD_FLOAT ndt;
			if (System->SqrMaxVelocity>0.0) {
    			ndt = maxp/(sqrt(System->SqrMaxVelocity)*1.10);
    			if(ndt>maxdt) ndt=maxdt;
			} else ndt=maxdt;
			ndt=TakeMIN(maxdt);
			System->Integrator->TimeStep=ndt;
			
		}
	}
};

#endif
