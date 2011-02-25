//? Name = DummyForce
//? Type = ForceKernel
//? Description = Force which does nothing
//? EmbedCommand = Forces->AddForce
//? HeaderFile = DummyForce.hpp

#ifndef _DUMMY_FORCE_HPP_
#define _DUMMY_FORCE_HPP_

#include <omd/forcekernel.hpp>

class DummyForce: public ForceKernel {
	public:	
	DummyForce(){CutRadius=1.0;}
	void Compute(Atom &at, Atom &to) {}
};

#endif
