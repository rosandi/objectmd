/* omdlib
 ********************************************************
 *
 * ObjectMD
 * Molecular Dynamics Class Library
 * Version 2.0 (2009)
 *
 * Released under:
 *       GNU GENERAL PUBLIC LICENSE Version 3
 *       Copyright (C) 2005  Yudi Rosandi
 *
 * Read license.txt file in the source's root directory.
 *
 ********************************************************
 * Object-MD Implementations
 *
 * Implementation of Integrator class
 *
*/

#include <cmath>
#include <omd/integrator.hpp>
#include <omd/iterator.hpp>
#include <omd/system.hpp>

//----------Integrator----------//

/**
 * The integrator accepts a double parameter, which is the integration time step.
 */

MDIntegrator::MDIntegrator(OMD_FLOAT time_step)
{
	Iterator = NULL;
	ForcID=0; // Force kernel counter
	TimeStep = time_step;
    MaxCutRadius  = 0.0;
	register_class("INTEGRATOR");
	set_name("VERLET_INTEGRATOR");
    for (OMD_INT i=0;i<MAXATOMTYPE;i++) 
    	for (OMD_INT j=0;j<MAXATOMTYPE;j++) Forces[i][j]=NULL;
}

MDIntegrator::~MDIntegrator() {
	for (OMD_SIZET i=0; i<ActForces.size(); i++) delete ActForces[i];
}

void MDIntegrator::PrintInfo(ostream& ost) {
    ost<<get_name()<< " - interaction forces:\n\n";
    
	for (OMD_SIZET i=0; i<(System->SystemAtoms.size()); i++)
		for (OMD_SIZET j=i; j<(System->SystemAtoms.size()); j++) {
			if (GetForce(i, j)!=NULL) {
				ost << System->SystemAtoms[i]->get_name() << "*" 
					<< System->SystemAtoms[j]->get_name() << " = ";
                GetForce(i,j)->PrintInfo(ost);
                ost << "\n";
				}
			}
}

/**
 * The initialization function. In the integrator, this function has a task
 * to initiate all force kernels. The forces are also refered in an array
 * (STL vector) ActForces to ease access. If there is only one force kernel
 * inserted to the system, the id of all atoms in the system will be set to
 * zero, corresponds to the force kernel of element (0,0).
 * 
 * If the time step is not given to th constructor (negative value), it will
 * be initialize with the default value 1e-3.
 * 
*/

void MDIntegrator::Init(MDSystem* WorkSys) {	
    MDGadget::Init(WorkSys);

	assert(ActForces.size()>0, "no defined force kernel");
    NType = WorkSys->SystemAtoms.size();
    
	// put forces into the kernel matrix
	for(OMD_SIZET i=0;i<ActForces.size();i++){
		OMD_SIZET from=ActForces[i]->A;
		OMD_SIZET to  =ActForces[i]->B;

		assert(from<NType&&to<NType, "wrong atom IDs in the interaction force ("+
		       as_string(from)+","+as_string(to)+")");
		
		if(Forces[from][to]) die("interaction force for "+
			WorkSys->SystemAtoms[from]->get_name()+"["+as_string(from)+"] and "+
			WorkSys->SystemAtoms[to]->get_name()+"["+as_string(to)+"] is already assigned by "+
			Forces[from][to]->get_name());

		Forces[from][to]=Forces[to][from]=ActForces[i];
	}

	// using one interaction force, all atoms has the same ID=0
	if (ActForces.size()==1) { 
		WorkSys->set_id(0); 
		Forces[0][0]->SetAtomID(0, 0);
	}
	
	MaxCutRadius=-1.0;
	for(OMD_SIZET i=0;i<ActForces.size();i++) {
		ActForces[i]->Init(WorkSys);
		if (ActForces[i]->CutRadius>MaxCutRadius)
			MaxCutRadius=ActForces[i]->CutRadius;
	}

	assert(MaxCutRadius>0.0, "bad cut radius: "+as_string(MaxCutRadius));

	if(TimeStep<=0.0) {
		warn("using default initial time step (dt=unit_time/1000)");
		TimeStep=0.001;
	}

}

/**
 * Check the status of the force matrix. If the matrix is not completely filled
 * by the required forces for all type of atoms, it will terminate the program.
 * The force for the same type may be empty if only one atom of that type
 * is available in the system. This is usefull, for instance, for an ion
 * projectile.
*/
	
bool MDIntegrator::Check() {
	for(OMD_SIZET i=0;i<NType;i++)
	  for(OMD_SIZET j=0;j<NType;j++)
	    if(Forces[i][j]==NULL) {
		  if(i==j && System->SystemAtoms[i]->GetNAtom()<=1)continue;		  
		  die("incomplete interaction between "+
		      System->SystemAtoms[i]->get_name()+" and "+
		      System->SystemAtoms[j]->get_name());
        }
	return true;
}

/**
 * Add a force kernel to the integrator using the ids of the atoms to detemine
 * the types of the interacting atoms. The force kernel (ForceKernel) is stored 
 * in the force matrix. The id of the atoms indicate the indices of the matrix.
 * 
 * <center> <table border=1>
 * 		 <tr><td> 00 </td> <td> 01 </td></tr>
 * 		 <tr><td> 10 </td> <td> 11 </td></tr>
 * </table></center>
 * 	
 * The interaction is commutative, so the upper and lower diagonal keep the 
 * reference to the same force.
 * 
 * that atoms of type i and j interact the same
 * way as between j and i. On the diagonal are forces for interaction between
 * the same type of atoms.
 * 
*/

MDIntegrator* MDIntegrator::AddForce(ForceKernel *forc) 
{
	OMD_INT from=forc->A; OMD_INT to=forc->B;
	assert(from<MAXATOMTYPE&&to<MAXATOMTYPE, "MAXATOMTYPE exceeded");
	forc->set_id(ForcID++);
	forc->SetUnit(Unit);
	ActForces.push_back(forc);
	return this;
}

/**
 * This function is the iteration node invoked by the iterator.
 * The force iteration calls the CheckCompute() function of the force
 * kernel. The indices of the interacting atoms are determined by the iterator
 * class.
 */

void MDIntegrator::IterationNode(OMD_SIZET at, OMD_SIZET to)
{	
	int atid=Atoms(at).id;
	int toid=Atoms(to).id;
	Forces[(OMD_SIZET)(atid)][(OMD_SIZET)(toid)]->CheckCompute(at,to,atid,toid);
}

/**
 * @brief Integrator iteration
 * 
 * This function trigger the force calculation iteration by the iterator
 * (MDIterator). Before the loop all accumulators are cleared.
 * After the loop the Correction() functions of all active forces
 * are executed.
 * 
 * The iteration is bracketed by the two conditioners:
 * pre-calculation and force-modifier. This must be taken care when the function
 * is reimplemented.
 * 
*/

void MDIntegrator::Iterate() {
	SyncData(SYNC_SPACE);
    System->ExecuteConditioners(COND_PRE_CALCULATION);
	OMD_INT na=GetNAtom();
	for(OMD_INT i=0; i<na; i++){
		Atom* a=AtomPtr(i);
		a->fx=a->fy=a->fz=a->virial=a->potential=0.0;		
	}
	for(OMD_SIZET i=0;i<ActForces.size();i++) ActForces[i]->ClearAccumulators();
	Iterator->Iterate(this);
	for (OMD_SIZET i=0; i<ActForces.size(); i++) ActForces[i]->Correction();
	SyncData(SYNC_FORCE);
    System->ExecuteConditioners(COND_FORCE_MODIFIER);
}

/**
 * This is the Verlet integration algorithm. The iterration through all atoms
 * is done by the iterator (MDIterator class). The task of this function is
 * to update the velocities and positions of atoms.
 * 
*/

void MDIntegrator::Integrate() {

	OMD_FLOAT half_dt=0.5*TimeStep;

	// first velocity and position update
	OMD_INT natom=System->GetLocalAtomNumber();

    for (OMD_INT i=0;i<natom;i++) {
    	if(CheckActive(i)) {
			Atom* a=AtomPtr(i);
			OMD_FLOAT Mass = GetMass(a);
			a->vx+=half_dt*a->fx/Mass;
        	a->vy+=half_dt*a->fy/Mass;
        	a->vz+=half_dt*a->fz/Mass;
			a->x+=TimeStep*(a->vx);
        	a->y+=TimeStep*(a->vy);
        	a->z+=TimeStep*(a->vz);
		}
    }
	
	// force iteration
    Iterate();

	// second velocity update
    for(OMD_INT i=0;i<natom;i++) {
    	if(CheckActive(i)){
			Atom* a=AtomPtr(i);
			OMD_FLOAT Mass = GetMass(a);
			a->vx+=half_dt*a->fx/Mass;
        	a->vy+=half_dt*a->fy/Mass;
        	a->vz+=half_dt*a->fz/Mass;
		}
    }

}
