#ifndef _NON_REFLECTING_HPP_
#define _NON_REFLECTING_HPP_

#include <omd/conditioner.hpp>
#include <omd/comhandler.hpp>

class NonReflecting: public Conditioner {
	CommunicationHandler* Communicator;
	
	OMD_FLOAT fconst;
	OMD_FLOAT imp;
	OMD_FLOAT termzone;
	OMD_FLOAT gridthick;
	
	OMD_INT    n_layer;
	OMD_INT    tzid;
	OMD_INT    tzflag;
	OMD_INT    bot_force;
	OMD_INT    top_force;
	bool   printforce;
	bool   evaluating;
	
	OMD_FLOAT *force_layer_average;
	OMD_INT    *layer_population;

public:

	NonReflecting(
	    OMD_FLOAT impedance,
        OMD_FLOAT force_constant,
        OMD_INT    terminating_zone_id,
        OMD_FLOAT grid_layer_thickness,
        bool print_force=false)
	{
		imp=impedance;
		fconst=force_constant;
		tzid=terminating_zone_id;
		gridthick=grid_layer_thickness;
		printforce=print_force;
		
		set_name("NON REFLECTING BOUNDARY (z bottom)");
		register_class("non_reflecting_boundary");
		
		SetConditionerType(COND_PRE_CALCULATION|COND_FORCE_MODIFIER);
		evaluating=false;
	}
	
	~NonReflecting(){
		MemFree(force_layer_average);
		MemFree(layer_population);
	}
	
	void Init(MDSystem* WorkSys){
		
		assert(WorkSys->type_of("SIMULATION_SYSTEM_GRID"),
			"this class requires SIMULATION_SYSTEM_GRID");
			
		Conditioner::Init(WorkSys);
		Communicator=dynamic_cast<MDSystemGrid*>(WorkSys)->GetCommunicator();
		
		tzflag=ClaimFlagBit();
		top_force=ClaimAuxVariable(printforce, "ftop");
		bot_force=ClaimAuxVariable(printforce, "fbot");
		
		OMD_FLOAT zmax=-DBL_MAX, zmin=DBL_MAX;
		for(OMD_INT i=0;i<GetNAtom(); i++) {   
			if(Atoms(i).xid==tzid) {
				if(zmax<Atoms(i).z)zmax=Atoms(i).z;
				if(zmin>Atoms(i).z)zmin=Atoms(i).z;
			}
		}
		
		zmax=Communicator->TakeMAX(zmax);
		zmin=Communicator->TakeMIN(zmin);
		
		termzone=zmax+0.5*gridthick;
		n_layer=(OMD_INT)((termzone-zmin)/gridthick);
		
		MemAlloc(force_layer_average, sizeof(OMD_FLOAT)*n_layer);
		MemAlloc(layer_population, sizeof(OMD_INT)*n_layer);
		
		for(OMD_INT i=0;i<GetNAtom(); i++) {
            if(Atoms(i).z<=termzone) {
            	
                OMD_INT ly=(OMD_INT)floor((Atoms(i).z-System->Box.z0)/gridthick);
                if(ly<0)ly=0;
                if(ly>=n_layer)ly=n_layer-1;
                Atoms(i).xid=ly;
                SetFlag(i, tzflag);
                
            } else UnsetFlag(i, tzflag);
        }
       		
	}
	
	void PreCalculation(){
		OMD_INT na=GetNAtom();   
        for(OMD_INT i=0; i<na; i++){
            Atoms(i).aux[top_force]=0.0;
            Atoms(i).aux[bot_force]=0.0;
        }

	}

	void ForceModifier(){
		if(!evaluating) die("not evaluating force!");
		
        OMD_INT na=GetNAtom();

        vector<OMD_INT> blist[n_layer];

        for(OMD_INT i=0;i<n_layer;i++) {
        	force_layer_average[i]=0.0;
        	layer_population[i]=0;
		}
        
        for(OMD_INT i=0;i<na;i++) {
            Atom* a=AtomPtr(i);
            if(a->flag&tzflag){
            	
            	if(a->xid<0||a->xid>=n_layer)
            		die("wrong member of terminating zone "+
            		    as_string(a->xid)+" index="+as_string(i));

                blist[a->xid].push_back(i);
                force_layer_average[a->xid]+=(a->aux[top_force]+imp*a->vz+fconst);
                layer_population[a->xid]++;

            }
        }

        for(OMD_INT l=0;l<n_layer;l++) {
        	OMD_FLOAT favg=force_layer_average[l]/(OMD_FLOAT)layer_population[l];
            for(OMD_INT i=0;i<blist[l].size();i++) {
                Atoms(blist[l].at(i)).fz=favg;
            }
        }

    }
	
	void EvaluateForce
		(Atom &a, Atom &b, OMD_FLOAT dx, OMD_FLOAT dy, OMD_FLOAT dz, OMD_FLOAT fr, OMD_FLOAT pot)
	{
		evaluating=true;
		
		if(a.z<=b.z) {
            a.aux[top_force]+=fr*dz;
            b.aux[bot_force]-=fr*dz;
        } else {
            a.aux[bot_force]+=fr*dz;
            b.aux[top_force]-=fr*dz;
        }

	}

};


#endif
