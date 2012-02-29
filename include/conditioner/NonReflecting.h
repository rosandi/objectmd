#ifndef _NON_REFLECTING_HPP_
#define _NON_REFLECTING_HPP_

#include <omd/conditioner.h>
#include <omd/comhandler.h>
using namespace omd;


class NonReflecting: public Conditioner, private ParallelGadget {
	
	double fconst;
	double imp;
	double termzone;
	double gridthick;
	
	int    n_layer;
	int    bot_force;
	int    top_force;
	bool   printforce;
	bool   evaluating;
	
	double *force_layer_average;
	int    *layer_population;
  
public:
  
	NonReflecting() {
		set_name("NON REFLECTING BOUNDARY (z bottom)");
		register_class("non_reflecting_boundary");		
		SetConditionerType(COND_PRE_CALCULATION|COND_FORCE_MODIFIER);
		evaluating=false;
	}
	
	~NonReflecting(){
		MemFree(force_layer_average);
		MemFree(layer_population);
	}
	
	void ReadParameter() {
    TargetName=SysParam->string_value("boundary.nrb.target");
		imp=SysParam->double_value("boundary.nrb.impedance");
		fconst=SysParam->double_value("boundary.nrb.force");
		gridthick=SysParam->double_value("boundary.nrb.layer");
	}
  	
	void Init(MDSystem* WorkSys){
		
		mdassert(WorkSys->type_of("SIMULATION SYSTEM GRID"),
           "this class requires SIMULATION SYSTEM GRID");
    
		Conditioner::Init(WorkSys);
		ParallelGadget::Init(WorkSys);
    
		top_force=ClaimAuxVariable(printforce, "ftop");
		bot_force=ClaimAuxVariable(printforce, "fbot");
    ter_layer=ClaimAuxVariable(false,"termlayer");
		
		double zmax=-DBL_MAX, zmin=DBL_MAX;

		for(int i=0;i<GetNAtom(); i++) {   
      if(zmax<Atoms(i).z)zmax=Atoms(i).z;
      if(zmin>Atoms(i).z)zmin=Atoms(i).z;
		}
		
		zmax=TakeMAX(zmax);
		zmin=TakeMIN(zmin);
		
		termzone=zmax+0.5*gridthick;
		n_layer=(int)((termzone-zmin)/gridthick);
		
		MemAlloc(force_layer_average, sizeof(double)*n_layer);
		MemAlloc(layer_population, sizeof(int)*n_layer);
		
		for(int i=0;i<GetNAtom(); i++) {
        int ly=(int)floor((Atoms(i).z-System->Box.z0)/gridthick);
        if(ly<0)ly=0;
        if(ly>=n_layer)ly=n_layer-1;
        Atoms(i).aux[ter_layer]=(double)ly;
    }
    
	}
	
	void PreCalculation(){
		int na=GetNAtom();   
    for(int i=0; i<na; i++){
      Atoms(i).aux[top_force]=0.0;
      Atoms(i).aux[bot_force]=0.0;
    }
	}
  
	void ForceModifier(){
    if(!evaluating) die("not evaluating force!");
    
    int na=GetNAtom();
    
    for(int i=0;i<n_layer;i++) {
      force_layer_average[i]=0.0;
      layer_population[i]=0;
    }
    
    for(int i=0;i<na;i++) {
      
      Atom* a=AtomPtr(i); 
      int nl=(int)a->aux[ter_layer];
      force_layer_average[nl]+=(a->aux[top_force]+imp*a->vz+fconst);
      layer_population[nl]++;
    
    }
    
    // reduce average...
    double* dtmp=new double[n_layer];
    int* itmp=new int[n_layer];
    
    MPI_Allreduce(force_layer_average,dtmp,n_layer,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(layer_population,itmp,n_layer,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    memcpy(force_layer_average,dtmp,n_layer*sizeof(double));
    memcpy(layer_population,itmp,n_layer*sizeof(int));
    
    delete[] dtmp;
    delete[] itmp;

    for(int i=0;i<n_layer;i++)
      force_layer_average[i]/=layer_population[i];
    
    for(int i=0;i<na;i++) {
      Atom* a=AtomPtr(i);
      a->fz=force_layer_average[(int)a->aux[ter_layer]];
    }
    
  }
	
	void EvaluateForce
  (Atom &a, Atom &b, double dx, double dy, double dz, double fr, double pot)
	{
		evaluating=true;
// TODO:		if(!(Target->Member(a) || Target->Member(b))) return;
    
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
