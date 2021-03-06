// Homogen means: homogen electron temperature
// SC means: scaling-corrected
// N means: use density factor on coupling
// A means: local density of atom is used in scalling

#ifndef _TTM_HOMOGEN_HPP_
#define _TTM_HOMOGEN_HPP_

#include <iostream>
#include <omd/paragadget.h>
#include <omd/modify.h>
#include <omd/treader.h>
#include <omd/param.h>
#include <detect/ThermoDetector.h>
using namespace omd;


/**
 * @ingroup modify
 * \brief Homogen Two Temperature Model to transfer energy 
 *        from electronic system to the lattice
 * 
 * No electron temperature diffusion (homogen electron temperature).
 * The energy is given to the electronic system and transfered to the latice via
 * the coupling function. Used in laser ablation simulation. The energy is given to the
 * lattice using corrected velocity-scalling method. The scalling factor is calculated
 * using:
 *  @f[ C_m \cdot c + \sqrt{(C_m\cdot c)^2 + c^2(\epsilon + V^2 - C_m^2)}/c^2 @f]
 * 
 * where @f$ \epsilon=2\Delta E/m @f$ and @f$ \Delta E=G(T_e-T_a)\Omega\Delta t @f$ is the 
 * input energy.
*/

// explain: table file format, pulse function file, unity

class TTM_Homogen: public PreModify, public ParallelGadget {
	
	TableReader G;
	TableReader Ce;
	TableReader Ft;
	double E0;
	double G0;
	double C0;
	double pulse_width;
	double pulse_offset;
	double fd_time_step;
	double electron_energy;
	double electron_temperature;
	double atom_volume;
	double e_source; // timestep source 
	double de_energy; // timestep electron energy
	double da_energy; // timestep one atom energy
	double int_Een; // integrated (total) electron energy
	double int_Aen; // integrated (total) input energy to lattice
	double stop_time;
	
	bool use_density;
	bool linear_density;
	double n_zero;
	double low_density;
	double low_density_factor;
	
	ThermoDetector* temp_detector;
	bool balance_el_lattice;
	bool do_coupling;
	bool firstrun;
	bool firstfd;
	bool stopped;
	bool stop_equ;
	
	int    allocated_length;
	int    atemp;
	int    adens;
	
	enum {none, function, pulsetemp, pulse} source_type;
	
public:
	
	TTM_Homogen():PreModify(){
		set_name("ttmsc");
		register_class("TWO_TEMPERATURE_MODEL");
		electron_temperature=-1.0; // undefined...
		electron_energy=-1.0; // undefined...
		e_source=0.0;
		int_Een=0.0;
		int_Aen=0.0;
		allocated_length=0;
		do_coupling=true;
		firstrun=true;
		firstfd=true;
		stop_equ=false;
		source_type=none;
		stopped=false;
		use_density=false;
		linear_density=false;
		n_zero=0;
	}

	virtual ~TTM_Homogen(){}
	
	void ReadParameter(){
		SysParam->peek(mytag("energy_density"), E0, -1.0); // energy / Volume (A^3)
		SysParam->peek(mytag("temperature_input"), E0, -1.0);
		SysParam->peek(mytag("atom_volume"), atom_volume, -1.0);
		SysParam->peek(mytag("pulse_width"), pulse_width, -1.0);
		SysParam->peek(mytag("pulse_offset"), pulse_offset, 0.0);

		if(SysParam->exist(mytag("electron_init_energy")))
			electron_energy=SysParam->double_value(mytag("electron_init_energy"));
		else
			SysParam->peek(mytag("electron_init_temperature"), electron_temperature, -1.0);

		SysParam->peek(mytag("balance_electron_lattice"), balance_el_lattice, false);
		SysParam->peek(mytag("stop_at"), stop_time, -1.0);
		SysParam->peek(mytag("stop_equ"), stop_equ, false);
		if(SysParam->exist(mytag("interpolate_density"))) {
		  use_density=true;
		  linear_density=true;
		} else {
		    SysParam->peek(mytag("low_density"), low_density, -1.0);
		    if(low_density>0) {
		        use_density=true;
		        SysParam->peek(mytag("low_density_factor"), low_density_factor, 1.0);
		    } 
        }
	}

	bool CheckParameter(){
		blog("two-temperature-model: checking parameters...", LOGCREATE);
		if(SysParam->exist(mytag("source_function_file"))) source_type=function;
		else if(SysParam->exist(mytag("temperature_input"))) source_type=pulsetemp;
		else if(pulse_width>0.0) source_type=pulse;
		return true;
	}

	// requires temperature and pressure detector
	void Init(MDSystem* WorkSys) {
		PreModify::Init(WorkSys);
		ParallelGadget::Init(WorkSys);
		
		RegisterMessageSlot(new DataSlot("e_temp:",2))
			->SetFormat("%0.3E")->SetData(electron_temperature);

		RegisterMessageSlot(new DataSlot("e_source:",2))
				->SetFormat("%0.3E")->SetData(e_source);

		RegisterMessageSlot(new DataSlot("int_Een:",2))
			->SetFormat("%0.5E")->SetData(int_Een); // total energy given from electron

		RegisterMessageSlot(new DataSlot("int_Aen:",2))
		->SetFormat("%0.5E")->SetData(int_Aen); // total average energy given to lattice

		G.open(SysParam->string_value(mytag("electron-phonon_file")), "ELECTRON_PHONON_COUPLING");
		G0=G.param.double_value("Constant");
		Ce.open(SysParam->string_value(mytag("electron-phonon_file")), "ELECTRON_HEAT_CAPACITY");
		C0=Ce.param.double_value("Constant");
		
		if(source_type==function)
			Ft.open(SysParam->string_value(mytag("source_function_file")));

		fd_time_step=WorkSys->Integrator->TimeStep;
		
		temp_detector=dynamic_cast<ThermoDetector*>(SearchGadgetType("thermo",false));
		mdassert(temp_detector, "ThermoDetector is needed by "+get_name());
		temp_detector->SetIntensive();
		
		int_Een=0.0;
		mdassert(atom_volume>0.0, "define atom_volume!");
	
	}

	void PrintInfo(ostream& ost){
		// FIXME! check units!
		ost << "ID." << id << " " << get_name() << "\n";
		if(source_type==pulsetemp) ost << "temperature input = " << E0 << "K\n";
		else if(E0>0.0) ost << "energy input = " << E0 << "eV/A^3\n";

		if(electron_temperature>0.0)
			ost<< "electron initial temperature = " << electron_temperature << "K\n";
		else if(electron_energy>0.0)
			ost<< "electron initial energy = " << electron_energy << "eV\n";
        
        if(linear_density) {
            ost<<"density factor linearly interpolated\n";
        } else {
            if(use_density) {
                ost<<"low density limit = " << low_density << "\n"
                   <<"low density factor = " << low_density_factor <<"\n";
            } else {
                ost << "low density limit is not used!\n";
            }
        }
        
        if(stop_equ) {
            ost<<"TTM stopped if T_e<T_a\n";
        }

		ost << "pulse offset = " << pulse_offset << " ps\n";
	}

private:

	void find_temperature(double eng) {

		if(GetRank()==MDROOT) {
			blog("finding temperature...");
			double tn=eng/3.0/MD_BOLTZMANN;
			double en=0.5*C0*Ce.read(tn)*tn*atom_volume;
			double tm=tn;
			double em;
			double d0=eng-en;
			double d1=d0;
			double tstep=d0/fabs(d0)*5.0;
			double tt,et;
			
			// FIXME! do better algorithm!
			while(d1*d0>0) {
				tn=tm;
				d1=d0;
				tm+=tstep;
				en=0.5*C0*Ce.read(tm)*tm*atom_volume;
				d0=eng-en;
			}
			
		    electron_temperature=tm;
		    blog("electron temperature = "+as_string(electron_temperature));
		}
		SyncProcesses();
		Broadcast(electron_temperature);
	}

	void fd_main_loop(){
		int na=GetNAtom();
		double atomtemp=temp_detector->GetTemperatureAvg();
		double atomdens=temp_detector->GetDensityAvg();
		
		// use constant volume!
		
		if(firstfd){
			if(electron_temperature<0.0) { // zero electron temperature is allowed!
				if(electron_energy<0.0) electron_temperature=atomtemp;
				else find_temperature(electron_energy);
			}

			int_Een=0.5*electron_temperature*
			                        C0*Ce.read(electron_temperature)*atom_volume;

			int_Aen=(2.0*System->Kinetic)/System->GetTotalAtom();
			firstfd=false;
			n_zero=temp_detector->GetDensityAvg();
		}

		// Te
		double de=-G0*G.read(electron_temperature)*(electron_temperature-atomtemp)*fd_time_step;
		
		if(linear_density) {
		    de*=atomdens/n_zero;
		} else if(use_density) {
		    if(atomdens<(low_density*n_zero)) de*=low_density_factor;
		}
		
		int_Een+=de*atom_volume;
		electron_temperature+=de/(C0*Ce.read(electron_temperature));
		
		// stop ttm if the electron temperature is lower than atom temperature
		if (stop_equ) {
			if (electron_temperature<atomtemp) stopped=true;
		}

	}

	void velocity_scalling() {
		int na=GetNAtom();
		double *cmx, *cmy, *cmz;
		temp_detector->GetCMVelocities(cmx, cmy, cmz);
		double volcut=temp_detector->GetDetectVolume();
		double coup=G0*G.read(electron_temperature);
		
		// \Delta t_FD = \Delta t_MD
		double sum_eng=0.0;
		int neng=0;
		double small_vel=System->GetMaxVelocity()*1e-6;
		small_vel*=small_vel;
		
		for(int i=0;i<na;i++){
			Atom* a=AtomPtr(i);
			if (a->flag&FLAG_GHOST) continue;
			if (!(a->flag&FLAG_ACTIVE)) continue;
			if (a->aux[adens]*volcut<=3.0) continue;

			double m=GetMass(a);
			double de=coup*(electron_temperature - a->aux[atemp])*atom_volume*fd_time_step;
		    
		    // here things happen...
		    if(linear_density) {
		      de*=a->aux[adens]/n_zero;
		    } else if(use_density) {
		      if(a->aux[adens]<(low_density*n_zero)) de*=low_density_factor;
		    }

			double cx=(a->vx-cmx[i]);
			double cy=(a->vy-cmy[i]);
			double cz=(a->vz-cmz[i]);
			
			double c2=cx*cx + cy*cy + cz*cz;
			if (c2<small_vel) continue;
			double cmc=(cmx[i]*cx    +cmy[i]*cy    +cmz[i]*cz);
			double cm2=(cmx[i]*cmx[i]+cmy[i]*cmy[i]+cmz[i]*cmz[i]);
			double v2=(a->vx*a->vx + a->vy*a->vy + a->vz*a->vz);
			double rr=cmc*cmc+c2*(2.0*de/m+v2-cm2);
			if (rr<0.0) continue;
			
			double sa=(-cmc+sqrt(rr))/c2;
			
			a->vx=(cmx[i]+sa*cx);
			a->vy=(cmy[i]+sa*cy);
			a->vz=(cmz[i]+sa*cz);

			sum_eng+=de;
			neng++;

		}

		sum_eng=TakeSUM(sum_eng);
		neng=TakeSUM(neng);
		int_Aen+=(sum_eng/(double)neng);
		
	}

public:
	
	/**
	 * At first run force the temperature detector to do measurement.
	 * If pulse_offset is un-zero, electron temperature takes the lattice temperature.
	 * 
	 * Steps:
	 *    # check system temperature, using TempPressDetector
	 *    # run finite difference loop
	 *    # scale the velocities using center of mass method
	 * 
	 */
	
	void PreIntegration(){
		
		if (stopped) return;

		if(stop_time>0.0) {
			if(System->ElapsedTime>stop_time) {
				stopped=true;
				return;
			}
		}
		
		if(firstrun){  // early Detector execution
			atemp=temp_detector->GetAuxTempIndex();
			adens=temp_detector->GetAuxDensIndex();
			temp_detector->Detect();
		
			if(balance_el_lattice) {
				electron_temperature-=(2.0*temp_detector->GetTemperatureAvg());
				blog(string("balancing temperature with lattice")+
				     "\n electron="+as_string(electron_temperature)+
				     "\n lattice="+as_string(temp_detector->GetTemperatureAvg())+
				     " (avg)");
			}
			firstrun=false;
			return;
		}

		if(do_coupling&&System->ElapsedTime>pulse_offset) {
			fd_main_loop();
			velocity_scalling();
		}
		
	}
};

#endif
