// F means: friction-force

#ifndef _TTM_HOMOGEN_HPP_
#define _TTM_HOMOGEN_HPP_

#include <iostream>
#include <omd/comhandler.hpp>
#include <omd/conditioner.hpp>
#include <omd/tablereader.hpp>
#include <omd/paramhandler.hpp>
#include <detector/TempPressDetector.hpp>

/*
 * changes: using atom density instead of atom_volume
 */

/**
 * @ingroup conditioner
 * \brief Homogen Two Temperature Model to transfer energy 
 *        from electronic system to the lattice
 * 
 * No electron temperature diffusion.
 * The energy is given to the electronic system and transfered to the latice via
 * the coupling function. Used in laser ablation simulation.
 * 
*/

// explain: table file format, pulse function file, unity

class TTM_Homogen: public Force_Conditioner {
	
	CommunicationHandler* Communicator;
	string parfile;
	TableReader G;
	TableReader Ce;
	TableReader Ft;
	OMD_FLOAT E0;
	OMD_FLOAT G0;
	OMD_FLOAT C0;
	OMD_FLOAT pulse_width;
	OMD_FLOAT pulse_offset;
	OMD_FLOAT fd_time_step;
	OMD_FLOAT electron_temperature;
	OMD_FLOAT atom_volume;
	OMD_FLOAT step_dtemperature;
	OMD_FLOAT e_source; // timestep source 
	OMD_FLOAT de_energy; // timestep electron energy
	OMD_FLOAT da_energy; // timestep one atom energy
	OMD_FLOAT int_Een; // integrated (total) electron energy
	OMD_FLOAT int_Aen; // integrated (total) input energy to lattice
	OMD_FLOAT stop_time;
	
	TempPressDetector* temp_detector;
	bool balance_el_lattice;
	bool do_coupling;
	bool firstrun;
	bool firstfd;
	bool stopped;
	
	// for non homogen temperature distribution
	int    allocated_length;
	int    atemp;
	int    adens;
	
	enum {none, function, pulsetemp, pulse} source_type;
	
public:
	
	TTM_Homogen(string parf=""):Force_Conditioner(){
		set_name("Homogen two temperature model");
		register_class("TWO_TEMPERATURE_MODEL");
		parfile=parf;
		electron_temperature=0.0;
		e_source=0.0;
		int_Een=0.0;
		int_Aen=0.0;
		allocated_length=0;
		do_coupling=true;
		firstrun=true;
		firstfd=true;
		source_type=none;
		stopped=false;
	}
	
	virtual ~TTM_Homogen(){}
	
	void ReadParameter(){
		param.peek("energy_density", E0, -1.0); // energy / Volume (A^3)
		param.peek("temperature_input", E0, -1.0);
		param.peek("atom_volume", atom_volume, -1.0);
		param.peek("pulse_width", pulse_width, -1.0);
		param.peek("pulse_offset", pulse_offset, 0.0);
		param.peek("electron_init_temperature", electron_temperature, -1.0);
		param.peek("balance_electron_lattice", balance_el_lattice, false);
		param.peek("stop_ttm", stop_time, -1.0);
	}

	void CheckParameter(){
		blog("two-temperature-model: checking parameters...", LOGCREATE);
		if(param.exist("source_function_file")) source_type=function;
		else if(param.exist("temperature_input")) source_type=pulsetemp;
		else if(pulse_width>0.0) source_type=pulse;
	}

	// requires temperature and pressure detector
	void Init(MDSystem* WorkSys) {
		Force_Conditioner::Init(WorkSys);
		
		assert(WorkSys->type_of("simulation_system_grid"), "MDSystemGrid required!");
		Communicator=dynamic_cast<MDSystemGrid*>(WorkSys)->GetCommunicator();

		
		if(parfile!="") param.read(parfile);
		temp_detector=dynamic_cast<TempPressDetector*>(SearchGadgetType("tpdetector",false));
		assert(temp_detector, "TempPressDetector is needed by "+get_name());
		
		RegisterMessageSlot(new DataSlot("e_temp",2))
			->SetFormat("%0.3E")->SetData(electron_temperature);

		RegisterMessageSlot(new DataSlot("e_source",2))
				->SetFormat("%0.3E")->SetData(e_source);

		RegisterMessageSlot(new DataSlot("int_Een",2))
			->SetFormat("%0.5E")->SetData(int_Een); // total energy given from electron

		RegisterMessageSlot(new DataSlot("int_Aen",2))
		->SetFormat("%0.5E")->SetData(int_Aen); // total average energy given to lattice
		
		ReadParameter();
		CheckParameter();

		G.open(param.string_value("electron-phonon_file"), "ELECTRON_PHONON_COUPLING");
		G0=G.param.double_value("Constant");
		Ce.open(param.string_value("electron-phonon_file"), "ELECTRON_HEAT_CAPACITY");
		C0=Ce.param.double_value("Constant");
		
		if(source_type==function)
			Ft.open(param.string_value("source_function_file"));


		fd_time_step=WorkSys->Integrator->TimeStep;
		
		// finds the right detector
		temp_detector=NULL;
		for(int i=0;i<WorkSys->Detectors.size();i++){
			if(WorkSys->Detectors[i]->type_of("tpdetector")){
				temp_detector=dynamic_cast<TempPressDetector*>(WorkSys->Detectors[i]);
				temp_detector->SetIntensive(true);
				break;
			}
		}

		assert(temp_detector, 
		       "The two temperature modell class requires 'TempPressDetector'");
		
		temp_detector->SetIntensive();
		
		int_Een=0.0;
		assert(atom_volume>0.0, "define atom_volume!");
	
	}

	void PrintInfo(ostream& ost){
		// FIXME! check units!
		ost << "ID." << id << " " << get_name() << "\n";
		if(source_type==pulsetemp) ost << "temperature input = " << E0 << "K\n";
		else ost << "energy input = " << E0 << "eV/A^3\n";
		if(electron_temperature>0.0)
			ost<< "electron initial temperature = " << electron_temperature << "K\n";
		ost << "pulse offset = " << pulse_offset << " ps\n";
	}

private:

	void fd_main_loop(){
		int na=GetNAtom();

		OMD_FLOAT atomtemp=temp_detector->GetTemperatureAvg();
		
		// use constant volume!
		
		if(firstfd){
			if(electron_temperature<=0.0) electron_temperature=atomtemp;
			int_Een=0.5*electron_temperature*
			                        C0*Ce.read(electron_temperature)*atom_volume;
			int_Aen=(2.0*System->Kinetic)/System->GetTotalAtom();
			firstfd=false;
		}
		
		// Te
		OMD_FLOAT de=-G0*G.read(electron_temperature)*(electron_temperature-atomtemp)*fd_time_step;
		int_Een+=de*atom_volume;
		electron_temperature+=de/(C0*Ce.read(electron_temperature));
		
		
		// stop ttm if the electron temperature is lower than atom temperature
		if (electron_temperature<atomtemp) stopped=true;
		
	}

	void friction_force() {
		int na=GetNAtom();
		OMD_FLOAT *sumvx, *sumvy, *sumvz;
		temp_detector->GetCMVelocities(sumvx, sumvy, sumvz);
		OMD_FLOAT volcut=temp_detector->GetDetectVolume();
		OMD_FLOAT coup=G0*G.read(electron_temperature);
		
		// \Delta t_FD = \Delta t_MD
		OMD_FLOAT sum_eng=0.0;
		int neng=0;
		
		for(int i=0;i<na;i++){
			Atom* a=AtomPtr(i);

			if (a->flag&FLAG_GHOST) continue;
			if (!(a->flag&FLAG_ACTIVE)) continue;
			if (a->aux[adens]*volcut<=3.0) continue;

			OMD_FLOAT vcx=(a->vx-sumvx[i]);
			OMD_FLOAT vcy=(a->vy-sumvy[i]);
			OMD_FLOAT vcz=(a->vz-sumvz[i]);
			OMD_FLOAT vcsq=(vcx*vcx+vcy*vcy+vcz*vcz);
			
			if (vcsq==0.0) continue;
			
			OMD_FLOAT de=coup*(electron_temperature - a->aux[atemp])*atom_volume;
			OMD_FLOAT chi=de/vcsq;
			
			a->fx+=chi*vcx;
			a->fy+=chi*vcy;
			a->fz+=chi*vcz;
			
			sum_eng+=de*fd_time_step;
			neng++;

		}
		
		sum_eng=Communicator->TakeSUM(sum_eng);
		neng=Communicator->TakeSUM(neng);
		int_Aen+=(sum_eng/(OMD_FLOAT)neng);
		
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
	
	void ForceModifier() {

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
			friction_force();
		}
		
	}
};

#endif
