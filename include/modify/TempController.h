
//-----------------------TempController------------------------//

#ifndef _TEMP_CONTROLLER_HPP_
#define _TEMP_CONTROLLER_HPP_

#include <omd/modify.h>
#include <detect/ThermoDetector.h>

namespace omd {


/**
 * @ingroup modify
 * @brief Temperature controller
 *
 * The class that controls crystal temperature by scalling the velocity 
 * according to the set temperature (Kelvin) which is stored as the first parameter 
 * of the double parameter. With DoRangeScaling() function user can choose 
 * the scalling scheme to be "all atom scalling" or "range scalling".
 * All atom scalling means the scalling is applied to all atoms, while
 * range scalling only do scalling to the atoms in the border ranges, defined
 * by SetRange(); 
 *
 * The scalling factor is calculated as follows,
 *
 * \f[ \lambda=\sqrt {1 + \frac{dt}{\tau} 
 *     \left( \frac{T_{set}}{T_{act} } - 1 \right )} \f]
 *
 * where \f$T_{set}\f$ and \f$T_{act}\f$ are the set temperature and 
 * the actual temperature, respectivelly.
 * 
 * The temperature is taken from TempPressDetector. If the detector is not
 * present the system kinetic energy will be used.
 * 
*/

class TempController: public PostModify {
	double TempSet, Tau;
	ThermoDetector* temp_detector;

public:
	
	TempController(double temp, double tau, string target="") {	
		TargetName=target;
		TempSet=temp;
		Tau=tau;
		temp_detector=NULL;
		set_name("tcontrol");
		mdassert(Tau>0.0, "zero coupling factor!");
	}
	
	void SetTemperature(double temp){TempSet=temp;}
		
	void Init(MDSystem *WorkSys) {
		PostModify::Init(WorkSys);		
		temp_detector=dynamic_cast<ThermoDetector*>(SearchGadgetType("thermo detector",false));
		if(!temp_detector) warn("using systems Kinetic variable for temperature measurement");
	}

	void PrintInfo(ostream& ost){
		ost << "ID." << id << " " << get_name()
			<< " -- Temperature=" << TempSet
			<< " Coupling=" << Tau << " temperature measurement: ";
		if(temp_detector) 
			ost <<"'"<<temp_detector->get_name()<<"'";
		else
			ost << " system Kinetic variable";
		
		ost<< "\n";
	}

	virtual double GetSystemTemperature() {
		if(temp_detector) return temp_detector->GetTemperature();
		return Unit->Temperature(System->Kinetic/System->GetTotalAtom());
	}

	void PostIntegration() {
		if(System->Step==0) return;
		double stemp=GetSystemTemperature();
		if(stemp<=0.0) return;
		double Factor=sqrt(1.+GetTimeStep()/Tau*(TempSet/stemp-1.));
		int na=GetNAtom();
		for(int i=0;i<na;i++) {
			if(CheckActive(i)) {
				Atoms(i).vx*=Factor;
				Atoms(i).vy*=Factor; 
				Atoms(i).vz*=Factor;
			}
		}
	}

};

}

#endif

