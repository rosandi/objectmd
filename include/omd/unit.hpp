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
 * ObjectMD header file
 *
 *     Unit conversion handler
 *
*/

#ifndef _MD_UNIT_H_
#define _MD_UNIT_H_

#include <cmath>
#include <string>
#include <omd/config.hpp>
#include <omd/conversion.hpp>

using std::string;

/**
 * @ingroup tool
 * @brief Unit Conversion Tool
 *
 * The concept of MDUnit is to make unit conversion easier. An application program
 * synchronizes its unit upon attachment of every MDGadget (AddConditioner/AddDetector),
 * by calling SetUnit() function. In the initialization stage an MDGadget rechecks again 
 * the Unit pointer. If it is unassigned, then the Unit of the main simulation system 
 * will be taken. This mechanism guaranties the unit is available and synchronized 
 * for the initialization.
 * 
 * To use another unit system this class must be inherited. The conversion 
 * factors must be initialized and the conversion methods must be overrided 
 * as needed. The new unit system must be attached to the simulation system class
 * before Initiate() function is invoked (before Run() or FirstRun()).
 * By default this class is initialized using Atomic MKS
 * (amu, eV, Angstrom, picosecond).
 * 
 * A physical quantity can be given to the simulation system using the provided
 * conversion methods with In prefix, like InAngstrom(), In PicoSecond(), etc.
 * To define how a number is converted in a different unit system, these 
 * functions must be overided. To GET a value from the unit system, simply use 
 * functions with the coresponding measurement name. e.g. Length().
 * For a collective measurements, like temperature and pressure, the value is converted
 * first by the "caller" function with the number of atoms. i.e. for temperature the conversion
 * will look like: Temperature(E_kinetic/Number_of_Atom).
 * 
 * The functions Format<Measurement>() are reserved for giving the system value in
 * a printable format (value+unit).
 * 
 */


class MDUnit {

protected:
	OMD_FLOAT fac_length;
	OMD_FLOAT fac_mass;
	OMD_FLOAT fac_force;
	OMD_FLOAT fac_energy;
	OMD_FLOAT fac_temperature;
	OMD_FLOAT fac_pressure;
	OMD_FLOAT fac_time;

public:
	string st_length;
	string st_mass;
	string st_force;
	string st_energy;
	string st_temperature;
	string st_pressure;
	string st_time;
	string default_format;
	
	void SetFactor(OMD_FLOAT ftime, OMD_FLOAT flength, OMD_FLOAT fmass, 
                   OMD_FLOAT fforce, OMD_FLOAT fenergy, 
                   OMD_FLOAT ftemp, OMD_FLOAT fpress)
	{
		fac_time=ftime;
		fac_length=flength;
		fac_mass=fmass;
		fac_force=fforce;
		fac_energy=fenergy;
		fac_temperature=ftemp;
		fac_pressure=fpress;
	}

	void SetUnit(const char* stime, const char* slength, const char* smass, 
                 const char* sforce, const char* senergy, 
                 const char* stemp, const char* spress)
	{
		st_length.assign(slength);
		st_mass.assign(smass);
		st_force.assign(sforce);
		st_energy.assign(senergy);
		st_temperature.assign(stemp);
		st_pressure.assign(spress);
	}

	// to be implemented
	void Set(const char* unitstr){}
	
	MDUnit(const char* unitstr=NULL) {
		default_format.assign("%0.5e");
		if(unitstr)Set(unitstr);
		else {
			SetFactor(
				1.0,                  // time
				1.0,                  // length
				NORMAL_AMU,           // mass
				1.0,                  // force
				1.0,                  // energy
				2.0/3.0/MD_BOLTZMANN, // temperature
				MD_GIGAPASCAL         // pressure
			);
			SetUnit("ps","Angstrom","amu","N","eV","K","GPa");
		}
	}
	
	virtual ~MDUnit(){}
	
	// from simulation to the defined units
	virtual OMD_FLOAT Time(OMD_FLOAT v){return fac_time*v;}
	virtual OMD_FLOAT Length(OMD_FLOAT v){return fac_length*v;}
	virtual OMD_FLOAT Mass(OMD_FLOAT v){return fac_mass*v;}
	virtual OMD_FLOAT Force(OMD_FLOAT v){return fac_force*v;}
	virtual OMD_FLOAT Energy(OMD_FLOAT v){return fac_energy*v;}
	virtual OMD_FLOAT Temperature(OMD_FLOAT v){return fac_temperature*v;}
	virtual OMD_FLOAT Pressure(OMD_FLOAT v){return fac_pressure*v;}
	
	// Give a value to simulation.. Inverse conversion
	// read: give data in PicoSecond to the system, etc.
	virtual OMD_FLOAT InPicoSecond(OMD_FLOAT v){return v/fac_time;}
	virtual OMD_FLOAT InAngstrom(OMD_FLOAT v){return v/fac_length;}
	virtual OMD_FLOAT InAmu(OMD_FLOAT v){return v/fac_mass;}
	virtual OMD_FLOAT InNewton(OMD_FLOAT v){return v/fac_force;}
	virtual OMD_FLOAT InEvolt(OMD_FLOAT v){return v/fac_energy;}
	virtual OMD_FLOAT InKelvin(OMD_FLOAT v){return v/fac_temperature;}
	virtual OMD_FLOAT InGigaPascal(OMD_FLOAT v){return v/fac_pressure;}

	
	virtual string Format(OMD_FLOAT v, const char* stformat=NULL){
		char st[128];
		if(stformat)sprintf(st,stformat,v);
		else sprintf(st,default_format.c_str(),v);
		return st;
	}
	
	void SetFormat(const char* stformat)
	{default_format.assign(stformat);}
	
	virtual string FormatTime(OMD_FLOAT v, const char* stformat=NULL)
	{return Format(Time(v),stformat)+" "+st_time;}

	virtual string FormatLength(OMD_FLOAT v, const char* stformat=NULL)
	{return Format(Length(v),stformat)+" "+st_length;}

	virtual string FormatMass(OMD_FLOAT v, const char* stformat=NULL)
	{return Format(Mass(v),stformat)+" "+st_mass;}

	virtual string FormatForce(OMD_FLOAT v, const char* stformat=NULL)
	{return Format(Force(v),stformat)+" "+st_force;}
	
	virtual string FormatEnergy(OMD_FLOAT v, const char* stformat=NULL)
	{return Format(Energy(v),stformat)+" "+st_energy;}

	virtual string FormatTemperature(OMD_FLOAT v, const char* stformat=NULL)
	{return Format(Temperature(v),stformat)+" "+st_temperature;}

	virtual string FormatPressure(OMD_FLOAT v, const char* stformat=NULL)
	{return Format(Pressure(v),stformat)+" "+st_pressure;}
	
	virtual string FormatVolume(OMD_FLOAT v, const char* stformat=NULL)
	{return Format(v*fac_length*fac_length*fac_length,stformat)+" "+st_length+"^3";}

	virtual string FormatArea(OMD_FLOAT v, const char* stformat=NULL)
	{return Format(v*fac_length*fac_length,stformat)+" "+st_length+"^2";}
	
};

#endif

