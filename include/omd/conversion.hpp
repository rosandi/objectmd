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
 * A collection of physical conversion units
 *
*/

#ifndef _CONVERSION_H_
#define _CONVERSION_H_

#define PLANK         6.62606896e-34
#define PLANKBAR      1.054571628e-34
#define BOLTZMANN     1.3806504e-23 // Joule/Kelvin
#define NORMAL_AMU    9.64853397731809e+3 // E_CHARGE*N_AVO*A/ps
#define E_CHARGE      1.602176487e-19
#define VACUUM_EPS    8.854187817e-12
#define N_AVO         6.02214179e+26
#define LIGHT_SPEED   2.99792458e8
#define AMU_IN_KG     1.66053886e-27
#define EVOLT_KELVIN  1.1604505e+04 // 1eV = x Kelvin
#define EVOLT_JOULE   E_CHARGE      // 1eV = x Joule

// EVOLT_JOULE/Angstrom^3/giga: 1.6022e-19/1e-30/1e9
#define MD_GIGAPASCAL 160.22
#define MD_MASS_IN_KG(MASS) (MASS*E_CHARGE*1e-4)

// MD_ prefix is used for constant that can be
// directly used in the simulation unity

// (BOLTZMANN/E_CHARGE) eV/K
#define MD_BOLTZMANN  8.61734279e-05

#define F_CONST       VACUUM_EPS

// eV.ps
#define MD_PLANK 4.13566733e-3
#endif
