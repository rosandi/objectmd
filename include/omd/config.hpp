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
 *
 */

#ifndef _OMD_CONFIG_HPP_
#define _OMD_CONFIG_HPP_

#define VERSION_INDICATOR "OBJECT-MD V.03.00 CALLSIGN GLIDE"

#define MAXATOMTYPE     8
#define MAXAUXVAR       6
#define MAXPROC         512
#define DEFAULT_CONFIG_FILENAME "omd-parameter"
#define EDGE_TOLE       0.1
#define OMD_EPSILON     1e-20

// types
#define OMD_FLOAT double
#define OMD_INT   int
#define OMD_SIZET uint
#define OMD_CHAR  char

#define mdrseed() srand(getpid())
#define mdrand() ((OMD_FLOAT)rand()/(OMD_FLOAT)RAND_MAX)

#endif
