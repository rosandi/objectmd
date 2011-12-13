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

#define VERSION_INDICATOR "OBJECT-MD V.03.02 CALLSIGN CHIRP"

#define MAXATOMTYPE     8
#define MAXAUXVAR       6
#define MAXPROC         512
#define DEFAULT_CONFIG_FILENAME "omd-parameter"
#define EDGE_TOLE       0.1
#define OMD_EPSILON     1e-20
#define OMD_PROXIMITY   1e-6

// types
#define OMD_FLOAT double

#define mdrseed() srand(getpid())
#define mdrand() ((OMD_FLOAT)rand()/(OMD_FLOAT)RAND_MAX)

/*
flag bits:
LOGMEMORY   1
LOGCREATE   2
LOGDESTROY  4
LOGINFO     8
LOGWARNING 16
*/

#define LOGFLAG (LOGCREATE|LOGDESTROY|LOGINFO|LOGWARNING)

#endif
