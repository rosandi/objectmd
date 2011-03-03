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
 * AtomContainer: Base class of all atom containers
 *
*/

#ifndef _CONTAINER_H_
#define _CONTAINER_H_

#include <fstream>
#include <omd/base.hpp>
#include <omd/conversion.hpp>

using std::ostream;
using std::ofstream;

// write modes
#define WM_APPEND     1
#define WM_VELOCITY   2
#define WM_FORCE      4
#define WM_DEADATOM   8
#define WM_BARE      16
#define WM_GHOST     32
#define WM_POTENTIAL 64
#define WM_VIRIAL    128
#define WM_ID        256
#define WM_XID       512
#define WM_NID       1024
#define WM_ZIP       2048

/**
 * @ingroup essential
 * @brief Base class of all atomic structures
 * 
 * AtomContainer creates, manages and manipulates atoms in the simulation box.
 * The simulation manager MDSystem is a descendant of this class.
 * The task of AtomContainer is managing the storage of simulation atoms: 
 * save, load, import, and dump the data to file. 
 * This class reserves an unimplemented method: Create(). This
 * method can be reimplemented in a child class to create any structure of atoms
 * (crystals). 
 * 
 * Atom data is stored in an AtomStorage object. The atoms can be accessed using
 * Atoms() and AtomPtr() function. To fetch the number of atoms, use GetNAtom().
 * 
 * Initially the class sets the mass and atomic number to negative,
 * which means "undefined". These variables is checked by the simulation 
 * application in the initialization stage. The application must first set the
 * values, which are the material properties, prior to this stage. This is
 * normally done by the constructor, by passing a material name to it (see
 * @link create_atom Creating Simulation Atoms @endlink section).
 * 
 * A container can be a refferal to parts of other container. The other
 * container is called the master container, as an example: MDSystem. 
 * In this case the container is not the owner of the atoms reffered by it.
 * The DumpAtoms() function is redirected the the master container.
 * 
*/

class AtomContainer: public MDClass {
protected:
	AtomKeeper AtomStorage;
	AtomContainer* MasterContainer;

// variables to save ---------------------
	string mat_file;
	OMD_INT posprec,valprec;
	OMD_INT write_mode;

public:	
    OMD_FLOAT   M;
    OMD_FLOAT   Z;
    OMD_FLOAT   Value;
	SysBox   Box;
	vector<string> StringInfo;
	string    filename;
//-----------------------------------------

public:

	bool created;

    //-------Methods---------//
    
	AtomContainer() { // empty container
		Value=0.0; 
		filename="";
		mat_file="";
		set_name("atom container");
		M=Z=-1.0;
		Value=0.0;
		posprec=valprec=5;
		MasterContainer=NULL;
		write_mode=0;
		created=false;
	}

	AtomContainer(string material_file); // initiate with material properties from file
	AtomContainer(AtomContainer& a, AtomContainer& b); // copy from two containers

	virtual ~AtomContainer() {blog("cleaning up",LOGDESTROY);}
	
	void set_logger(MDToolkit* mdlog){
		logger=mdlog;
		AtomStorage.set_logger(mdlog);
	}

	AtomContainer* copy_data(AtomContainer* a);

	/** Access an atom data from the atom structure. This function may be used
	 *  both as left or right value on an expression **/
	virtual Atom& Atoms(OMD_INT idx){
		return AtomStorage.Atoms(idx,this);
	}
	
	/** Take the pointer of an atom **/
	virtual Atom* AtomPtr(OMD_INT idx){
		return AtomStorage.AtomPtr(idx);
	}
	
	virtual MDClass* set_id(OMD_INT nid);
	virtual AtomContainer*       SetID(OMD_INT nid){set_id(nid); return this;}
	virtual AtomContainer*       SetXID(OMD_INT nid);

	/** Sets the name of the AtomContainer. By default a newly created 
	 *  class is named "AtomContainer" **/
	AtomContainer* SetName(string nm){
		set_name(nm);
		AtomStorage.set_name("storage@"+get_name());
		return this;
	}

	/** Assigning the master. e.g. used when unificating atoms;
	 *  the simulation main class is the master of all unificated containers **/	
	AtomContainer* SetMaster(AtomContainer* master){
		MasterContainer=master;
		return this;
	}
	
	/** Reads material file **/
	AtomContainer* ReadMaterial(string material_file);
	
	/** set the write mode in DumpAtoms() function.**/
	AtomContainer* SetWriteMode(OMD_INT wm){write_mode|=wm;return this;}
	AtomContainer* UnsetWriteMode(OMD_INT wm){write_mode&=(~wm);return this;}
	AtomContainer* ClearWriteMode(){write_mode=0;return this;}
	OMD_INT GetWriteMode(){return write_mode;}

	/** Get the number of atoms stored in or referred by the container **/
	virtual OMD_SIZET GetNAtom(){return AtomStorage.GetNAtom();}

	// calling compatibility...
	virtual OMD_FLOAT GetMass(OMD_SIZET idx){return M;}
	virtual OMD_FLOAT GetZ(OMD_SIZET idx){return Z;}

	/** Take the reference to the primitive atom keeper, AtomKeeper **/
	virtual AtomKeeper& GetAtomStorage(){return AtomStorage;}
	
	/** Releases the atom storage by calling the AtomKeeper::Release() **/
	void           Release();

	/** prints the information of the class **/
 	virtual void           PrintInfo(ostream& ost);
	virtual void           Allocate(OMD_INT na, bool clear=true, AtomKeeper::KeeperType type=AtomKeeper::Storage);
	
	virtual AtomContainer* Create(){
		if(CanImport())Import(filename);
		created=true;
		return this;
	}	
	
    virtual AtomContainer* Import(string fname,OMD_INT aid=-1);
    
	virtual AtomContainer* DumpAtoms(AtomKeeper& ak, 
									 string fname, 
									 OMD_INT mode=0, 
									 bool* AuxPrintable=NULL, 
									 OMD_CHAR* AuxFormat[]=NULL, 
									 string AuxNames="");
									
	virtual AtomContainer* DumpAtoms(string fname="",
	                                 OMD_INT mode=0, 
	                                 bool* AuxPrintable=NULL, 
	                                 OMD_CHAR* AuxFormat[]=NULL, 
	                                 string AuxNames="");
	                                 
	AtomContainer*         Combine(AtomContainer& a);
	AtomContainer*         Combine(AtomContainer& a, AtomContainer& b);
	
	virtual AtomContainer* Save(string binname, string mode="w");
	virtual AtomContainer* Load(string binname, string blockname="");
	
	AtomContainer*         Shift(OMD_FLOAT dx, OMD_FLOAT dy, OMD_FLOAT dz);
	AtomContainer*         AddVelocity(OMD_FLOAT vx, OMD_FLOAT vy, OMD_FLOAT vz){
								for(OMD_SIZET i=0;i<GetNAtom();i++){
									Atoms(i).vx+=vx;
									Atoms(i).vy+=vy;
									Atoms(i).vz+=vz;
								}
								return this;
							}

	AtomContainer*         SetMZ(OMD_FLOAT mass, OMD_FLOAT za) {M=mass/NORMAL_AMU; Z=za; return this;}
	AtomContainer*         SetFilename(string fname) {filename.assign(fname);return this;}
	AtomContainer*         SetKineticEnergy(OMD_FLOAT ek_per_atom);
	AtomContainer*         SetTemperature(OMD_FLOAT temperature);
	bool                   CanImport() {return (filename!="" && filename!="-");}
	SysBox& CalcBox();    
	SysBox& GetBox(){return Box;}
	AtomContainer* Deactivate() {for(OMD_SIZET i=0; i<GetNAtom(); i++)SetActive(i,false);return this;}
	AtomContainer* Activate() {for (OMD_SIZET i=0;i<GetNAtom();i++)SetActive(i,true); return this;}
	AtomContainer* ActivateAtom(OMD_INT idx){SetActive(idx,true); return this;}
	AtomContainer* DeactivateAtom(OMD_INT idx){SetActive(idx,false);return this;}
	AtomContainer* SetPrecision(OMD_INT prec, OMD_INT vprec=0){posprec=prec;valprec=vprec?vprec:prec;return this;}
	
    AtomContainer* PopInfo(OMD_INT cnt=1){while(0<cnt--)StringInfo.pop_back();return this;}
    AtomContainer* PushInfo(string sp){StringInfo.push_back(sp);return this;}
    
    string GetMaterialFile(){return mat_file;}

	/** To be implemented **/
	/*
	AtomContainer*         Center(){die("not implemented yet");return this;}
	AtomContainer*         Rotate(){die("not implemented yet");return this;}
	AtomContainer*         RestorePosition(){die("not implemented yet");return this;}
    */
};

#endif
