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
 * Base structures:
 *   Atom, Box, MDClass
 * 
 *
*/

#ifndef _BASE_STRUCT_H_
#define _BASE_STRUCT_H_

#include <vector>
#include <omd/config.hpp>
#include <omd/toolkit.hpp>
#include <omd/paramhandler.hpp>

using std::vector;
using std::string;
using std::ostream;

#define ROOT 0

/**
 * @ingroup basestruct
 * @brief The structure to hold the system box information 
 * 
 * This structure keeps the coordinate and dimension of the box containing the
 * simulation atom. The coordinate is defined as following,
 * 
 *    - (x0,y0,z0) = (west,south,bottom)
 *    - (x1,y1,z1) = (east,north,top)
 * 
 * The dimension is stored as side lengths of the box (lx,ly,lz) and the half length
 * is kept in (hlx,hly,hlz).
 *
*/

class SysBox {
public:
	OMD_FLOAT x0, x1; ///< west--east
	OMD_FLOAT y0, y1; ///< south--north
	OMD_FLOAT z0, z1; ///< bottom--top
	OMD_FLOAT lx, ly, lz; ///< length
	OMD_FLOAT hlx, hly, hlz; ///< half length

	SysBox(){ // set box as undefined...
		x0=y0=z0=DBL_MAX;
		x1=y1=z1=-DBL_MAX;
		lx=ly=lz=hlx=hly=hlz=0.0;
	}

	bool undefined() {
		return x0==DBL_MAX&&y0==DBL_MAX&&z0==DBL_MAX&&x1==-DBL_MAX&&y1==-DBL_MAX&&z1==-DBL_MAX;
	}

};

enum stage_type {
		stage_undefined,
		stage_prepare, 
		stage_create,
		stage_init,
		stage_run,
		stage_finalize
};

/** 
 * @ingroup basestruct
 * @brief Base structure of an atom 
 * 
 * Atom class keeps all information of an atom, i.e. 
 * the position (x,y,z), velocity (vx,vy,vz), and force (fx,fy,fz).
 * The "potential" and "virial" variables are used to accumulate potential energy 
 * and the virial from the interatomic interaction, that is calculated by the ForceKernel.
 * These variables store double amount of the actual value due to double counting
 * in the force calculation loop. The calculation of total potential (virial) is 
 * performed by the main class (SimSystem), in function MeasurePotential().
 * 
 * Beside those variables, an atom has a general purpose auxilary array (aux).
 * This array can be used by the gadgets, by claiming one index of the array using
 * ClaimAuxVariable() function. The access the the aux variable is performed
 * via this index.
 * 
 * The "id" variable is used for identifying the type of an atom. This is also
 * used to determine the interaction force between atoms, and membership status 
 * of the container. The "xid" is reserved for multi purpose usage, i.e. 
 * the identification of group membership, atom tagging, etc.
 * 
 * The flag variable determine the status of an atom. The convention and the bit
 * position (value) and meaning is:
 *  - 1 --> 1=active  0=inactive
 *  - 2 --> 1=outside 0=inside box
 *  - 4 --> 1=ghost   0=real particle
 * 
 * These first three flags are claimed by MDSystem.
 * 
 */

class Atom {
public:
    char  id;              ///< atom's type id: interaction and container membership
    char  xid;             ///< extended id, tagging, grouping, etc
    int nid;			   ///< enumerated id
    int flag;            ///< the status and multi-purpose flag
    OMD_FLOAT x,   y,  z;
    OMD_FLOAT vx, vy, vz;
    OMD_FLOAT fx, fy, fz;
    OMD_FLOAT potential,virial;
    OMD_FLOAT aux[MAXAUXVAR];  ///< reserved multi-purpose vector
};

/**
 * @ingroup basestruct
 * @brief A multi-purpose linked list structure of integer (indices)
 *
 * This class is a garbage collecting 1D linked list of integer. It may be used
 * to book-keep the indices of atoms. The fetch() function is prepared to make
 * iteration through the elements. To put the iteration pointer at the begining
 * of the list reset() is used. fetch() will advance the pointer to the next
 * element, and when all element has been iterated a negative value is returned.
 *
 */

class IndexList: public MDToolkit {
public:
	struct ListElement{
		int value;
		ListElement* next;
	} *head, *tail, *iter;

	int length, current;

	IndexList() {
		head=tail=iter=NULL;
		length=current=0;
	}

	~IndexList() {
		while(head) {
			ListElement *le=head;
			head=head->next;
			delete le;
		}
	}

	// garbage collecting list...
	void push(int value) {

		if(head==NULL) {
			// new list...
			ListElement *newel=new ListElement;
			newel->next=NULL;
			newel->value=value;
			iter=tail=head=newel;
		} else {
			// head!=null, some memory is already allocated...
			if(tail==NULL){
				// list is empty, this is the first element...
				iter=tail=head;
				tail->value=value;
			} else {
				// there are already data in the list...
				if(tail->next==NULL) {
					// new element is appended to the tail...
					ListElement *newel=new ListElement;
					newel->next=NULL;
					newel->value=value;
					tail->next=newel;
					tail=tail->next;
				} else {
					// use garbage collection...
					tail=tail->next;
					tail->value=value;
				}
			}
		}

		length++;
	}

	void clear() {
		tail=iter=NULL;
		length=current=0;
	}

	void reset() {iter=head;current=0;}

	bool step()  {
		if(!iter) return false;
		if(current>=length) return false;
		iter=iter->next; current++;
		return true;
	}

	int fetch() {
		ListElement* a=iter;
		if(step()) return a->value;
		return -1;
	}

	void dump(string fname) {
		ofstream fl(fname.c_str());
		assert(fl.good(), "fail opening file "+fname+" for writing...");
		reset();
		int a;
		while((a=fetch())>=0) {
			fl << a << std::endl;
		}
		fl.close();
	}

};

// These flags are claimed by SimSystem
#define FLAG_ACTIVE  1
#define FLAG_OUTSIDE 2
#define FLAG_GHOST   4

/**
 * @ingroup baseclass
 * @brief Base class of all Object-MDs
 * 
 * MDClass is the base class of all Object-MDs classes.
 * It provides basic functions:
 *   - memory allocation
 *   - abstraction of general functions
 *   - class identity
 *   - managing multi-purpose parameters
 *   - atom checking and basic managements
 * 
 */

class MDClass:public MDToolkit {
	
protected:
	int mem_alloc_cnt[3]; // free, alloc, realloc counter

public:

	ParamHandler param;
	int id;
	MDClass();
	virtual ~MDClass(){}

	int* get_alloc_count() {return mem_alloc_cnt;}
	int  get_id() {return id;}

	virtual void mem_free(void*& ptr);
	virtual void mem_alloc(void*& ptr, size_t size);
	virtual void mem_realloc(void*& ptr, size_t size);	

	string           search_path(string path, string fl);
	virtual MDClass* set_id(int nid) {id=nid;return this;}

	virtual Atom& Atoms(int idx)=0;
	virtual Atom* AtomPtr(int idx)=0;
	virtual int GetNAtom()=0;
	virtual bool Check(){return true;}
	
	//---- Atom checking & management ---//
	bool CheckActive(int idx){return(Atoms(idx).flag&FLAG_ACTIVE);}
	bool CheckInside(int idx){return(!(Atoms(idx).flag&FLAG_OUTSIDE));}
	bool CheckGhost(int idx) {return((Atoms(idx).flag&FLAG_GHOST));}
	bool CheckFlag(int idx, const unsigned int f){return(Atoms(idx).flag&f);}

	void SetFlag(int idx, int f){Atoms(idx).flag|=f;}
	void UnsetFlag(int idx, int f){Atoms(idx).flag&=(~f);}
	void SetActive(int idx, bool a=true);
	void SetOutside(int idx, bool a=true);	
	void SetGhost(int idx, bool a=true);

	virtual void PrintInfo(ostream& ost);	
};


//-------------MEMORY-ALOCATION-MACROS-----------------------//
#define MemFree(PTR) mem_free((void*&)PTR)
#define MemAlloc(PTR,SIZE) mem_alloc((void*&)PTR,SIZE)
#define MemRealloc(PTR,SIZE) mem_realloc((void*&)PTR,SIZE)
#define MemNew(PTR, OBJTYPE) \
{PTR=new OBJTYPE;mem_alloc_cnt[1]++;logmem(PTR, sizeof(OBJTYPE)*SIZE);}
#define MemDelete(PTR) \
if(PTR){delete PTR;mem_alloc_cnt[0]++;logmem(PTR);} \
else{warn("trying to delete unallocated object");}
#define MemNewArray(PTR, OBJTYPE, SIZE) \
{PTR=new OBJTYPE[SIZE];mem_alloc_cnt[1]++;logmem(PTR, sizeof(OBJTYPE)*SIZE);}
#define MemDeleteArray(PTR) \
if(PTR){ delete[] PTR;mem_alloc_cnt[0]++;logmem(PTR);} \
else{warn("trying to delete unallocated array");}
//-------------------------------------------------------------//

/** 
 * @ingroup baseclass
 * @brief The primitive container of atoms
 * 
 * This class has duty to do basic memory management of atoms: allocating and
 * freeing memory. The class is flexible that is able to extend
 * or shrink the size of atom array. There are three types of AtomKeeper:
 *   - Storage, manages the atom data array and the index array.
 *   - Referral, does not manage the atom data. The index array points to atoms
 *     of other AtomKeeper.
 *   - Clone, both AtomArray and AtomIndex belong to other AtomKeeper.
 * 
 * The class uses garbage collector. The memory will not be released upon 
 * srinking and trucating the size, disposing or clearing the atom keeper. 
 * The pointers will be freed only when Release() is called explicitly.
 * To reserve an empty atom keeper the Allocate() function must be followed
 * by Clear(). Then atoms may be attached sequentially.
 * 
 * The basic functions of the class are the followings:
 *   - Expand, expands or shrinks the size of the container (atom array)
 *   - Append, appends another container to the tail of the array. If the type
 *     of the container is storage, the atoms of the other container are memory
 *     copied. Otherwise, the AtomIndex points continously to the AtomArray of 
 *     the appended container. A clone can not append.
 *   - Copy, copies atom structure from other container. Only a storage can do
 *     copy.
 *   - Cut, cuts the tail of the array.
 *   - Attach, attaches an atom to the array.
 *   - Dispose, dispose an atom from the array.
 * 
 */

class AtomKeeper:public MDClass {
	
	Atom* AtomArray;
	struct Ptr{Atom *p;} *AtomIndex;
	AtomKeeper* MyClone;
	int NBlock;
	int NAtom;
	int NAlloc;

public:

	enum KeeperType{Storage, Referral, Clone} Type;

	// FIXME! updating indices if required...
	IndexList IndexBook;

	//----Constructors/Destructor----//
	
	AtomKeeper(KeeperType kty=Storage);
	AtomKeeper(AtomKeeper& ak);	
	void Release();
	virtual ~AtomKeeper(){Release();}
	
	//----Atom acces methods----/
		
	Atom& operator[](int idx);	
	Atom& Atoms(int idx, MDClass *caller);
	Atom& Atoms(int idx) {return Atoms(idx,this);}	
	Atom* AtomPtr(int idx);
	
	//---- storage management methods ----//
	
	void Allocate(int sz, KeeperType t=Storage);
	void ChangeType(KeeperType tp);	
	bool GetType() {return Type;}
	void Expand(int sz);
	void Cut(int Offset) {Expand(Offset);}
	void Clear() {NAtom=0;}	
	
	void  Append(AtomKeeper& ak);
	void  Attach(Atom* pa);
	void  Attach(Atom& a){Attach(&a);}
	void  AssignByIndex(AtomKeeper& ak);
	void  Dispose(int idx);
	void  Copy(AtomKeeper& ak);
	void  CloneArray(AtomKeeper& ak);
	void  ReferArray(AtomKeeper& ak);
	Atom* GetArrayPtr(){return AtomArray;}
	Ptr*  GetIndexPtr(int offset=0){return(AtomIndex+offset);}
	int GetNAtom(){return NAtom;}
	int GetAllocSize(){return NAlloc;}	

};

#endif
