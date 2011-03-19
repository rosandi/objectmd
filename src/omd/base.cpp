/* omdlib
 ********************************************************
 *
 * ObjectMD
 * Molecular Dynamics Class Library
 * Version 3.0 (2009, 2011)
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
 * MDClass & AtomKeeper implementation
 *
*/

#include <omd/base.hpp>

//-------------------MD-CLASS----------------------//

MDClass::MDClass() {
	id=-1;
	set_name("MDCLASS");
	mem_alloc_cnt[0]=mem_alloc_cnt[1]=mem_alloc_cnt[2]=0;
	LoadEnv();
}

void MDClass::LoadEnv() {
	// load OMD environment variables
	char *p;
	if((p=getenv("OMD_LIB"))) param.push_pair("$OMD_LIB", p);
	else param.push_pair("$OMD_LIB", ".");
	if((p=getenv("OMD_CLASS"))) param.push_pair("$OMD_CLASS", p);
	else param.push_pair("$OMD_CLASS", ".");
	if((p=getenv("OMD_TABLE"))) param.push_pair("$OMD_TABLE", p);
	else param.push_pair("$OMD_TABLE", ".");	
}

/**
 * free memory
 * log entry: address
 */

void MDClass::mem_free(void*& ptr){
	if(!ptr){warn("trying to free unallocated memory");return;}
	free(ptr); mem_alloc_cnt[0]++;
	logmem(ptr);
}

/**
 * allocate memory
 * log entry: address size
 */

void MDClass::mem_alloc(void*& ptr, size_t size){
	void* p=NULL;p=malloc(size);mem_alloc_cnt[1]++;
	assert(p, "memory allocation error... (malloc)", get_name());
	ptr=p;logmem(p, size);
}

/** 
 * reallocate memory
 *  log entry: old_address new_address new_size
 */

void MDClass::mem_realloc(void*& ptr, size_t size){
	void *p=ptr;p=realloc(p, size);mem_alloc_cnt[2]++;
	assert(p, "memory allocation error... (realloc)", get_name());
	logmem(ptr,p,size);ptr=p;
}

/**
 * Check the existence of a filename and return the path to that file.
 * A file exists in current directory has the highest priority, otherwise
 * the function will search in the path stated by "path" parameter. If no file
 * is found, the program terminates.
 */

string MDClass::search_path(string path, string fl) {
	if(file_exist(fl)) return fl;
	if(param.exist(path)) {
		if(file_exist(param.string_value(path)+"/"+fl)) 
			return (param.string_value(path)+"/"+fl);
	}
	die("can not find file '"+fl+"' in current directory and in path "+path);
	return ""; // avoids warning..
}

void MDClass::PrintInfo(ostream& ost){
	ost <<"MDClass: "<<get_name()<<"(id:"<<id<<")\n";
}

//---- Atom checking & management ---//

void MDClass::SetActive(int idx, bool a){
	if(a)Atoms(idx).flag|=FLAG_ACTIVE;
	else Atoms(idx).flag&=(~FLAG_ACTIVE);
}

void MDClass::SetOutside(int idx, bool a){
	if(a)Atoms(idx).flag|=FLAG_OUTSIDE;
	else Atoms(idx).flag&=(~FLAG_OUTSIDE);
}

void MDClass::SetGhost(int idx, bool a){
	if(a)Atoms(idx).flag|=FLAG_GHOST;
	else Atoms(idx).flag&=(~FLAG_GHOST);
}


//------------------ATOM-KEEPER----------------------//

AtomKeeper::AtomKeeper(KeeperType kty){
	NAtom=0;
	NAlloc=0;
	Type=kty;
	AtomArray=NULL;
	AtomIndex=NULL;
	MyClone=NULL;
	set_name("storage");
}

/** 
 * Take an AtomKeeper and clone it. The type of new allocated AtomKeeper
 *  is a clone.
 */

AtomKeeper::AtomKeeper(AtomKeeper& ak){
	Type=Clone;
	AtomArray=ak.AtomArray;
	AtomIndex=ak.AtomIndex;
	NAtom=ak.NAtom;
	NAlloc=ak.NAlloc;
	MyClone=&ak;
	set_name("storage");
}


/**
 * If the structure is not a clone, the pointers will be released.
 * NAtom and NAlloc will be zeroed in any case. This function only
 * release the pointers if it is not NULL. After releasing all pointers
 * are reset to NULL.
 */

void AtomKeeper::Release(){
	switch(Type){
		case Clone: return;
		case Storage:
			if(AtomArray) MemFree(AtomArray);
			if(AtomIndex){
				MemFree(AtomIndex);
				AtomIndex=NULL;
			}
			break;
		case Referral:
			if(AtomIndex) MemFree(AtomIndex);
			break;
	}

	AtomArray=NULL;
	AtomIndex=NULL;
	NAtom=NAlloc=0;
}

/**
 * This operator give array like access to the elements.
 */

Atom& AtomKeeper::operator[](int idx) {
	if(idx>=NAtom)
		die("index out of range. index="+as_string(idx)+
			", n_atom="+as_string(GetNAtom()));

	return *(AtomIndex[idx].p);
}

/**
 * Returns reference to the atom pointed by idx. In the application
 * program, this function is called via the AtomContainer::Atoms() function.
 */

Atom& AtomKeeper::Atoms(int idx, MDClass *caller) {
	if(idx>=NAtom){
		die("index out of range. index="+as_string(idx)+
			", n_atom="+as_string(GetNAtom())+
			", called by "+caller->get_name());
	}
	return *(AtomIndex[idx].p);
}

/**
 * Returns the pointer of an atom.
 */

Atom* AtomKeeper::AtomPtr(int idx) {
	if(idx>=NAtom)
		die("index out of range. index="+as_string(idx)+
			", n_atom="+as_string(GetNAtom()));
	return AtomIndex[idx].p;
}

//---- storage management methods ----//

/**
 * Allocating memory of the container. A clone can not allocate memory.
 * The allocated atoms are active. This function will release first the
 * memories to make sure that they point nowhere. The safety of releasing
 * the pointer is guarantied by Release() function.
 * 
 * After allocation:
 *   - NAtom=NAlloc if type=Storage
 *   - NAtom=0 if type=Refferal
 * 
 */

void AtomKeeper::Allocate(int sz, KeeperType t){

	Type=t; Release();
	if(sz==0)return;
	
	Atom* a;
	switch(t) {
		case Storage:
			MemAlloc(AtomArray,sizeof(Atom)*sz);
			a=AtomArray;
			MemAlloc(AtomIndex,sizeof(Ptr)*sz);
			
			// Initialize as active ~outside ~ghost ~mirror
			for(int i=0;i<sz;i++) {
				a->id=a->xid=-1;a->flag=1;AtomIndex[i].p=a++;
			}

			NAlloc=NAtom=sz;
			break;
			
		case Referral:
			MemAlloc(AtomIndex,sizeof(Ptr)*sz);
			NAlloc=sz;
			NAtom=0;
			break;

		case Clone:
			die("attempt to allocate memory for a clone atom keeper");
	}
}

void AtomKeeper::ChangeType(KeeperType tp) {
	if(Type==tp) return;
	Type=tp;
	int n=NAtom;
	Release();
	if(tp!=Clone)Allocate(n,tp);
}

/**
 * @brief Expands the allocated storage
 * 
 * AtomArray will be expanded only if it is previously allocated (not NULL). 
 * Otherwise this function assumes AtomIndex points to anothed array, and expands only
 * AtomIndex. The function makes sure that the elements of AtomIndex points to the 
 * right position in the atom array.
 * 
 * Shrinking is possible by expanding to a smaller number than the previous size.
 * The garbage collection will be used when available. The variables NAtom and NAlloc 
 * are updated.
 * 
 * When the storage is expanded, the new atoms have undefined id and are unactive.
 * To prepare an empty storage (buffer), this function must be followed by Clear().
 */

/* todo: Use minimum expand size!
 *       expanding only one particle can make the program slow.
 * 
 */

void AtomKeeper::Expand(int sz){
	assert(Type!=Clone, "attempt to expand a clone atom keeper");
	int nold=NAtom;
	
	if(sz>NAlloc){ // reallocate storage
		MemRealloc(AtomIndex, sizeof(Ptr)*sz);
		assert(AtomIndex, "unable to reallocate atoms");
		blog("expanding storage from "+as_string(nold)+" to "+as_string(sz));

		// if storage: sync the atom index!
		if(Type==Storage){
			MemRealloc(AtomArray,sizeof(Atom)*sz);
			Atom* a=AtomArray;
			assert(AtomArray, "unable to reallocate atoms");
			for(int i=0;i<sz;i++)AtomIndex[i].p=a++;
			for(int i=nold;i<NAtom;i++)
				{AtomArray[i].id=-1;AtomArray[i].flag=0;}
		}
		NAlloc=NAtom=sz;
	} else { //use garbage collection
		NAtom=sz;
	}
}


/** 
 * @brief Append an Atom array to the structure
 * 
 * FIX THIS! If the structure is not a clone it copies the data 
 * and append it at the tail of the atom array. Otherwise this function 
 * only assign the Index array to points to the appended atom array (to be fixed!). 
 * 
 * If length is 0 or negative, the function will ignore it and return.
 * 
 */

void AtomKeeper::Append(AtomKeeper& ak) {
	if(ak.GetNAtom()<=0) return;
	int lastn=NAtom;
	int na=ak.GetNAtom();
	Expand(NAtom+na);
	Atom* a;
	switch(Type) {
		case Storage:
			a=AtomArray+lastn;
			memcpy(a,ak.GetArrayPtr(),sizeof(Atom)*na);
			break;
		case Referral:
			a=ak.GetArrayPtr();
			for(int i=lastn;i<NAtom;i++)AtomIndex[i].p=a++;
			break;
		case Clone:
			die("attempt to append to a clone");
	}
}

/**
 * @brief Attach one atom to the structure
 *
 * Attache an atom to the tail of a structure. 
 * If the sctructure is of Refferal type, the refference to the data will
 * be taken.
 * 
 */

void AtomKeeper::Attach(Atom* pa){
	if(Type==Clone) return;
	int ip=NAtom;
	Expand(NAtom+1);
	if(Type==Storage)Atoms(ip)=*pa;
	else AtomIndex[ip].p=pa; // referral
}

void AtomKeeper::Copy(AtomKeeper& ak) {
	assert(Type==Storage, "attempt to copy to a non storage keeper");
	int na=ak.GetNAtom();
	Expand(na);
	if(ak.GetType()==Storage) {
		memcpy(AtomArray, ak.GetArrayPtr(), na*(sizeof(Atom)));
	} else {
		for(int i=0;i<na;i++) AtomArray[i]=ak[i];
	}
}

void AtomKeeper::AssignByIndex(AtomKeeper& ak) {
	int a;
	Expand(IndexBook.length);
	Clear();
	IndexBook.reset();
	while((a=IndexBook.fetch())>=0) Attach(ak[a]);
}

/**
 * @brief Dispossing an atom from array
 *
 * The atom at index 'idx' will be replaced by the last atom in the
 * structure, then the number of atom decreases by one. Important: 
 * te method disregards the sequence of atoms in the structure!
 * The 'hole' of the atom pointer is always at the end of the structure,
 * and is kept as a garbage collection.
 * 
 */

void AtomKeeper::Dispose(int idx){
	if(NAtom==0||idx>=NAtom) return;
	NAtom--;
	AtomIndex[idx]=AtomIndex[NAtom];
}

void AtomKeeper::CloneArray(AtomKeeper& ak) {
	Release();
	Type=Clone;
	MyClone=&ak;
	AtomArray=ak.GetArrayPtr();
	AtomIndex=ak.GetIndexPtr();
	NAtom=ak.GetNAtom();
	NAlloc=ak.GetAllocSize();
}

void AtomKeeper::ReferArray(AtomKeeper& ak) {
	Release();
	Type=Referral;
	int sz=ak.GetNAtom();
	Expand(sz);Clear();
	for(int i=0;i<sz;i++) Attach(ak[i]);
}
