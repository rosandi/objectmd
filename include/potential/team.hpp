/* omdlib
 ********************************************************
 *
 * ObjectMD
 * Molecular Dynamics Class Library
 * Version 2.0 (2009)
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
 *       Table eam potential
 *
*/
 
#ifndef _T_EAM_H_
#define _T_EAM_H_

#include <vector>
#include <string>
#include <omd/forcekernel.hpp>
#include <omd/conditioner.hpp>
#include <omd/tablereader.hpp>

using std::vector;
using std::string;

// A rather overestimated length. But safe!
// Use even number of MAX_ALLOWED_SPECIES!!

#define MAX_ALLOWED_SPECIES 4

class TEmbedding: public Conditioner {
	friend class TForceEAM;
	int id_size;
	int NAlloc;
	int ED_idx;

	// REMEMBER!! Not all the atoms use EAM force!!! AtomID list is needed!
	int          AtomID[MAX_ALLOWED_SPECIES][MAX_ALLOWED_SPECIES]; 	
	TableReader  edens[MAX_ALLOWED_SPECIES*(MAX_ALLOWED_SPECIES/2 + 1)]; 
	TableReader  embed[MAX_ALLOWED_SPECIES];
	OMD_FLOAT       CutRadius[MAX_ALLOWED_SPECIES];
	public:
		
		TEmbedding();
		virtual ~TEmbedding() {}
		virtual void Init(MDSystem* WorkSys);
		
		virtual void AddTable(int, string);
		virtual void AddTable(int, int, string); //for different atom type

		OMD_FLOAT GetRho(Atom &at);
		OMD_FLOAT GetEmb(Atom &at);
		OMD_FLOAT GetRhoDeriv(OMD_FLOAT r, Atom &at);
		OMD_FLOAT GetEmbDeriv(Atom &at);
		
		virtual void Dump();
		virtual void IterationNode(Atom &at, Atom &to);
		virtual void PreCalculation();
		
		// Force modifier is used to do correction of potential
		virtual void ForceModifier();
		virtual void PrintInfo(ostream& ost);
};

/**
 * @brief Embedded Atom Method
 * 
 * This class is the force kernel for embedded atom method potential (EAM).
 * The potential is read from an eam table. Three tables is required:
 * PAIR, EDENS, and EMBED table. See @link functable the Object-MD table @endlink
 * chapter.
 * 
*/

class TForceEAM: public ForceKernel {
protected:
	string tablefile;
	TableReader phi;
	TEmbedding *emb;
public:
	TForceEAM(string PhiFile, TEmbedding *EM=NULL);

	TForceEAM* SetEmbeddingClass(TEmbedding *EM) {emb=EM;return this;} 
	TEmbedding* GetEmbeddingClass(){return emb;}
	virtual TForceEAM* SetTable(const char* PhiFile);

	virtual bool SearchAttachEmbeddingClass();
	virtual void Init(MDSystem* WorkSys); 
	virtual void ComputeHalf(Atom &at, Atom &to);
	virtual void PrintInfo(ostream& ost);
};

#endif
