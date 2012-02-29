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

#ifndef _ATOM_GROUP_H_
#define _ATOM_GROUP_H_

#include <omd/container.h>
namespace omd {
  class AtomGroup: public AtomContainer {
    friend class MDSystem;
    class MDSystem* System;
    int group_flagmask;
    AtomContainer* source;
    AtomKeeper scratch_ak;
    
    // called only by MDSystem...
    void SetGroupMask(int bitpos){
      mdassert(bitpos<=32,"too many groups");
      group_flagmask=1<<bitpos;
    }
    
    void SyncAtomGroupMask();
    
  public:
    
    // key=tag -> a group member is tagged.
    // key=geometry -> membership is checked every step according
    //     to the defined range.
    enum GroupKey {tag, geometry} key;
    
    AtomGroup(string group_name, AtomContainer* ac, MDSystem* WorkSys);
    AtomGroup* Create();
    AtomGroup* Source(AtomContainer* a);
    AtomGroup* Union(AtomContainer& a);
    AtomGroup* Insert(Atom& a);
    AtomGroup* Select(int* idx_array);
    AtomGroup* SelectType(int type_id);
    AtomGroup* SelectGID(int gid);
    AtomGroup* SelectGT(double x, double y, double z); // greater than
    AtomGroup* SelectGE(double x, double y, double z); // greater or equal
    AtomGroup* SelectLT(double x, double y, double z); // less than
    AtomGroup* SelectLE(double x, double y, double z); // less or equal
    
    AtomGroup* SelectBox(double x0, double x1, 
                         double y0, double y1, 
                         double z0, double z1);
    
    AtomGroup* SelectInverseBox(double x0, double x1, 
                                double y0, double y1, 
                                double z0, double z1);
    
    AtomGroup* Commit(){return Create();}
    
    int  GetGroupMask() {return group_flagmask;}
    bool Member(Atom& a) {return (a.gid|group_flagmask);}
    
    double GetMass(int idx);
  };
}

#endif
