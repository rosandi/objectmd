//-----------------------FreeAtom----------------------//

#ifndef _FREE_ATOM_HPP_
#define _FREE_ATOM_HPP_

#include <omd/container.hpp>

/**
 * @ingroup atom
 * @brief The class of one free atom
 *
 * This class is a single atom container. It will be usefull in creating
 * single atom with different type, such as projectile. The constructor
 * accept coordinate, and velocity of the atom. User must provide mass
 * and atomic number, as well.
*/

class FreeAtom:public AtomContainer {
public:
    //------------Methods---------------//
    FreeAtom(OMD_FLOAT x,  OMD_FLOAT y,  OMD_FLOAT z, 
             OMD_FLOAT vx, OMD_FLOAT vy, OMD_FLOAT vz,
             string mat_file)
    :AtomContainer(mat_file)
    {                    
        Allocate(1);
        Atoms(0).x=x;
        Atoms(0).y=y;
        Atoms(0).z=z;
        Atoms(0).vx=vx;
        Atoms(0).vy=vy;
        Atoms(0).vz=vz;
        Atoms(0).fx=Atoms(0).fy=Atoms(0).fz=0.0;
        created=true;
    }
};

#endif
