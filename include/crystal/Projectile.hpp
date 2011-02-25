
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <omd/container.hpp>

using std::ifstream;
using std::istringstream;

/**
 * @ingroup atom
 * @brief Projectile atom
 * 
 * The atom coordinate is the x,y coordinate of the target at z=0 surface.
 * The z coordinate of the projectile is defined as the height above the surface.
 * Theta is the incidence angle to the surface normal direction. and phi is the
 * direction angle measured from x axis, anti-clockwise.
 * 
 */

class Projectile:public AtomContainer {
	public:
	Projectile(string mat, OMD_FLOAT tx, OMD_FLOAT ty, OMD_FLOAT h,
               OMD_FLOAT theta, OMD_FLOAT phi, OMD_FLOAT energy)
	:AtomContainer(mat) {

		Allocate(1);		
		SetOutside(0); // put always outside box
		SetActive(0,true);

		OMD_FLOAT spi = M_PI/180.0;
		OMD_FLOAT vv=sqrt(2.0*energy/M);
		OMD_FLOAT vxy=vv*sin(theta*spi);
		
		Atoms(0).x=tx;
		Atoms(0).y=ty;
		Atoms(0).z=h;
		Atoms(0).vz = -vv*cos(theta*spi);
		Atoms(0).vx = vxy*cos(phi*spi);
		Atoms(0).vy = vxy*sin(phi*spi);

		OMD_FLOAT pxy =  Atoms(0).z*tan(theta*spi);
		OMD_FLOAT sx  =  pxy*cos(phi*spi);
		OMD_FLOAT sy  =  pxy*sin(phi*spi);		
		Atoms(0).x-=sx;
		Atoms(0).y-=sy;
		
	}

};
