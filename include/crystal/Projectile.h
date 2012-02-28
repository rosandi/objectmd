
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <omd/container.h>

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
	Projectile(string mat, double tx, double ty, double h,
               double theta, double phi, double energy)
	:AtomContainer(mat) {

		Allocate(1);		
		SetOutside(0); // put always outside box
		SetActive(0,true);

		double spi = M_PI/180.0;
		double vv=sqrt(2.0*energy/M);
		double vxy=vv*sin(theta*spi);
		
		Atoms(0).x=tx;
		Atoms(0).y=ty;
		Atoms(0).z=h;
		Atoms(0).vz = -vv*cos(theta*spi);
		Atoms(0).vx = vxy*cos(phi*spi);
		Atoms(0).vy = vxy*sin(phi*spi);

		double pxy =  Atoms(0).z*tan(theta*spi);
		double sx  =  pxy*cos(phi*spi);
		double sy  =  pxy*sin(phi*spi);		
		Atoms(0).x-=sx;
		Atoms(0).y-=sy;
		
	}

};
