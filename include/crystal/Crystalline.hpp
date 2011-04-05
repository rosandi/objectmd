//----------------------Crystal FCC(111)----------------------//

#ifndef _CRYSTALLINE_HPP_
#define _CRYSTALLINE_HPP_

#include <omd/container.hpp>

/**
  @ingroup atom
  @brief Crystalline structure
  
  Geometry (right-hand system):
 
   <pre>
   xy:                        xz: 
             north   x1,y1               top      x1,z1
          +---------+                +-----------+
          |         |                |           |
   west   |   (.)   | east      west |    (x)    | east
          |   +z    |                |    +y     |
          +---------+                +-----------+
      x0,y0  south               x0,z0  bottom
   </pre>
 
	The orientation of the crystal defines the +z crystalline surface.
*/

class Crystalline:public AtomContainer {
protected:
	string orientation;
	int xml, yml, zml;
	OMD_FLOAT lattice_constant;
	OMD_FLOAT XMLDist, YMLDist, ZMLDist;

public: 
	Crystalline(string orien, int XMLayer, int YMLayer, int ZMLayer, string mat_file):
	AtomContainer(mat_file) {
		set_name("CRYSTALLINE");
		register_class(get_name());
		assert(param.exist("lattice_constant"), 
			   "lattice_constant needed in material definition file: "+mat_file);
		lattice_constant=param.double_value("lattice_constant");
		orientation=orien;
		xml=XMLayer;
		yml=YMLayer;
		zml=ZMLayer;
		XMLDist=YMLDist=ZMLDist=-1.0;
	}
	
	virtual AtomContainer* Create()=0;
	
	virtual SysBox& CalcBox(){
		
		assert(XMLDist>0.0&&YMLDist>0.0&&ZMLDist>0.0, 
			   "the monolayer distances is required for a crystalline structure. "
			   "these variables should be defined in Create() function.");
		
		AtomContainer::CalcBox();

		Box.x0-=0.5*XMLDist;
		Box.y0-=0.5*YMLDist;
		Box.z0-=0.5*ZMLDist;
		Box.x1+=0.5*XMLDist;
		Box.y1+=0.5*YMLDist;
		Box.z1+=0.5*ZMLDist;
		Box.lx=fabs(Box.x1-Box.x0);
		Box.ly=fabs(Box.y1-Box.y0);
		Box.lz=fabs(Box.z1-Box.z0);	
		Box.hlx=Box.lx/2.0;
		Box.hly=Box.ly/2.0;
		Box.hlz=Box.lz/2.0;
		return Box;		
	}
};


#endif
