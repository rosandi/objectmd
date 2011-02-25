//----------------------Crystal FCC(111)----------------------//

#ifndef _CRYSTAL_FCC_111_HPP_
#define _CRYSTAL_FCC_111_HPP_

#include <omd/container.hpp>

/**
  @ingroup atom
  @brief FCC(111) Crystal
  
  The class to create a bulk of crystal.
  This class is initiated, by default, to create Platinum FCC(111) crystal.
  Another FCC(111) crystals can be made by providing the proper parameters
  to the constructor.
  To make another structure of crystal, one should overide Create 
  function (virtual).
  
  Its also capable to import crystal from a file with the following format,
  
  [x] [y] [z] [value]
  
  [value] field is arbitrary and will be ignored.
  
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
*/

class CrystalFCC111:public AtomContainer {
protected:
	OMD_INT xml, yml, zml;
	OMD_FLOAT lattice_constant;

public: 
	CrystalFCC111(OMD_INT XMLayer, OMD_INT YMLayer, OMD_INT ZMLayer, string mat_file);
	CrystalFCC111(const OMD_CHAR* fname, string mat_file);
	virtual AtomContainer* Create();
};

#endif
