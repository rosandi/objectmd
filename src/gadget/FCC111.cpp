//----------------------Crystal FCC(111)----------------------//

#include <crystal/FCC111.hpp>

//------------------Crystal-------------------//

CrystalFCC111::CrystalFCC111(OMD_INT XMLayer, OMD_INT YMLayer, OMD_INT ZMLayer, string mat_file)
:AtomContainer(mat_file) {
    
    if(param.exist("lattice_constant"))
      lattice_constant=param.double_value("lattice_constant");
    else die("lattice_constant needed in material file: "+mat_file);
    
    xml=XMLayer;
    yml=YMLayer;
    zml=ZMLayer;
	Box.x0=Box.x1=Box.y0=Box.y1=Box.z0=Box.z1=0.0;
}

/**
 This function creates FCC(111) lattice. This
 code is taken from impact code, AG Urbassek (Thomas J. Colla).
 The function takes monolayer numbers, x, y, z, respectivelly, and the
 lattice constant.
 
 The boundary box of the crystal is not defined here.
 The caller object must take responsible for it. 
*/

// FIXME! This function is ugly... must be replaced!

AtomContainer* CrystalFCC111::Create()
{

#define S2                sqrt(2.)
#define S6                sqrt(6.)
#define V2S3              (2. * sqrt(3.))
#define V01S6             (1./6. * sqrt(6.))
#define V03S6             (1./3. * sqrt(6.))
#define V05S2             (0.5 * sqrt(2.))
#define V05S6             (0.5 * sqrt(6.))
#define V06S3             (2./3. * sqrt(3.))
#define V06S6             (2./3. * sqrt(6.))
#define V08S6             (5./6. * sqrt(6.))
#define V13S3             (4./3. * sqrt(3.))
    
    OMD_INT     ATOMS_PER_UC=6;
     
    OMD_FLOAT  XUCOffset, YUCOffset, ZUCOffset,
		XRelPos[ATOMS_PER_UC], YRelPos[ATOMS_PER_UC], ZRelPos[ATOMS_PER_UC],
		XOrg, YOrg, ZOrg, XPos, YPos, ZPos, Deviation, hlc;

    OMD_INT     XMLayerPerUC, YMLayerPerUC, ZMLayerPerUC,
            XUCs, YUCs, ZUCs, Count, i, j, k, l;
	
	assert(xml&&yml&&zml, "crystal monolayers are not defined");
        
    XMLDist= lattice_constant*sqrt(6.)/12.;
    YMLDist= lattice_constant*sqrt(2.)/4.;
    ZMLDist= lattice_constant/sqrt(3.);
    hlc    = lattice_constant/2.;
        
    /* Offset in X-, Y-, Z-Richtung in Einheiten der halben Gitterkonstanten */
    XUCOffset = S6;    YUCOffset =S2;    ZUCOffset = V2S3;
    /* Anzahl der Monolagen in X-, Y-, Z-Richtung pro Einheitszelle          */
    XMLayerPerUC = 6;  YMLayerPerUC = 2;  ZMLayerPerUC = 3;

    /* Relativkoordinaten der vier Atome einer Einheitszelle                 */
    XRelPos[0] = 0.;      YRelPos[0] = 0.;      ZRelPos[0] = 0.;
    XRelPos[1] = V05S6;   YRelPos[1] = V05S2;   ZRelPos[1] = 0.;
    XRelPos[2] = V06S6;   YRelPos[2] = 0.;      ZRelPos[2] = V06S3;
    XRelPos[3] = V01S6;   YRelPos[3] = V05S2;   ZRelPos[3] = V06S3;
    XRelPos[4] = V03S6;   YRelPos[4] = 0.;      ZRelPos[4] = V13S3;
    XRelPos[5] = V08S6;   YRelPos[5] = V05S2;   ZRelPos[5] = V13S3;

    // scale with half lattice constant
    for (i = 0; i < ATOMS_PER_UC; i++) {
        XRelPos[i] *= hlc;   YRelPos[i] *= hlc;   ZRelPos[i] *= hlc;
    }
 
    XUCOffset *= hlc;  YUCOffset *= hlc;  ZUCOffset *= hlc;
    /* Abstand der Monolagen                      */ 
    XMLDist = XUCOffset / (OMD_FLOAT) XMLayerPerUC;
    YMLDist = YUCOffset / (OMD_FLOAT) YMLayerPerUC;
    ZMLDist = ZUCOffset / (OMD_FLOAT) ZMLayerPerUC;
  	
    // Obergrenzen in X-, Y- und Z-Richtung
    // Positif Z direction is open surface
    // top surface @ z=0

    Box.x1=Box.x0+XMLDist*((OMD_FLOAT)xml-1.);
    Box.y1=Box.y0+YMLDist*((OMD_FLOAT)yml-1.);
    Box.z0=Box.z1-ZMLDist*((OMD_FLOAT)zml-1.);

    /* Anzahl der Einheitszellen in X-, Y- und Z-Richtung    */
    XUCs = (xml+XMLayerPerUC-1)/XMLayerPerUC;
    YUCs = (yml+YMLayerPerUC-1)/YMLayerPerUC;
    ZUCs = (zml+ZMLayerPerUC-1)/ZMLayerPerUC;
    
    Deviation = 0.0001 * hlc;
    
    /* First check crystal size */
    Count = 0;
    for (i = 0; i < ZUCs; i++) {
        ZOrg = (OMD_FLOAT) i * ZUCOffset + Box.z0;
        for (j = 0; j < YUCs; j++) {
            YOrg = (OMD_FLOAT) j * YUCOffset + Box.y0;
            for (k = 0; k < XUCs; k++) {
                XOrg = (OMD_FLOAT) k * XUCOffset + Box.x0;

                for (l = 0; l < ATOMS_PER_UC; l++) {
                    XPos = XOrg + XRelPos[l];
                    YPos = YOrg + YRelPos[l];
                    ZPos = ZOrg + ZRelPos[l];
                	if(XPos<=Box.x1+Deviation&&
                	   YPos<=Box.y1+Deviation&&
                	   ZPos<=Box.z1+Deviation)Count++;
                }
            }
        }
    }

    OMD_INT na = Count;
    Allocate(na);
    Count=0;
    /* Then create crystal */
    for (i = 0; i < ZUCs; i++) {
        ZOrg = (OMD_FLOAT) i * ZUCOffset + Box.z0;
        for (j = 0; j < YUCs; j++) {
            YOrg = (OMD_FLOAT) j * YUCOffset + Box.y0;
            for (k = 0; k < XUCs; k++) {
                XOrg = (OMD_FLOAT) k * XUCOffset + Box.x0;
                for (l = 0; l < ATOMS_PER_UC; l++) {
                    XPos = XOrg + XRelPos[l];
                    YPos = YOrg + YRelPos[l];
                    ZPos = ZOrg + ZRelPos[l];
                	if (XPos<=Box.x1+Deviation&&
                	    YPos<=Box.y1+Deviation&&
                	    ZPos<=Box.z1+Deviation){
                        Atoms(Count).x = XPos;
                        Atoms(Count).y = YPos;
                        Atoms(Count).z = ZPos;
                        Count++;
                    }
                }
            }
        }
    }
    created=true;
    return this;
}

SysBox& CrystallFCC111::CalcBox() {
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
}

