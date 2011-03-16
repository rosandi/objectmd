//--------------FCC 110 Crystal------------------//

#include <crystal/FCC100.hpp>

AtomContainer* CrystalFCC100::Create() { 
  int     ATOMS_PER_UC=4;
  OMD_FLOAT  XUCOffset, YUCOffset, ZUCOffset,
		  XRelPos[ATOMS_PER_UC], YRelPos[ATOMS_PER_UC], ZRelPos[ATOMS_PER_UC],
		  XOrg, YOrg, ZOrg,
		  XPos, YPos, ZPos,
		  Deviation, hlc;

  int     XMLayerPerUC, YMLayerPerUC, ZMLayerPerUC,
		  XUCs, YUCs, ZUCs,
		  Count,
		  i, j, k, l;

	assert(xml&&yml&&zml, "crystal monolayers are not defined");

	hlc=lattice_constant/2.;
  
  /* Offset in X-, Y-, Z-Richtung in Einheiten der halben Gitterkonstanten */
	XUCOffset = 2.;    YUCOffset = 2.;    ZUCOffset = 2.;

  /* Anzahl der Monolagen in X-, Y-, Z-Richtung pro Einheitszelle          */
	XMLayerPerUC = 2;  YMLayerPerUC = 2;  ZMLayerPerUC = 2;

  /* Relativkoordinaten der vier Atome einer Einheitszelle                 */
	XRelPos[0] = 0.;   YRelPos[0] = 0.;   ZRelPos[0] = 0.;
	XRelPos[1] = 1.;   YRelPos[1] = 1.;   ZRelPos[1] = 0.;
	XRelPos[2] = 1.;   YRelPos[2] = 0.;   ZRelPos[2] = 1.;
	XRelPos[3] = 0.;   YRelPos[3] = 1.;   ZRelPos[3] = 1.;

  /*** scale with half lattice constant ***/
	for (i = 0; i < ATOMS_PER_UC; i++) {
		XRelPos[i] *= hlc;   YRelPos[i] *= hlc;   ZRelPos[i] *= hlc;
	}

	XUCOffset *= hlc;  YUCOffset *= hlc;  ZUCOffset *= hlc;

  /* Abstand der Monolagen                      */
	XMLDist = hlc; //XUCOffset / (double) XMLayerPerUC;
	YMLDist = hlc; //YUCOffset / (double) YMLayerPerUC;
	ZMLDist = hlc; //ZUCOffset / (double) ZMLayerPerUC;

	Box.x1=Box.x0+XMLDist*((OMD_FLOAT)xml-1.);
	Box.y1=Box.y0+YMLDist*((OMD_FLOAT)yml-1.);
	Box.z0=Box.z1-ZMLDist*((OMD_FLOAT)zml-1.);

  /* Anzahl der Einheitszellen in X-, Y- und Z-Richtung    */
	XUCs = (xml+XMLayerPerUC-1)/XMLayerPerUC;
	YUCs = (yml+YMLayerPerUC-1)/YMLayerPerUC;
	ZUCs = (zml+ZMLayerPerUC-1)/ZMLayerPerUC;

	Deviation = 0.0001 * hlc;
	// First check the crystal size
  
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
					if (    XPos <= Box.x1 + Deviation
						 && YPos <= Box.y1 + Deviation
						 && ZPos <= Box.z1 + Deviation ) 
						 Count++;
				}
			} 
		}
	} 

	int na=Count;
	Allocate(na);
	
	Count=0;
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
					if (   XPos <= Box.x1 + Deviation
						&& YPos <= Box.y1 + Deviation
						&& ZPos <= Box.z1 + Deviation ) {
						Atoms(Count).x = XPos;
						Atoms(Count).y = YPos;
						Atoms(Count).z = ZPos;
						Count++;
					}
				} 
			}
		}
	}

	CalcBox();
	created=true;
	return this;
}
