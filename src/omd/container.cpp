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
 *
 * Implementation of Crystal class
 *
*/

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <omd/base.hpp>
#include <omd/container.hpp>
#include <omd/forcekernel.hpp>
#include <omd/config.hpp>

using std::ios;
using std::string;
using std::ostream;
using std::ifstream;
using std::istringstream;

#define NOLEADINGSPACES(f) {\
	while(f.peek()==' ' || f.peek()=='\t' || f.peek()=='\n')f.get();\
}

AtomContainer::AtomContainer(AtomContainer& a, AtomContainer& b) {
	Value=0.0; 
	filename="";
	mat_file="";
	set_name("ATOM CONTAINER");
	register_class("atom_container");
	Value=0.0;
	posprec=valprec=5;
	MasterContainer=NULL;
	M=Z=-1.0;
	write_mode=0;
	Combine(a);
	Combine(b);
}

AtomContainer::AtomContainer(string material_file) {
	Value=0.0;
	filename="";
	mat_file="";
	set_name("ATOM CONTAINER");
	register_class("atom_container");
	Value=0.0;
	posprec=valprec=5;
	MasterContainer=NULL;
	write_mode=0;
	
	if(material_file!=""&&material_file!="-") ReadMaterial(material_file);

}

/**
 * Reads material definition file. The first read file has higher priority, so
 * the second read will just be ignored. This behavior is the property of
 * ParamHandler which reads the first parameter name occurance.
 */

AtomContainer* AtomContainer::ReadMaterial(string material_file) {
	assert(material_file!="", "empty material file");
	mat_file=material_file;	
	param.read(search_path("$OMD_TABLE", "def."+material_file));
	SetName(param.string_value("element"));
	blog("reading material: "+mat_file, LOGCREATE);

	try {
		M=param.double_value("mass")/NORMAL_AMU;
		Z=param.double_value("number");
	} catch(...) {
		throw string("incomplete material file: ")+
		      search_path("$OMD_TABLE", material_file);
	}
	
	return this;
}

MDClass* AtomContainer::set_id(int nid){
	id=nid;
	for (int i=0; i<GetNAtom(); i++) Atoms(i).id=(char)nid;
	return this;
}

AtomContainer* AtomContainer::SetXID(int nid){
	for (int i=0; i<GetNAtom(); i++) Atoms(i).xid=(char)nid;
	return this;
}

void AtomContainer::PrintInfo(ostream& ost) {
    ost  <<"id."<<id<<" "<<get_name()<<"\n"
         <<"number_of_atoms= "<<GetNAtom()<< "\n"
         <<"number="<<(int)Z<< "\n"
         <<"mass="<<M<<" (sim) "<<M*NORMAL_AMU<<"(amu)\n";
    if(CanImport()) ost<<"imported_from_file "<<filename<<'\n';
}

void AtomContainer::Release() {
		AtomStorage.Release();
}

/**
 * Shifts the atom positions to the desired position. 
 * dx, dy, and dz are the distance in the corresponding axis.
*/

AtomContainer* AtomContainer::Shift(OMD_FLOAT dx, OMD_FLOAT dy, OMD_FLOAT dz)
{
	if(!created) Create();
    for(int i=0;i<GetNAtom();i++) {
        Atoms(i).x+=dx;
        Atoms(i).y+=dy;
        Atoms(i).z+=dz;
    }
    CalcBox();
    return this;
}

/**
 * Allocates the atom keeper class. All atoms are zero initiated, and have 
 * the same id with the container.
*/

void AtomContainer::Allocate(int na, bool clear, AtomKeeper::KeeperType type)
{
	if(na){
		AtomStorage.Allocate(na,type);
		if (clear) {
    		for (int i=0;i<na;i++) {
        		Atoms(i).x=Atoms(i).vx=Atoms(i).fx=0.0;
        		Atoms(i).y=Atoms(i).vy=Atoms(i).fy=0.0;
        		Atoms(i).z=Atoms(i).vz=Atoms(i).fz=0.0;
        		Atoms(i).id=id;    
    		}
		}
	}
}

/**
 * Calculates the box dimension surrounding the AtomContainer, and returns the
 * SysBox object. The results is saved in variable Box.
 * The box boundary is calculated by taking the maximum and minimum atom coordinates. 
 * In the crystalline structures, this function must be called in the overriding 
 * descendant function and the offsets of periodic boundary must be added.
*/

SysBox& AtomContainer::CalcBox(){
	OMD_FLOAT MinX, MaxX, MinY, MaxY, MinZ, MaxZ;        
	int natom=GetNAtom();
	bool checkin=false;
	
	MinX=MinY=MinZ= DBL_MAX;
	MaxX=MaxY=MaxZ=-DBL_MAX;

	for(int i=0; i<natom; i++) {
		if(!CheckInside(i)) continue;
		checkin=true;
		// Minimums
		if (MinX>Atoms(i).x)MinX=Atoms(i).x;
		if (MinY>Atoms(i).y)MinY=Atoms(i).y;
		if (MinZ>Atoms(i).z)MinZ=Atoms(i).z;
    	// Maximums
		if (MaxX<Atoms(i).x)MaxX=Atoms(i).x;
		if (MaxY<Atoms(i).y)MaxY=Atoms(i).y;
		if (MaxZ<Atoms(i).z)MaxZ=Atoms(i).z;
	}
	
	if(!checkin) warn("no atom is inside simulation box");

	Box.x0=MinX;Box.x1=MaxX;
	Box.y0=MinY;Box.y1=MaxY;
	Box.z0=MinZ;Box.z1=MaxZ;

	Box.lx=fabs(Box.x1-Box.x0);
	Box.ly=fabs(Box.y1-Box.y0);
	Box.lz=fabs(Box.z1-Box.z0);
	
	Box.hlx=Box.lx/2.0;
	Box.hly=Box.ly/2.0;
	Box.hlz=Box.lz/2.0;
	
	return Box;
}

SysBox& AtomContainer::GetBox() {
	if(Box.undefined())
		return CalcBox();
	else
		return Box;
}

/**
 * Combining/appending a container. The container takes the parameters from
 * the material file of the first combined one.
 */

AtomContainer* AtomContainer::Combine(AtomContainer& a) {
	blog("Combining container: "+a.get_name(), LOGCREATE);
	if(!a.created)a.Create();
	set_name(a.get_name());

	if(a.M>0.0&&a.Z>0.0) {
		if(Z<0||M<0){ // if not initialized yet
			ReadMaterial(a.GetMaterialFile());
		} else assert(a.M==M&&a.Z==Z, "attempt to combine different type of atom");
	} else a.SetMZ(M,Z);

	AtomStorage.Append(a.GetAtomStorage());
	param.copy(a.param);

	if(Box.undefined())
		Box=a.Box;
	else {
		SysBox ab=a.GetBox();
		if(Box.x0>ab.x0)Box.x0=ab.x0;
		if(Box.y0>ab.y0)Box.y0=ab.y0;
		if(Box.z0>ab.z0)Box.z0=ab.z0;
		if(Box.x1<ab.x1)Box.x1=ab.x1;
		if(Box.y1<ab.y1)Box.y1=ab.y1;
		if(Box.z1<ab.z1)Box.z1=ab.z1;
		Box.lx=fabs(Box.x1-Box.x0);
		Box.ly=fabs(Box.y1-Box.y0);
		Box.lz=fabs(Box.z1-Box.z0);
		Box.hlx=Box.lx/2.0;
		Box.hly=Box.ly/2.0;
		Box.hlz=Box.lz/2.0;
	}
	created=true;
	return this;
}

AtomContainer* AtomContainer::copy_data(AtomContainer* a) {
	set_name(a->get_name());
	mat_file=a->mat_file;
	filename=a->filename;
	param.copy(a->param);
	id=a->id;
	M=a->M;
	Z=a->Z;
	Value=a->Value;
	posprec=a->posprec;
	valprec=a->valprec;
	write_mode=a->write_mode;
	set_logger(a->get_logger());
	return this;
}

/** 
 * Imports atoms data from a file. This function reads a plain text file
 * containing the fields of position, velocity, and some auxilary values. and some
 * pseudo-commands. A pseudo command begins with "#$" at the begining of the line.
 * If not defined, the format of the file is "x y z", and only this values will be read.
 * Otherwise, a pseudo-command "fields" must exist to define the fields.
 * The field definition must at least contains "x y z" for the position, and may
 * have "vx vy vz" for the velocities. The position of this characters, defines
 * the field number of the corresponding value.
 * The atomic properties may be determined using the following pseudo tags:
 * 
**/

AtomContainer* AtomContainer::Import(string fname, int aid) {

	string material_file;
	
	if(fname!="")filename=fname;
	
	assert(filename!="", "no file to import");
	assert(file_exist(filename), "import file name "+filename+" not found");
	
	ParamHandler p;
	p.read_pseudo(fname);
	
	if(p.peek("Material", material_file)){
		if(Z<0.0||M<0.0){
			ReadMaterial(material_file);
			warn("reading material data from "+material_file);
		} else warn("material properties was set before import");
	}

	int ix=0, iy=1, iz=2, ivx=-1, ivy=-1, ivz=-1, iid=-1, ixid=-1;
	if(p.exist("Fields")) {
		// take x,y,z and vx,vy,vz,xid fields if exist
		// fields line ends with "--"
		int idx=p.index_of("Fields")+1;
		for(int i=idx;i<p.size();i++){
			string sfl=p.string_value(i);	
			if(sfl=="x")ix=i-idx;
			if(sfl=="y")iy=i-idx;
			if(sfl=="z")iz=i-idx;
			if(sfl=="vx") ivx=i-idx;
			if(sfl=="vy") ivy=i-idx;
			if(sfl=="vz") ivz=i-idx;
			if(sfl=="id") iid=i-idx;
			if(sfl=="xid") ixid=i-idx;
			if(sfl=="--") break;
		}
	}
	
	if(aid>=0) assert(iid>=0, "id field not found");

	int it=0,lino=0;
	ifstream posf(filename.c_str());	
	assert(posf.good(), "failed to open file "+filename);

	vector<Atom> atom_register;

    while (posf.good()) {
	    char stbuff[512];
    	NOLEADINGSPACES(posf);
        posf.getline(stbuff, 511, '\n');
        if(posf.eof())break;

        lino++;
        if(stbuff[0]=='#') continue;
        
		vector<OMD_FLOAT> vline;
		istringstream ist(stbuff);
		OMD_FLOAT dd;
		Atom a;
		a.id=a.xid=-1;a.flag=1; //active
		a.x=a.y=a.z=a.vx=a.vy=a.vz=a.fx=a.fy=a.fz=0.0;

		try {
			while(ist>>dd) vline.push_back(dd);
						
			if(aid>=0&&vline.at(iid)!=aid) continue;
			
			a.x=vline.at(ix);
			a.y=vline.at(iy);
			a.z=vline.at(iz);
		
			if(ivx>=0&&ivy>=0&&ivz>=0) {
				a.vx=vline.at(ivx);
				a.vy=vline.at(ivy);
				a.vz=vline.at(ivz);
			}

			if(ixid>=0) a.xid=(int)vline.at(ixid);
			atom_register.push_back(a);

		} catch(...) {
			die("insufficient number of fields to import: "+filename+
				" line="+as_string(lino)+
				" particle_index="+as_string(it)+
				" number_of_fields="+as_string((int)(vline.size())));
		}

    }

    Allocate(atom_register.size());
    for(int i=0;i<(int)atom_register.size();i++) Atoms(i)=atom_register[i];

	created=true;
    posf.close();
    blog("read "+as_string(GetNAtom())+" atoms from "+filename, LOGCREATE);

    return this;
}

/**
* This function will dump atoms position data to a file reffered by 
* file name (fname). The mode parameter is used to control the writing mode.
* The default mode is "v", and is in the truncate mode.
* 
* To explain:
*   - This function dumps any atom keeper, unrelated to the owner class.
*   - if it has a master, then the master will do the dumping.
* 
*/

AtomContainer* AtomContainer::DumpAtoms(AtomKeeper& ak,
                                        string fname,
                                        int mode,
                                        bool*  AuxPrintable,
                                        char* AuxFormat[],
                                        string AuxNames)
{
	
	if(MasterContainer) {
		MasterContainer->DumpAtoms(ak,fname,mode,AuxPrintable,AuxFormat,AuxNames);
		return this;
	}

	ofstream f;
	if(!mode) mode=write_mode;

	if(mode&WM_APPEND) f.open(fname.c_str(), ios::app);
	else f.open(fname.c_str(), ios::trunc);
	
	assert(f.good(), "unable to create file "+fname);
	
	//-----Header stuff------//
	// header printed only if not WM_BARE and not WM_APPEND
	if(!(mode&WM_BARE)) {
		if (!(mode&WM_APPEND)) { 
			f 	<< "#####################################\n"
				<< "#  Object-MD data file\n"
				<< "#  (c) 2005 Rosandi\n#\n"
			    << "#$ Box "
			    << Box.x0<<" "<<Box.y0<<" "<<Box.z0<<" "
		        << Box.x1<<" "<<Box.y1<<" "<<Box.z1<<"\n";

			for (int i=0;i<(int)StringInfo.size();i++)f<<'#'<<StringInfo[i]<<'\n';
			
			if(mat_file!="")f<<"#$ Material "<<mat_file<<"\n";
		
			//---- field description ----//
			f << "#$ Fields x y z";

			if(AuxNames!="") f<<AuxNames;
			else if(AuxPrintable){
				for(int i=0;i<MAXAUXVAR;i++) 
					if(AuxPrintable[i])f << " aux"<<i;
			}
			
			if(mode&WM_VELOCITY) f<<" vx vy vz";
			if(mode&WM_FORCE)    f<<" fx fy fz";
			if(mode&WM_POTENTIAL)f<<" pot";
			if(mode&WM_VIRIAL)   f<<" vir";
			if(mode&WM_ID)       f<<" id";
			if(mode&WM_XID)      f<<" xid";
			if(mode&WM_NID)      f<<" nid";
								 f<<" --\n";
		}
	}
	//---------------------------//
	
	for (int i=0;i<ak.GetNAtom();i++) {
		
		if(!(mode&WM_GHOST)) if(ak.CheckGhost(i)) continue; // skip ghost atoms
		if(!ak.CheckActive(i)) if(!(mode&WM_DEADATOM)) continue; // skip dead atoms

		f << std::fixed << std::setprecision(posprec);
		f << ak[i].x << " "<< ak[i].y << " " << ak[i].z << " ";
		
		// Print aux variables
		if(AuxPrintable) {
		  for(int a=0;a<MAXAUXVAR;a++)
			if(AuxPrintable[a]) {
				if(AuxFormat) f<<as_string(ak[i].aux[a],AuxFormat[a])<<" ";
				else f<<as_string(ak[i].aux[a],"%0.5e")<<" ";
			}
		}
		f << std::scientific << std::setprecision(valprec);
		
		if(mode&WM_VELOCITY) f<<ak[i].vx<<" "<<ak[i].vy<<" "<<ak[i].vz<<" ";
		if(mode&WM_FORCE)    f<<ak[i].fx<<" "<<ak[i].fy<<" "<<ak[i].fz<<" ";
		if(mode&WM_POTENTIAL)f<<ak[i].potential<<" ";
		if(mode&WM_VIRIAL)   f<<ak[i].virial<<" ";
		if(mode&WM_ID)       f<<(int)ak[i].id<<" ";
		if(mode&WM_XID)      f<<(int)ak[i].xid<<" ";
		if(mode&WM_NID)      f<<(int)ak[i].nid<<" ";
                             f << "\n";
	}

	f.close();
	// gzip it!!
	if (mode&WM_ZIP) {
		string cmd("gzip "+fname);
		assert(system(cmd.c_str())!=-1, "cannot execute gzip program...");
	}

	return this;
}

/**
 * This function will write the atoms owned by the calling class.
 * The default values: name is the container class name, mode uses the class's
 * write_mode.
 * 
 */

AtomContainer* AtomContainer::DumpAtoms(
					string fname,
					int mode,
					bool* AuxPrintable,
					char* AuxFormat[],
					string AuxNames)
{
	string dumpname(get_name());
	if(fname!="")dumpname.assign(fname);
	return DumpAtoms(AtomStorage,dumpname,mode,AuxPrintable,AuxFormat,AuxNames);
}

#define WRITEST(STR, FL) { \
	int sz=(STR).size()+1; \
	char longst[sz]; \
	memset(longst,0,sz); \
	(STR).copy(longst, sz-1); \
	fwrite(&sz, sizeof(int), 1, FL); \
	fwrite(longst, sizeof(char), sz, FL); \
}

/**
 * Save all containers data to binary file. THe mode parameter is the file
 * writing mode:
 *   - "w", overwrite
 *   - "a", append
 * 
 * Inside the binary file, the data block fo one container starts with "$ATOM:"
 * string followed by its name of 32 character length. Binary header:
 * 
 *  All string tags are written in 32 character long.
 */

AtomContainer* AtomContainer::Save(string binname, string mode) {

	FILE* fl = fopen(binname.c_str(), mode.c_str());

	char cname[32];
	
	assert(fl, "cannot open file for writing "+binname);
	
	memset(cname,0,32);
	strcpy(cname,"$ATOMCONTAINER:");
	fwrite(cname, sizeof(char), 32, fl);

	WRITEST(get_name(), fl);
	WRITEST(mat_file, fl);
	WRITEST(filename, fl);
	WRITEST(param.raw_string(), fl);

	fwrite(&id, sizeof(int), 1, fl);
	fwrite(&M, sizeof(OMD_FLOAT), 1, fl);
	fwrite(&Z, sizeof(OMD_FLOAT), 1, fl);
	fwrite(&Value, sizeof(OMD_FLOAT), 1, fl);
	fwrite(&posprec, sizeof(int), 1,fl);
	fwrite(&valprec, sizeof(int), 1, fl);
	fwrite(&write_mode, sizeof(int), 1, fl);

	int na=GetNAtom();
	fwrite(&na, sizeof(int), 1, fl);

	for(int i=0;i<na;i++) {
		int dr=fwrite(&(Atoms(i)), sizeof(Atom), 1, fl);
		assert(dr, "failure in writting file "+binname);
	}
	
	assert(!ferror(fl), "error writing to file "+binname);
	fclose(fl);

	return this;
}

// strings are previously saved ending with \0

#define READST(STR, FL) { \
	int sz; \
	fread(&sz, sizeof(int), 1, FL); \
	char longst[sz]; \
	fread(longst, sizeof(char), sz, FL); \
	(STR).assign(longst); \
}

/**
 * Load data from previously saved binary file. This function search first the
 * starting mark of a container block, and compare the "blockname" to the
 * read name in the binary file after "$ATOMCONTAINER:" mark.
 * If no "blockname" is defined the class name will be used.
 */

AtomContainer* AtomContainer::Load(string binname, string blockname) {
	char cname[32];
	
	FILE* fl = fopen(binname.c_str(), "r");
	assert(fl,"LOAD: can not open file '"+binname+"' for reading");
	if(blockname=="")blockname.assign(get_name());
	assert(blockname!="", "undefined block name to read binary from "+binname);
	
	bool found=false;

	while(!found) {
		int pos;
		char ch=0x00;
		while((ch=getc(fl))!='$') if(feof(fl)) break;
		pos=ftell(fl);
		memset(cname, 0, 32);
		fread(cname, sizeof(char), 31, fl); // minus '$'
		if(feof(fl)) break;

		if(string("ATOMCONTAINER:")==cname){
			string ss;
			READST(ss, fl);
			if(blockname==ss) {found=true;break;}
		}

		fseek(fl, pos, SEEK_SET);
	}

	assert(found, "can not find atom data block '"+blockname+"' inside file '"+binname+"'");
	set_name(blockname);
	param.clear();
	string stake;

	READST(mat_file, fl);
	READST(filename, fl);
	READST(stake, fl); param.assign(stake);

	fread(&id, sizeof(int), 1, fl);
	fread(&M, sizeof(OMD_FLOAT), 1, fl);
	fread(&Z, sizeof(OMD_FLOAT), 1, fl);
	fread(&Value, sizeof(OMD_FLOAT), 1, fl);
	fread(&posprec, sizeof(int), 1,fl);
	fread(&valprec, sizeof(int), 1, fl);
	fread(&write_mode, sizeof(int), 1, fl);
	
	// read atom
	int na;
	fread(&na, sizeof(int), 1, fl);	

	blog("Loaded data length="+as_string(na)+" blocks "+
	    as_string((int)(na*sizeof(Atom)))+" bytes", LOGCREATE);

	Allocate(na, false);
	int dr=0;
	for(int i=0;i<na;i++) dr+=fread(&(Atoms(i)), sizeof(Atom), 1, fl);
	assert(dr==na, "Fail in reading file '"+binname+"'");
	assert(!ferror(fl), "error loading file '"+binname+"'");
	fclose(fl);
	blog("Loaded: "+get_name()+"@"+binname, LOGCREATE);

	return this;
}

/**
 * @brief Give kinetic energy to atoms in the container.
 * 
 * This function gives uniform distributed velocity to all atoms in the
 * container. The energy is given in energy unit per atom.
 * 
 */

AtomContainer* AtomContainer::SetKineticEnergy(OMD_FLOAT ek_per_atom) {
	if(!created) Create();
	int natom=GetNAtom();

	if(ek_per_atom<=0.0) {
		for(int i=0;i<natom;i++) {
			Atoms(i).vx=0.0;
			Atoms(i).vy=0.0;
			Atoms(i).vz=0.0;
		}
		return this;
	}

	mdrseed();

	// Give same energy, random direction
	OMD_FLOAT svx=0.0,svy=0.0,svz=0.0;
	for(int i=0;i<natom;i++) {
		OMD_FLOAT v=sqrt(2.0*ek_per_atom/GetMass(i));
		OMD_FLOAT sx=mdrand()*v;
		OMD_FLOAT sy=mdrand()*(v-sx);
		OMD_FLOAT sz=v-(sx+sy);

		svx+=(Atoms(i).vx=(mdrand()>0.5)?-sqrt(sx):sqrt(sx));
		svy+=(Atoms(i).vy=(mdrand()>0.5)?-sqrt(sy):sqrt(sy));
		svz+=(Atoms(i).vz=(mdrand()>0.5)?-sqrt(sz):sqrt(sz));

	}
	
	svx/=(OMD_FLOAT)natom;
	svy/=(OMD_FLOAT)natom;
	svz/=(OMD_FLOAT)natom;
	
	for(int i=0;i<natom;i++) {
		Atoms(i).vx-=svx;
		Atoms(i).vy-=svy;
		Atoms(i).vz-=svz;
	}

	OMD_FLOAT checkek=0.0;
	for(int i=0;i<natom;i++) {
		checkek+=0.5*GetMass(i)*(
			Atoms(i).vx*Atoms(i).vx+
			Atoms(i).vy*Atoms(i).vy+
			Atoms(i).vz*Atoms(i).vz);
	}
	
	OMD_FLOAT fact=sqrt((ek_per_atom*(OMD_FLOAT)natom)/checkek);
	
	for(int i=0;i<natom;i++) {
		Atoms(i).vx*=fact;
		Atoms(i).vy*=fact;
		Atoms(i).vz*=fact;
	}

	return this;
}

/** 
 * @brief Set kinetic temperature
 * 
 * This function gives kinetic temperature to the atom container.
 * The initial temperature is twice of set value.
 * Considering all atoms at are in their equilibrium position, the system
 * needs one or two pico seconds to reach at the desired temperature with
 * moderate fluctuation. 
 * 
 */

AtomContainer* AtomContainer::SetTemperature(OMD_FLOAT temperature) {
    return SetKineticEnergy(3.0*temperature*MD_BOLTZMANN);
}
