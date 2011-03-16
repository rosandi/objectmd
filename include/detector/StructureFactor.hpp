/**
 * @ingroup detector
 * @brief Structure Factor
 * 
 * This class calculates the structure factor of a crystal. The fields of the
 * output file are:
 * 
 * <pre>
 *   {r} {pair correlation, g(r)} {q} {structure function, s(q)}
 * </pre>
 *
 * where, r is radius and q is a function of incidence angle and wave length: 
 * \f$q=4\pi sin(\theta)/\lambda\f$
 * 
 * This class should be used in STATIC_MODE due to long range iteration 
 * in radial distribution calculation.
 * 
 * ref: Zhibin Lin, L. V. Zhigilei, Phys. Rev. B, 2006, 73, 184113
 */
 
#include <detector/DataDumper.hpp>

class StructureFactor: public DataDumper {
	OMD_FLOAT wbin;
	OMD_FLOAT rho;
	int nbin;
	OMD_FLOAT qstep;
	OMD_FLOAT delta, rmax;
	OMD_FLOAT *gr, *sq;
	OMD_FLOAT z_up_limit;
	OMD_FLOAT z_down_limit;
	bool enumerate;
	int  contrib_na;
	SysBox box;

public:

	StructureFactor(
		string outf,
		bool enumerate_outf,
		int ibin=1024,
		OMD_FLOAT iq=10.0,
		OMD_FLOAT r_max=-1.0,
		OMD_FLOAT detect_time=0.0):
	DataDumper("",detect_time, outf) {
		enumerate=enumerate_outf;
		nbin=ibin;
		rmax=r_max;
		qstep=iq/(OMD_FLOAT)ibin; // starts from 1
	}

	virtual ~StructureFactor() {MemFree(gr);MemFree(sq);}

	void calc_box() {
		OMD_FLOAT MinX, MaxX, MinY, MaxY, MinZ, MaxZ;        
		int natom=GetNAtom();
		
		MinX=MinY=MinZ= DBL_MAX;
		MaxX=MaxY=MaxZ=-DBL_MAX;

		for(int i=0; i<natom; i++) {
			// Minimums
			if (MinX>Atoms(i).x)MinX=Atoms(i).x;
			if (MinY>Atoms(i).y)MinY=Atoms(i).y;
			if (MinZ>Atoms(i).z)MinZ=Atoms(i).z;
			// Maximums
			if (MaxX<Atoms(i).x)MaxX=Atoms(i).x;
			if (MaxY<Atoms(i).y)MaxY=Atoms(i).y;
			if (MaxZ<Atoms(i).z)MaxZ=Atoms(i).z;
		}

		box.x0=MinX;box.x1=MaxX;
		box.y0=MinY;box.y1=MaxY;
		box.z0=MinZ;box.z1=MaxZ;

		box.lx=fabs(box.x1-box.x0);
		box.ly=fabs(box.y1-box.y0);
		box.lz=fabs(box.z1-box.z0);
		
		box.hlx=box.lx/2.0;
		box.hly=box.ly/2.0;
		box.hlz=box.lz/2.0;
	}
	
	string GetFilename(){
		if(enumerate) return (FilenamePrefix+Filename+FilenamePostfix);
		return (FilenamePrefix+Filename);
	}

	void Init(MDSystem* WorkSys) {
		DataDumper::Init(WorkSys);
		calc_box();
		if(rmax<0.0||rmax>box.hlx) rmax=box.hlx;
		wbin = rmax/(OMD_FLOAT)nbin;
		rho = (OMD_FLOAT)GetNAtom()/(box.lx*box.ly*box.lz);
		delta = 2.0*M_PI/rmax;
		MemAlloc(gr, nbin*sizeof(OMD_FLOAT));
		MemAlloc(sq, nbin*sizeof(OMD_FLOAT));
		
		if(param.exist("skip_surface")) {
			OMD_FLOAT dsurf=param.double_value("skip_surface");
			z_down_limit=box.z0+dsurf;
			z_up_limit=box.z1-dsurf;
			blog("skiping "+as_string(dsurf)+" from surface");
		} else {
			z_down_limit=0.0;
			z_up_limit=0.0;
		}
		
		contrib_na=0;
		for(int i=0;i<GetNAtom();i++) {
			if((Atoms(i).z<=z_up_limit)&&
			   (Atoms(i).z>=z_down_limit))contrib_na++;
		}
		
		if(System->GetMode()!=STATIC_MODE)
			warn("using structure factor detector in non STATIC_MODE");

	}
	
	// the contribution of an atom is defined by its z position compared
	// to up and down limit (measured from surface).
	// the free surfaces must be in z direction.
	
	// FIXME! parallelize this...
	
	void calc_distribution() {
		int na=GetNAtom();
		for (int i=0; i<nbin; i++) gr[i]=0.0;
		char prog[]="...........oooooooOOOOOOOOOOooooooo...........:::::::* **";
		int progc=0;

		for (int i=0; i<na-1; i++) {
			
			int contrib=0;
			OMD_FLOAT z=Atoms(i).z;
			contrib=((z<=z_up_limit)&&(z>=z_down_limit))?1:0;

			for (int j=i+1; j<na; j++) {
				OMD_FLOAT r=CalcDistance(i, j);
				if (r<rmax) {
					
					z=Atoms(j).z;
					int ig = (int)floor(r/wbin);

					gr[ig]+=(((z<=z_up_limit)&&(z>=z_down_limit))?contrib+1:contrib);
					
				}
			}
			if(progc>strlen(prog))progc=0;
			std::cerr << prog[progc++]<<"\r";
		}
		std::cerr<<" \r";

	}

	void calc_struc_factor() {
		int na=GetNAtom();
		OMD_FLOAT pi_rmax=M_PI/rmax;
		for(int q=0;q<nbin;q++) {
			OMD_FLOAT qq=(OMD_FLOAT)q*qstep;
			sq[q]=0.0;
			for(int i=0;i<nbin;i++) {
				OMD_FLOAT r=(OMD_FLOAT)i*wbin;
				if((r==0.0)||(qq==0.0))	sq[q]=0.0; else
				sq[q]+=gr[i]*(sin(qq*r)/(qq*r))*(sin(pi_rmax*r)/(pi_rmax*r));
			}
			sq[q]*=(wbin/(OMD_FLOAT)contrib_na);
			sq[q]+=1.0;
		}
	}

	void Measure() {		
		calc_box();
		calc_distribution();
		calc_struc_factor();
		int na=GetNAtom();
		ofstream fout(GetFilename().c_str());
		for(int i=0; i<nbin; i++) {
			int j=i+1;
			OMD_FLOAT nv=(4./3.)*M_PI*(OMD_FLOAT)((j*j*j)-(i*i*i))*wbin*wbin*wbin;
			gr[i]/=(nv*contrib_na*rho);
			fout << wbin*(OMD_FLOAT)i<<" "<<gr[i]<<" "<<qstep*(OMD_FLOAT)i<<" "<<sq[i]<<"\n";
		}
		fout.close();
	}
};
