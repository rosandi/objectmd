//-------------------------Ackland Structure Detector -----------------------//

#include <cmath>
#include <algorithm>
#include <omd/iterator.hpp>
#include <detector/StructureDetector.hpp>

#define STRUCT_NON 0.0
#define STRUCT_FCC 1.0
#define STRUCT_BCC 2.0
#define STRUCT_HCP 3.0
#define STRUCT_ICO 4.0

#define X0  -1.0
#define X1  -0.945
#define X2  -0.915
#define X3  -0.755
#define X3A -0.705
#define X4  -0.195
#define X5   0.195
#define X6   0.245
#define X7   0.795
#define X7A  1.0

bool CmpDist(struct_neig a, struct_neig b){return (a.rsq<b.rsq);}

void StructureDetector::Init(MDSystem* WorkSys){
	DataDumper::Init(WorkSys);
	sqrcut=WorkSys->GetMaxCutRadius();
	sqrcut*=sqrcut;
	stidx=ClaimAuxVariable(true,"struct(F:1,B:2,H:3,I:4)","%1.0f"); // to store the structure
	nalloc=GetNAtom();
	MemNewArray(neig, vector<struct_neig>, nalloc);
}

void StructureDetector::IterationNode(OMD_SIZET at, OMD_SIZET to){
	OMD_FLOAT sqr=CalcSqrDistance(at,to);
	if(sqr<sqrcut){
		struct_neig sat={to,sqr};
		struct_neig sto={at,sqr};
		neig[at].push_back(sat);
		neig[to].push_back(sto);
	}
}

void StructureDetector::FindStructure(){
	OMD_INT na=GetNAtom();
	
	// reallocate if needed
	if(nalloc<na){
		nalloc=na;
		MemDeleteArray(neig);
		MemNewArray(neig, vector<struct_neig>, nalloc);
	}

	for(OMD_INT i=0;i<na;i++)neig[i].clear();
	Iterator->Iterate(this);
	
	for(OMD_INT i=0;i<na;i++){
		Atom* Ai=AtomPtr(i);
		
		sort(neig[i].begin(),neig[i].end(),CmpDist);
		OMD_FLOAT Rsq0=0.0;
		for(OMD_INT n=0;n<6;n++){
			if(n<(OMD_INT)neig[i].size()) Rsq0+=neig[i].at(n).rsq;
		}
		Rsq0/=6.0;
		
		OMD_FLOAT dist0=1.45*Rsq0;
		OMD_FLOAT dist1=1.55*Rsq0;
		vector<OMD_INT> neig_n0;
		
		OMD_INT n0,n1=0;
		for(OMD_INT n=0;n<(OMD_INT)neig[i].size();n++){
			if(neig[i].at(n).rsq<dist0)
				neig_n0.push_back(neig[i].at(n).idx);
			if(neig[i].at(n).rsq<dist1) n1++;
		}
		n0=neig_n0.size();
		
		// evaluate angles
		OMD_INT c[8];
		memset(c,0,sizeof(OMD_INT)*8);
		
		for(OMD_INT j=0;j<n0-1;j++){
			Atom* Aj=AtomPtr(neig_n0[j]);
			OMD_FLOAT dxj=Ai->x-Aj->x;
			OMD_FLOAT dyj=Ai->y-Aj->y;
			OMD_FLOAT dzj=Ai->z-Aj->z;
			OMD_FLOAT dj=sqrt(dxj*dxj+dyj*dyj+dzj*dzj);
			for(OMD_INT k=j+1;k<(OMD_INT)neig_n0.size();k++){
				Atom* Ak=AtomPtr(neig_n0[k]);
				OMD_FLOAT dxk=Ai->x-Ak->x;
				OMD_FLOAT dyk=Ai->y-Ak->y;
				OMD_FLOAT dzk=Ai->z-Ak->z;
				OMD_FLOAT dk=sqrt(dxk*dxk+dyk*dyk+dzk*dzk);
				OMD_FLOAT angle=(dxj*dxk+dyj*dyk+dzj*dzk)/dj/dk;
				
				// angle regions
				if     (X0<=angle&&angle<X1 )c[0]++;
				else if(X1<=angle&&angle<X2 )c[1]++;
				else if(X2<=angle&&angle<X3 )c[2]++;
				else if(X3<=angle&&angle<X3A)c[3]++;
				else if(X4<=angle&&angle<X5 )c[4]++;
				else if(X5<=angle&&angle<X6 )c[5]++;
				else if(X6<=angle&&angle<X7 )c[6]++;
				else if(X7<=angle&&angle<X7A)c[7]++;
			}
		}

// Gerolf's version: referred to lammps
		OMD_FLOAT dbcc=0.35*(OMD_FLOAT)c[4]/(OMD_FLOAT)(c[5]+c[6]-c[4]);
		OMD_FLOAT dcp =fabs(1.0-(OMD_FLOAT)c[6]/24.0);
		OMD_FLOAT dfcc=0.61*(OMD_FLOAT)(abs(c[0]+c[1]-6)+c[2])/6.0;
		OMD_FLOAT dhcp=(OMD_FLOAT)(abs(c[0]-3)+abs(c[0]+c[1]+c[2]+c[3]-9))/12.0;
		
		if(c[0]==7)dbcc=0.0;
		else if(c[0]==6) dfcc=0.0;
		else if(c[0]<=3) dhcp=0.0;

		OMD_FLOAT str=STRUCT_NON;			
		if(c[7]<=0.0) {
			if(c[4]<3.0){if(11<=n1&&n1<=13)str=STRUCT_ICO;}
			else if(dbcc<=dcp) {if(n1>=11)str=STRUCT_BCC;}
			else if(n1==12||n1==11){
				if(dfcc<dhcp)str=STRUCT_FCC;
				else str=STRUCT_HCP;
			}
		}
		Ai->aux[stidx]=str;
		//pos vector??
	}
}
