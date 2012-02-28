#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <climits>
#include <muParser/muParser.h>
#include <omd/omdtool.h>
#include <misc/fdtd.h>
// For 2D problems...

using std::vector;
using std::ifstream;
using std::ofstream;
using std::istringstream;
using namespace omd;

// equation parser
FDTDSolver* myptr;
double var_x=0.0;
double var_y=0.0;
double var_z=0.0;
double var_t=0.0;
double var_lx=0.0;
double var_ly=0.0;
double var_lz=0.0;
double var_dt=0.0;
double var_dl=0.0;
double var_cell=0.0;
mu::Parser equpar;

// table enumeration starting at 1
double read_table(double tnum, double x) {
  int tabidx=(int)tnum-1;
  if(tabidx<0||tabidx>=(int)myptr->table.size())
    throw "equ: read to undefined table";
  return myptr->table[tabidx]->read(x);
}

double dread_table(double tnum, double x) {
  int tabidx=(int)tnum-1;
  if(tabidx<0||tabidx>=(int)myptr->table.size())
    throw "equ: read to undefined table";
  return myptr->table[tabidx]->dread(x);
}

//----------------------------------------------------------------
// on return str will be cleared if out of time range or no range

#ifdef FDTD_NO_MPI
FDTDSolver::FDTDSolver(int _nx, int _ny, int _nz, const char* _symbol, const char* _map)
#else
FDTDSolver::FDTDSolver(int _nx, int _ny, int _nz, const char* _symbol, const char* _map, MPI_Comm wrld)
#endif
{  
  nx=_nx;
  ny=_ny;
  nz=_nz; // if nz<=2 then 2D simulation
  if(nz<=2) nz=1;
  
  dim=nx*ny*nz;
  
  map=NULL;
  T=A=S=E=NULL;
  cap=con=NULL;
  me=0;
  nproc=1;
  step=0;
  timestep=1.0;
  spacestep=1.0;
  name="";
  neqvar=0;
  step_update_src=step_update_con=step_update_cap=false;
  
  SetSymbolMap(_symbol);
  SetMap(_map);
  
  // equation solver
  equpar.DefineFun("Read", read_table);
  equpar.DefineFun("DRead", dread_table);
  equpar.DefineConst("Pi",M_PI);
  equpar.DefineConst("kb",8.6173324E-5); // eV/K
  equpar.DefineVar("LX",&var_lx);
  equpar.DefineVar("LY",&var_ly);
  equpar.DefineVar("LZ",&var_lz);
  equpar.DefineVar("DT",&var_dt);
  equpar.DefineVar("DL",&var_dl);
  equpar.DefineVar("x",&var_x);
  equpar.DefineVar("y",&var_y);
  equpar.DefineVar("z",&var_z);
  equpar.DefineVar("t",&var_t);
  equpar.DefineVar("cell",&var_cell);

#ifndef FDTD_NO_MPI
  world=wrld;
#endif
                 
  initiate_cells();
  initiate_comm();
  
}

void FDTDSolver::check_expression(string& str, int& t0, int& t1, bool& isequ) {
  
  if(str.find("{")!=str.npos) isequ=true;
  else isequ=false;
  
  istringstream sequ(remove_char(remove_char(replace_char(str,':',' '),'{'),'}') );
  
  sequ>>str;
  if(sequ>>t0) {
    if(!(sequ>>t1)) t1=t0;
  } else {
    t0=-INT_MAX;
    t1=INT_MAX;
  }
  
}

double FDTDSolver::evaluate(int ix, int iy, int iz) {
  double retv;
  try {
    myptr=this;      
    var_x=(double)ix*spacestep;
    var_y=(double)iy*spacestep;
    var_z=(double)iz*spacestep;
    var_t=(double)step*timestep;
    var_cell=A[index(ix,iy,iz)];
    retv=equpar.Eval();
      
    // @mac: must be carefull using brackets
    // std::cerr<<"equ:"<<equ<<" eval: "<<retv<<"\n";
      
  } catch(mu::Parser::exception_type &e){
      std::cout << e.GetMsg() << std::endl;
      throw "equation parsing error";
  }
  
  return retv;
}

double FDTDSolver::evaluate(string& str) {
  double retv=0.0;
  try {
    // non sense to use the following variable here
    myptr=this;
    var_x=0.0;
    var_y=0.0;
    var_z=0.0;
    var_t=(double)step*timestep;
    var_cell=0.0;
    equpar.SetExpr(str);
    retv=equpar.Eval();
  } catch(mu::Parser::exception_type &e){
    std::cout << e.GetMsg() << std::endl;
    throw "equation parsing error";
  }
  return retv;
}

void FDTDSolver::AddEquVar(string vname, double value) {
  if(neqvar>=(FD_MAX_VAR-1)) throw "qllowed number of variable exceeded";
  memset(eqvar[neqvar].name,0,32);
  vname.copy(eqvar[neqvar].name,31);
  eqvar[neqvar].value=value;
  equpar.DefineVar(eqvar[neqvar].name,&(eqvar[neqvar].value));
  neqvar++;
}

void FDTDSolver::SetEquVar(string expr) { // expr -> list of var value pairs
  istringstream ss(replace_char(expr,'=',' '));
  string vn,vv;
  while(ss.good()) {
    if(!(ss>>vn>>vv)) break;
    
    bool found=false;
    
    for(int i=0;i<neqvar;i++) {
      if(vn==eqvar[i].name) {
        eqvar[i].value=evaluate(vv);
        found=true;
        break;
      }
    }
    
    if(not found) AddEquVar(vn,evaluate(vv));
  
  }
}

void FDTDSolver::initiate_cells() {
  T=new double[dim];
  A=new double[dim];
  S=new double[dim];
  E=new double[dim];
  con=new double[dim];
  cap=new double[dim];
  
  for(int i=0;i<dim;i++) {
    T[i]=0.0;
    A[i]=0.0;
    S[i]=0.0;
    E[i]=0.0;
    con[i]=1.0;
    cap[i]=1.0;
  }
  
}

bool FDTDSolver::mycell(int ix, int iy, int iz) {
  int ip=index(ix,iy,iz)/(dim/nproc);
  if(ip==nproc) ip--;
  if(me==ip) return true;
  return false;
}

void FDTDSolver::select_ncell(int cond, int idx, char** nmap) {
  if(cond) {
    int ip=idx/(dim/nproc);
    if(ip==nproc) ip--; // the tail atoms stay in last proc
    if(ip!=me) nmap[ip][idx]=1;
  }
}

void FDTDSolver::send_receive_list() {
#ifndef FDTD_NO_MPI
  MPI_Request req[2*nproc];
  MPI_Status  stat[2*nproc];
  int rq=0;
              
  for(int p=0;p<nproc;p++) { // send to p
    if(p==me) continue;
    if(recvnum[p]) {
      MPI_Send_init(glist[p], recvnum[p], MPI_INT, p, 1, world, &req[rq++]);
    }
  }
  
  // fetch index list to send
  for(int p=0;p<nproc;p++) { // receive from p
    if(p==me) continue;
    if(sendnum[p]) {
      MPI_Recv_init(nlist[p], sendnum[p], MPI_INT, p, 1, world, &req[rq++]);
    }
  }
  
  MPI_Startall(rq,req);
  MPI_Waitall(rq,req,stat);
#endif
}

void FDTDSolver::get_neigcell_index(int idx, int& l, int& r, int& t, int& b, int& n, int& f) {
  int i,j,k;

  get_cellcoord(idx,i,j,k);
  CellBoundary *cc=&(bound[map[idx]]);
  l=i-1; r=i+1; t=j-1; b=j+1; n=k-1; f=k+1;

  
  if(l<0) {if(cc->l==FDCELL_PERIODIC) l=nx-1; else cc->l=FDCELL_EMPTY;}
  if(r>=nx) {if(cc->r==FDCELL_PERIODIC) r=0; else cc->r=FDCELL_EMPTY;}
  if(t<0) {if(cc->t==FDCELL_PERIODIC) t=ny-1; else cc->t=FDCELL_EMPTY;}
  if(b>=ny) {if(cc->b==FDCELL_PERIODIC) b=0; else cc->b=FDCELL_EMPTY;}
  
  if(nz!=1) {
    if(n<0) {if(cc->n==FDCELL_PERIODIC) n=nz-1; else cc->n=FDCELL_EMPTY;}
    if(f>=nz){if(cc->f==FDCELL_PERIODIC) f=0; else cc->f=FDCELL_EMPTY;} 
  } else n=f=0;

  l=index(l,j,k);
  r=index(r,j,k);
  t=index(i,t,k);
  b=index(i,b,k);

  if(nz!=1) {
    n=index(i,j,n);
    f=index(i,j,f);
  }
}

void FDTDSolver::initiate_comm() {

#ifndef FDTD_NO_MPI
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nproc);
#endif
  
  counts=new int[nproc];
  displ=new int[nproc];
  displ[0]=0;

  if(nproc==1) {
    ilo=0;
    ihi=dim;
    counts[0]=ihi-ilo;
  } 
#ifndef FDTD_NO_MPI
  else {
    
    int nl=dim/nproc;
    ilo=nl*me;
    ihi=nl*(me+1);
    if((dim-ihi)<nl) ihi=dim; // take the rest by last proc
    int nc=ihi-ilo;
 
    MPI_Allgather(&nc,1,MPI_INT,counts,1,MPI_INT,world);
    for(int i=1;i<nproc;i++) displ[i]=counts[i-1]+displ[i-1];
    
    char** nmap;
    nmap=new char*[nproc];
    for(int p=0;p<nproc;p++) {
      nmap[p]=new char[dim];
      for(int i=0;i<dim;i++) nmap[p][i]=0;
    }
    
    recvnum=new int[nproc];
    sendnum=new int[nproc];
    nlist=new int*[nproc];
    glist=new int*[nproc];
    ndata=new double*[nproc];
    gdata=new double*[nproc];
    
    for(int p=0;p<nproc;p++) {
      nlist[p]=NULL;
      glist[p]=NULL;
      ndata[p]=NULL;
      gdata[p]=NULL;
    }

    // check neighbor owner
    for(int idx=ilo;idx<ihi;idx++) {
      if(map[idx]<0) continue;
      int l,r,t,b,n,f;
      CellBoundary *cc=&(bound[map[idx]]);
      get_neigcell_index(idx,l,r,t,b,n,f);

      select_ncell(cc->l,l,nmap);
      select_ncell(cc->r,r,nmap);
      select_ncell(cc->t,t,nmap);
      select_ncell(cc->b,b,nmap);
      
      if(nz!=1) {
        select_ncell(cc->n,n,nmap);
        select_ncell(cc->f,f,nmap);
      }
      
    }
    
    for(int p=0;p<nproc;p++) {
      recvnum[p]=0;
      sendnum[p]=0;
      for(int i=0;i<dim;i++)recvnum[p]+=(int)nmap[p][i];
    }
    
    for(int p=0;p<nproc;p++) {
      if(recvnum[p]==0||p==me) continue;
      glist[p]=new int[recvnum[p]];
      gdata[p]=new double[recvnum[p]];
      int n=0;
      for(int i=0;i<dim;i++)
        if(nmap[p][i]) glist[p][n++]=i;
    }
    
    // free...
    for(int p=0;p<nproc;p++) delete nmap[p];
    delete[] nmap;
    
    // ready and sync
    
    int allnum[nproc*nproc];
    MPI_Allgather(recvnum,nproc, MPI_INT, allnum, nproc, MPI_INT, world);
    
    for(int p=0;p<nproc;p++) {
      if(p==me) continue;
      int* ptr=allnum+(p*nproc);
      sendnum[p]=ptr[me];
    }
    
    // allocate send and index buffer
    for(int p=0;p<nproc;p++) {
      if(!sendnum[p]) continue;
      nlist[p]=new int[sendnum[p]];
      ndata[p]=new double[sendnum[p]];
    }
    
    send_receive_list();
        
  }
  MPI_Barrier(world);
#endif
}                       
  
FDTDSolver::~FDTDSolver() {
  if(map) delete[] map;
  if(T) delete[] T;
  if(A) delete[] A;
  if(S) delete[] S;
  if(E) delete[] E;
  if(cap) delete[] cap;
  if(con) delete[] con;
  delete[] counts;
  delete[] displ;
  for (int i=0;i<(int)table.size();i++) delete table[i];

  if(nproc>1) {
    for (int p=0;p<nproc;p++) {
      if(ndata[p]) delete[] ndata[p];
      if(gdata[p]) delete[] gdata[p];
      if(nlist[p]) delete[] nlist[p];
      if(glist[p]) delete[] glist[p];
    }
    
    delete[] ndata;
    delete[] gdata;
    delete[] nlist;
    delete[] glist;
    delete[] sendnum;
    delete[] recvnum;
  }

}

void FDTDSolver::SetTimestep(double dt){
  timestep=dt;
  var_dt=dt;
}
  
void FDTDSolver::SetSpacestep(double h){
  spacestep=h;
  var_dl=h;
  var_lx=(double)nx*h;
  var_ly=(double)ny*h;
  h2inv=1.0/(h*h);
}

int FDTDSolver::search_symbol(char c) {
  for(int i=0;i<(int)bound.size();i++) {
    if(bound[i].sym==c) return i;
  }
  return -1;
}

char FDTDSolver::translate_symbol(int num) {
  if(num<0) return '.';
  if(num>=(int)bound.size()) return '.';
  return bound[num].sym;
}
 
// example:
//  a ++++++

void FDTDSolver::SetSymbolMap(string str) {
  istringstream ss(str.c_str());
  
  while(ss.good()) {
    string tok, sb;
    CellBoundary cc;

    if(ss>>tok>>sb) {
      if(sb.size()<4)
        throw "invalid map symbol entry";
    
      for(int i=0;i<(int)tok.size();i++) {
        cc.sym=tok[i];
        cc.l=sb[0]=='+'?FDCELL_OCCUPIED:(sb[0]=='p'?FDCELL_PERIODIC:0);
        cc.r=sb[1]=='+'?FDCELL_OCCUPIED:(sb[1]=='p'?FDCELL_PERIODIC:0);
        cc.t=sb[2]=='+'?FDCELL_OCCUPIED:(sb[2]=='p'?FDCELL_PERIODIC:0);
        cc.b=sb[3]=='+'?FDCELL_OCCUPIED:(sb[3]=='p'?FDCELL_PERIODIC:0);
        if(nz>1 && sb.size()==6) { // 3D
          cc.n=sb[4]=='+'?FDCELL_OCCUPIED:(sb[4]=='p'?FDCELL_PERIODIC:0);
          cc.f=sb[5]=='+'?FDCELL_OCCUPIED:(sb[5]=='p'?FDCELL_PERIODIC:0);
        } else {
          cc.n=cc.f=0; // FDCELL_EMPTY
        }
        bound.push_back(cc);
      }
    }
  }
}

void FDTDSolver::SetMap(string str) {
  
  // map must be assigned before...
  map = new int[dim];
  if(str.size()!=dim) 
    throw "map string count inconsistent";
  
  for(int i=0;i<dim;i++)
    map[i]=search_symbol(str[i]);
  
}

void FDTDSolver::AddTable(string str) {
  istringstream ss(str);
  while(ss.good()) {
    string stab,tabname;
    if(!(ss>>stab)) break;
    int at=stab.find('@');
    if(at!=stab.npos){
      tabname=stab.substr(0,at);
      stab.erase(0,at+1);
    }
    TableReader *tab=new TableReader(stab,tabname);
    table.push_back(tab);
  }
}

void FDTDSolver::communicate(double* ptr) {
#ifndef FDTD_NO_MPI
  if(nproc>1) {
    MPI_Request req[2*nproc];
    MPI_Status  stat[2*nproc];
    int rq=0;
    
    for(int p=0;p<nproc;p++) {
      if(!sendnum[p]) continue;
      for(int i=0;i<sendnum[p];i++) ndata[p][i]=ptr[nlist[p][i]];
      MPI_Send_init(ndata[p],sendnum[p],MPI_DOUBLE,p,1,world,&req[rq++]);
    }
    for(int p=0;p<nproc;p++) {
      if(recvnum[p])
        MPI_Recv_init(gdata[p],recvnum[p],MPI_DOUBLE,p,1,world,&req[rq++]);
    }
    MPI_Startall(rq,req);
    MPI_Waitall(rq,req,stat);
    
    for(int p=0;p<nproc;p++) {
      if(recvnum[p]) for(int i=0;i<recvnum[p];i++) ptr[glist[p][i]]=gdata[p][i];
    }
  }    
#endif
}

void FDTDSolver::AssignArray(double* ptr, int ix, int iy, int iz, double val) {
  int idx=index(ix,iy,iz);
  if(idx>=dim) throw ("AssignArray out of range: "+
                      as_string(idx)+
                      " ("+as_string(ix)+","+as_string(iy)+","+as_string(iz)+")"
                      ).c_str();
  ptr[idx]=val;
}

// some token in str is removed if already out of range
void FDTDSolver::AssignArray(double* ptr, string& str) {
  if(str.empty()) return;
  
  string tok,val;
  double dv;
  int t0,t1;
  bool isequ;
  istringstream ss(str.c_str());
  str.clear(); // keep rests of expression: accumulated
    
  while(ss.good()) {
    if(ss>>tok>>val) {
      check_expression(val,t0,t1,isequ);
      if(step>=t0 && step<=t1) { // only if in time range
        // if only one time: skip
        if(t0!=t1) {
          if(isequ) 
            str.append(tok+" {"+val+":"+as_string(t0)+":"+as_string(t1)+"} ");
          else
            str.append(tok+" "+val+":"+as_string(t0)+":"+as_string(t1)+" ");
        }
        

        if(tok[0]=='@') {
          int p=val.find('=');
          istringstream scor(replace_char(val.substr(0,p),',',' '));
          val.erase(0,p+1);
          if(isequ) equpar.SetExpr(val);
          int ix,iy,iz,jx,jy,jz;        
          if(!(scor>>ix>>iy>>iz)) throw "wrong coordinate";
          if(scor>>jx>>jy>>jz) { // rectangle block
            if(nz==1){iz=0;jz=0;}
            for(int i=ix;i<=jx;i++)
              for(int j=iy;j<=jy;j++)
                for(int k=iz;k<=jz;k++) 
                  if(mycell(i,j,k)) AssignArray(ptr,i,j,k,(isequ?evaluate(i,j,k):as_double(val))); 
          } else if(mycell(ix,iy,iz)) 
            AssignArray(ptr,ix,iy,iz,(isequ?evaluate(ix,iy,iz):as_double(val)));          
        } else { // fill all
          if(isequ) equpar.SetExpr(val);
          for(int t=0;t<(int)tok.size();t++) {
            for(int i=0;i<nx;i++)
              for(int j=0;j<ny;j++)
                for(int k=0;k<nz;k++)
                  if(tok[t]=='*' ||translate_symbol(map[index(i,j,k)])==tok[t])
                    if(mycell(i,j,k)) AssignArray(ptr,i,j,k,(isequ?evaluate(i,j,k):as_double(val)));
          }
        }
      }
    }
  }
  communicate(ptr);
}

void FDTDSolver::sync_data() {

  AssignSource();
  AssignConductivity();
  AssignCapacity();
  
}

void FDTDSolver::Solve() {
  memcpy(T,A,sizeof(double)*dim); 
  for(int idx=ilo;idx<ihi;idx++) {
    if(map[idx]<0) continue;
    int i,j,k;
    get_cellcoord(idx,i,j,k);
    CellBoundary *cc=&(bound[map[idx]]);
    int l,r,t,b,n,f;
    get_neigcell_index(idx,l,r,t,b,n,f);
    A[idx]=0.0;
    
    if(cc->r) {
      double K=con[r]+con[idx];
      if(K!=0.0) {
        K=con[r]*con[idx]/K;
        A[idx] +=K*(T[r]-T[idx]);
      }
    } 
    
    if(cc->l) {
      double K=con[l]+con[idx];
      if(K!=0.0) {
        K=con[l]*con[idx]/K;
        A[idx] +=K*(T[l]-T[idx]);
      }
    }
    
    if(cc->b) {
      double K=con[b]+con[idx];
      if(K!=0.0) {
        K=con[b]*con[idx]/K;
        A[idx] +=K*(T[b]-T[idx]);
      }
    }
    
    if(cc->t) {
      double K=con[t]+con[idx];
      if(K!=0.0) {
        K=con[t]*con[idx]/K;
        A[idx] +=K*(T[t]-T[idx]);
      }
    }
    
    if(nz!=1) {
      if(cc->n) {
        double K=con[n]+con[idx];
        if(K!=0.0) {
          K=con[n]*con[idx]/K;
          A[idx] +=K*(T[n]-T[idx]);
        }
      }
      if(cc->f) {
        double K=con[f]+con[idx];
        if(K!=0.0) {
          K=con[f]*con[idx]/K;
          A[idx] +=K*(T[f]-T[idx]);
        }
      }            
    }
  }
  
  // source value is valid only for one step.
  // next step must be recreated
  // factor 2 is from averaging of conductivity
  
  for(int i=ilo;i<ihi;i++) {
    A[i]=T[i]+(2.0*A[i]*h2inv+S[i]+E[i])*timestep/cap[i];
    E[i]=S[i]=0.0;
  }
  
}

int FDTDSolver::Step(int n, bool gather) {

  for(int loop=0;loop<n;loop++) {

    sync_data();
    Solve();
    if(!(n==1&&gather)) communicate(A);
    step++;     
  }

#ifndef FDTD_NO_MPI
  // gather after loop
  if(gather) {
    MPI_Allgatherv(&(A[ilo]),(ihi-ilo), MPI_DOUBLE, T, counts, displ, MPI_DOUBLE, world);
    memcpy(A,T,sizeof(double)*dim);
  }
  
  MPI_Barrier(world);
#endif
  
}  

double FDTDSolver::GetAverage() {
  double avg=0.0,gavg;;
  for(int i=ilo;i<ihi;i++) avg+=A[i];
#ifndef FDTD_NO_MPI
  MPI_Reduce(&avg,&gavg,1,MPI_DOUBLE,MPI_SUM,0,world);
  gavg/=(double)dim;
  MPI_Bcast(&gavg,1,MPI_DOUBLE,0,world);
#else
  gavg/=(double)dim;
#endif
  return gavg;
}

void FDTDSolver::Dump(std::ostream &ofl) {
  if(me==0) {    
    ofl<<"# FDTD Solver map file\n"
       <<"# (c) 2012, Yudi Rosandi\n"
       <<"# rosandi@gmail.com\n\n";
    
    ofl<<"dim "<<nx<<" "<<ny<<" "<<nz<<"\n"
       <<"timestep "<<timestep<<"\n"
       <<"spacestep "<<spacestep<<"\n"
       <<"symbol\n";
    
    for(int i=0;i<(int)bound.size();i++) {
      
      ofl<<bound[i].sym<<" ";         
      char l=bound[i].l==FDCELL_PERIODIC?'p':(bound[i].l==FDCELL_OCCUPIED?'+':'-');
      char r=bound[i].r==FDCELL_PERIODIC?'p':(bound[i].r==FDCELL_OCCUPIED?'+':'-');
      char t=bound[i].t==FDCELL_PERIODIC?'p':(bound[i].t==FDCELL_OCCUPIED?'+':'-');
      char b=bound[i].b==FDCELL_PERIODIC?'p':(bound[i].b==FDCELL_OCCUPIED?'+':'-');
      if(nz!=1) {      
        char n=bound[i].n==FDCELL_PERIODIC?'p':(bound[i].n==FDCELL_OCCUPIED?'+':'-');
        char f=bound[i].f==FDCELL_PERIODIC?'p':(bound[i].f==FDCELL_OCCUPIED?'+':'-');
        ofl<<l<<r<<t<<b<<n<<f<<"\n";
      } else  ofl<<t<<l<<r<<b<<"\n";
    }
    
    ofl<<"\n";
    
    int m=0;
    
    ofl<<"map "<<std::endl;
    for(int k=0;k<nz;k++) {
      for(int j=0;j<ny;j++) {
        for(int i=0;i<nx;i++) 
          ofl<<translate_symbol(map[m++]);
        ofl<<std::endl;
      }
    }
    
    ofl<<"\n";
    ofl.flush();
  }

#ifndef FDTD_NO_MPI
  MPI_Barrier(world);
#endif

}

