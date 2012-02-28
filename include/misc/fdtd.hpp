
#ifndef __FDTD_H__
#define __FDTD_H__

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cstdlib>
#include "omd/treader.h"

#ifndef FDTD_NO_MPI
#include <mpi.h>
#endif


// For 2D problems...

using std::vector;
using std::ifstream;
using std::istringstream;
using namespace omd;

#define FDCELL_EMPTY 0
#define FDCELL_OCCUPIED 1
#define FDCELL_PERIODIC 2

#define FD_MAX_VAR 50

struct CellBoundary {
  char sym;
  // left right top bottom near far
  int  l,r,t,b,n,f;
};


class FDTDSolver {
  friend double read_table(double, double);
  friend double dread_table(double, double);

  double*  A;
  double*  T;
  double*  S;
  double*  E;
  double*  con;
  double*  cap;
  string str_src;
  string str_cap;
  string str_con;
  
  int neqvar;
  struct {char name[32]; double value;} eqvar[FD_MAX_VAR];
  
  int me; // mpi rank
  int nproc; // number of proc
  int* counts;
  int* displ;
  int* sendnum;
  int* recvnum;
  int** nlist; // index to send
  int** glist; // index to receive
  double** ndata; // data to send
  double** gdata; // data to receive
  
  int dim,nx,ny,nz;
  int ilo,ihi; // my block
  int step;
  double timestep;
  double spacestep;
  double h2inv; // = 1/spacestep^2
  bool step_update_src,step_update_con,step_update_cap;
  
  // map contains index to the coeficients
  int* map;
  vector<CellBoundary> bound;
  vector<TableReader*> table;

  int search_symbol(char c);
  char translate_symbol(int num);
  
  void get_cellcoord(int idx, int& i, int& j, int& k){
    int xy=nx*ny;
    int a=idx%xy;
    k=idx/xy;
    j=a/nx;
    i=a%nx;
  }
  
  int index(int i, int j, int k) {return ((k*nx*ny)+(i+j*nx));}
  void initiate_cells();
  void initiate_comm();
  void sync_data();
  void get_neigcell_index(int,int&,int&,int&,int&,int&,int&);
  void send_receive_list();
  void select_ncell(int,int,char**);
  bool mycell(int,int,int);
  double evaluate(string&);
  double evaluate(int,int,int);
  void check_expression(string&,int&,int&,bool&);
  void communicate(double* ptr);
  
public:
  string name;
#ifdef FDTD_NO_MPI
  FDTDSolver(int, int, int, const char*, const char*);
#else
  MPI_Comm world;
  FDTDSolver(int, int, int, const char*, const char*, MPI_Comm);  
#endif
  
  virtual ~FDTDSolver();
  
  void AssignArray(double*,int,int,int,double);
  void AssignArray(double*,string&);
  double* GetSourceVector() {return S;}
  double* GetExtVector() {return E;}
  double* GetVector(){return A;}
  double* GetConductivity(){return con;}
  double* GetCapacity(){return cap;}
  double GetTimestep(){return timestep;}
  double GetSpacestep(){return spacestep;}
  int  GetStep(){return step;}
  double GetAverage();
  int GetDim(){return dim;}
  void GetDim(int& _nx, int& _ny, int& _nz){_nx=nx;_ny=ny;_nz=nz;}
  void AddTable(string);  
  void SetSymbolMap(string str);
  void AddEquVar(string, double);
  void SetEquVar(string);
  
  // these methods may be altered runtime
  void SetMap(string str);
  
  void AssignCell(){if(str_src.size()) AssignArray(A,str_src);}
  void SetCell(string str){AssignArray(A,str);} // only for initialization

  void AssignSource(){if(str_src.size()) AssignArray(S,str_src);}
  void SetSource(string str){str_src=str; AssignArray(S,str);}
  
  void AssignConductivity(){if(str_con.size()) AssignArray(con,str_con);}
  void SetConductivity(string str){str_con=str; AssignArray(con,str);}  
  
  void AssignCapacity(){if(str_cap.size()) AssignArray(cap,str_cap);}
  void SetCapacity(string str){str_cap=str; AssignArray(cap,str);}
  
  void SetTimestep(double dt);
  void SetSpacestep(double h);
  virtual void Solve();
  int Step(int n=1, bool gather=true);
    
  void Dump(std::ostream& ofl);
};

#endif
