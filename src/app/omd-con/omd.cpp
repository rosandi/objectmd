// console program for ObjectMD
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cstdlib>
// #include <muParser/muParser.h>
#include <mpi.h>
#include <omd/systemgrid.h>
#include <omd/omdtool.h>
#include <omd/param.h>
#include <omd/piper.h>

using namespace omd;

class simulation {
  MDSystemGrid* sim;
  
  // commands.....
#include "cmd_create.h"
#include "cmd_shift.h"
#include "cmd_temperature.h"
#include "cmd_velocity.h"
#include "cmd_integrator.h"
#include "cmd_interaction.h"
#include "cmd_group.h"
#include "cmd_detect.h"
#include "cmd_control.h"
  
public:
  
  bool parse(string cmd) {
    bool retval=true;
    int sep;
    vector<string> cmds;
    int pcmt=cmdpar.find('#');
    if(pcmt!=cmdpar.npos) cmdpar.erase(pcmt,cmdpar.npos);
    
    do {
      sep=cmdpar.find(";");
      cmds.push_back(cmdpar.substr(0,sep));
      cmdpar.erase(0,sep+1);
    } while(sep!=cmdpar.npos);
    
    for(int n=0;n<(int)cmds.size();n++) {
      if(looplevel) cmds[n]=string("loop ")+cmds[n];
      istringstream ss(cmds[n]);
      string cmd;
      if(!(ss>>cmd)) continue;
      
      try {
        if(cmd=="create") cmd_create(ss);
        else if(cmd=="shift") cmd_shift(ss);
        else if(cmd=="temperature") cmd_temperature(ss);
        else if(cmd=="velocity") cmd_velocity(ss);
        else if(cmd=="interaction") cmd_interaction(ss);
        else if(cmd=="group") cmd_group(ss);
        else if(cmd=="detect") cmd_detect(ss);
        else if(cmd=="control") cmd_control(ss);
        else {
          if(me==0) sayer("undefined command: "+cmd);
          break;
        }
      } catch(const char* sterr) {
        sayer(sterr);
        loading=false;
        loopbreak=true;
        break;
      } catch(...) {
        sayer("undefined error");
        loading=false;
        loopbreak=true;
        break;
      }
      
    }
    
    return retval;
  }
  
  void inputline() {
    bool run=true;
    while(run) {
      char sline[2048];
      if(me==0) {
        std::cin.getline(sline,2047);
        if(std::cin.eof()) sprintf(sline,"exit");
      }
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(sline,2048,MPI_CHAR,0,MPI_COMM_WORLD);
      run=parse(sline);
    }
  }
  
  
};


int main(int argc, char* argv[]) {
  
  

  return 0;
}
