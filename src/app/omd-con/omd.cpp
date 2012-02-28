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
  MDSystemGrid sim;
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
        else if(cmd=="exit"||cmd=="quit"||cmd=="q") retval=0;
        else if(cmd=="step") cmd_step(ss);
        else if(cmd=="reset") cmd_reset(ss);
        else if(cmd=="use") cmd_use(ss);
        else if(cmd=="del"||cmd=="delete") cmd_del(ss);
        else if(cmd=="set") cmd_set(ss);
        else if(cmd=="fill") cmd_fill(ss);
        else if(cmd=="buffer"||cmd=="buf") cmd_buffer(ss);
        else if(cmd=="loop") cmd_loop(ss);
        else if(cmd=="load") cmd_load(ss);
        else if(cmd=="quiet") quiet=true;
        else if(cmd=="verbose") quiet=false;
        else if(cmd=="table") cmd_table(ss);
        else if(cmd=="show") {ROOT(cmd_show(ss));}
        else if(cmd=="plot") {ROOT(cmd_plot(ss));}
        else if(cmd=="sleep") {ROOT(cmd_sleep(ss));}
        else if(cmd=="dump") {ROOT(cmd_dump(ss));}
        else if(cmd=="echo") {ROOT(cmd_echo(ss));}
        else if(cmd=="ping") {ROOT(cmd_ping(ss));}
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
        std::cerr<<getpid()<<"unexpected too\n";
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
