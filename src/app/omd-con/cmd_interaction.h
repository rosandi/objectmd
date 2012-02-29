/*
 *  cmd_interaction.h
 *  omd
 *
 *  Created by Yudi Rosandi on 2/28/12.
 *  Copyright 2012 TU-Kaiserslautern. All rights reserved.
 *
 */

void inter_eam(istringstream& ss, string& a, string& b) {
  string mat,tfile;
  mat=sim->SearchContainer(a)->GetMaterialFile();
  mdassert(ss>>tfile,"eam table name required");
  mdassert(mat!=sim->SearchContainer(b)->GetMaterialFile(), 
           "incompatible interaction for "+a+" and "+b);
  sim->AddForce(new TForceEAM(mat),a,b);
}

void inter_table(istringstream& ss, string& a, string& b) {
  string tfile;
  mdassert(ss>>tfile,"pair table name required");
  sim->AddForce(new TForcePair(tfile),a,b);
}

void cmd_interaction(istringstream& ss) {
  string a,b,sty;
  mdassert(ss>>a>>b>>sty, "invalid interaction options");
  
  if(sty=="eam") inter_eam(ss,a,b);
  else if(sty=="table") inter_table(ss,a,b);
  else die("unimplemented interaction: "+sty);
  
}
