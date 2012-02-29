/*
 *  cmd_create.h
 *  omd
 *
 *  Created by Yudi Rosandi on 2/28/12.
 *  Copyright 2012 TU-Kaiserslautern. All rights reserved.
 *
 */

AtomContainer* create_fcc(istringstream& ss) {
  string ori,mat;
  int xml,yml,zml;
  mdassert(ss>>ori>>xml>>yml>>zml>>mat, "invalid create options");
  AtomContainer* A=new FCC(ori,xml,yml,zml,mat);
  sim->AddAtom(A);
  return A;
}

AtomContainer* create_bcc(istringstream& ss) {
  string ori,mat;
  int xml,yml,zml;
  mdassert(ss>>ori>>xml>>yml>>zml>>mat, "invalid create options");
  AtomContainer* A=new BCC(ori,xml,yml,zml,mat);
  sim->AddAtom(A);
  return A;
}

AtomContainer* create_diamond(istringstream& ss) {
  die("not implemented yet");
  return NULL;
}

AtomContainer* create_atom(istringstream& ss) {
  double x,y,z,vx,vy,vz;
  string mat;
  mdassert(ss>>x>>y>>z>>vx>>vy>>vz>>mat, "invalid create options");
  AtomContainer* A=new FreeAtom(x,y,z,vx,vy,vz,mat);
  sim->AddAtom(A);
  return A;
}

AtomContainer* create_projectile(istringstream& ss) {
  string mat;
  double x,y,z,en,th,ph;
  mdassert(ss>>x>>y>>z>>th>>ph>>en>>mat, "invalid create options");
  AtomContainer* A=new Projectile(mat,x,y,z,th,ph,en);
  sim->AddAtom(A);
  return A;
}

AtomContainer* create_import(istringstream& ss) {
  string fname;
  int aid;
  
  mdassert(ss>>fname, "file name required");
  AtomContainer* A=new AtomContainer();
  if(ss>>aid) A->Import(fname,aid);
  else A->Import(fname);
  sim->AddAtom(A);
  return A;
}

AtomContainer* create_combine(istringstream& ss) {
  AtomContainer* A=new AtomContainer;
  while(ss.good()) {
    string c;
    if(ss>>c) {
      AtomContainer* C=sim->SearchContainer(c);
      A->Combine(C);
      sim->DeleteAtom(c);
    }
  }
  return A;
}

void cmd_create(istringstream& ss) {
  string name;
  string str;
  mdassert(ss>>name>>str, "invalid create options");
  
  if(str=="fcc") create_fcc(ss)->SetName(name);
  else if(str=="bcc") create_bcc(ss)->SetName(name);
  else if(str=="diamond") create_diamond(ss)->SetName(name);
  else if(str=="import") create_import(ss)->SetName(name);
  else if(str=="atom") create_atom(ss)->SetName(name);
  else if(str=="projectile") create_projectile(ss)->SetName(name);
  else if(str=="combine") create_combine(ss)->SetName(name);
  else die("unimplemented create command");
  
}
