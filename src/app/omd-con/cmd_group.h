/*
 *  cmd_group.h
 *  omd
 *
 *  Created by Yudi Rosandi on 2/29/12.
 *  Copyright 2012 TU-Kaiserslautern. All rights reserved.
 *
 */

void group_box(istringstream& ss, string& name) {
  double x0,y0,z0,x1,y1,z1;
  mdassert(ss>>x0>>y0>>z0>>x1>>y1>>z1, "invalid box coordinate");
  sim->AddAtomGroup(name)->SelectBox(x0,x1,y0,y1,z0,z1);
}

void group_ibox(istringstream& ss, string& name) {
  double x0,y0,z0,x1,y1,z1;
  mdassert(ss>>x0>>y0>>z0>>x1>>y1>>z1, "invalid box coordinate");
  sim->AddAtomGroup(name)->SelectInverseBox(x0,x1,y0,y1,z0,z1);
}

void group_type(istringstream& ss, string& name) {
  int ty;
  mdassert(ss>>ty, "invalid atom type");
  sim->AddAtomGroup(name)->SelectType(ty);
}

void group_index(istringstream& ss, string& name) {
  vector<int> vidx;
  while(ss.good()) {
    int ii;
    if(ss>>ii) vidx.push_back(ii);
  }
  int* ii=new int[vidx.size()];
  for(int i=0;i<(int)vidx.size();i++) ii[i]==vidx[i];
  sim->AddAtomGroup(name)->Select(ii);
  delete[] ii;
}

void group_g(istringstream& ss, string& name) {
  double x,y,z;
  mdassert(ss>>x>>y>>z, "invalid coordinate");
  sim->AddAtomGroup(name)->SelectGT(x,y,z);
}

void group_ge(istringstream& ss, string& name) {
  double x,y,z;
  mdassert(ss>>x>>y>>z, "invalid coordinate");
  sim->AddAtomGroup(name)->SelectGE(x,y,z);
}

void group_l(istringstream& ss, string& name) {
  double x,y,z;
  mdassert(ss>>x>>y>>z, "invalid coordinate");
  sim->AddAtomGroup(name)->SelectLT(x,y,z);
}

void group_le(istringstream& ss, string& name) {
  double x,y,z;
  mdassert(ss>>x>>y>>z, "invalid coordinate");
  sim->AddAtomGroup(name)->SelectLE(x,y,z);
}

void group_union(istringstream& ss, string& name) {
  string c;
  mdassert(ss>>c, "container name required");
  sim->AddAtomGroup(name)->Union(*(sim->SearchContainer(c)));
}

void cmd_group(istringstream& ss) {
  string name,alg;
  mdassert(ss>>name>>alg,"invalid group options");
  
  if(alg=="box") group_box(ss,name);
  else if(alg=="ibox") group_ibox(ss,name);
  else if(alg=="type") group_type(ss,name);
  else if(alg=="index") group_index(ss,name);
  else if(alg=="greater") group_g(ss,name);
  else if(alg=="greater-equal") group_ge(ss,name);
  else if(alg=="less") group_l(ss,name);
  else if(alg=="less-equal") group_le(ss,name);
  else if(alg=="union") group_union(ss,name);
  else die("unimplementeed grouping algorithm: "+alg);
  
}
