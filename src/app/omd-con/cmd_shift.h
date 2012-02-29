/*
 *  cmd_shift.h
 *  omd
 *
 *  Created by Yudi Rosandi on 2/28/12.
 *  Copyright 2012 TU-Kaiserslautern. All rights reserved.
 *
 */

void cmd_shift(istringstream& ss) {
  string cname;
  double sx,sy,sz;
  mdassert(ss>>cname>>sx>>sy>>sz,"invalid shift command");
  AtomContainer *A=sim->SearchContainer(cname);
  A->Shift(sx,sy,sz);
}
