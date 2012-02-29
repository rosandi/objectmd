/*
 *  cmd_velocity.h
 *  omd
 *
 *  Created by Yudi Rosandi on 2/28/12.
 *  Copyright 2012 TU-Kaiserslautern. All rights reserved.
 *
 */

void cmd_velocity(istringstream& ss) {
  string cname;
  double vx,vy,vz;
  sim->SearchContainer(cname)->SetVelocity(vx,vy,vz);
}
