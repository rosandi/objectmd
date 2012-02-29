/*
 *  cmd_temperature.h
 *  omd
 *
 *  Created by Yudi Rosandi on 2/28/12.
 *  Copyright 2012 TU-Kaiserslautern. All rights reserved.
 *
 */

void cmd_temperature(istringstream& ss) {
  string cname;
  double t;
  mdassert(ss>>cname>>t, "temperature value required");
  sim->SearchContainer(cname)->SetTemperature(t);
}
