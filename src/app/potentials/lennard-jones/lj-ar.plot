#!/usr/bin/gnuplot

print "# Argon Lennard-Jones potential"
print "#$ Format SETFL"

pi=3.141592653589793
eps=0.010322
sig=3.404
rcut=2.5*sig
td=1.0
trad=rcut-2.*td
tf=pi/(2.*td)
tp=pi*(td-rcut)/(2.*td)
lj(x)=x*4.0*eps*((sig/x)**12.-(sig/x)**6.)
tr(x)=(0.5-0.5*sin(tf*x+tp))

set xrange [0:rcut]
set table
set format y "%.15e"
set format x "%.15e"
set sample 1000

p (x<trad)?lj(x):lj(x)*tr(x)
