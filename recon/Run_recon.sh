#!/bin/csh

set d=param
set fprior=prior_qggmrf.txt

set f=recon
../src/ct ../data/headscan/$d/geom_recon.txt ../data/headscan/$d/info_recon.txt ../data/headscan/$d/$fprior {$f} 15 


