#!/bin/sh -l

for i in {0..199}
do
	recon_directory="/global/cscratch1/sd/wang1698/AAPM_2022/recon/dcm$(printf %03d $i)/"
	mkdir ${recon_directory}

done


echo
echo " ============================ END DEMO ==========================="

