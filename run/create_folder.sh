#!/bin/sh -l

for i in {0..199}
do
	recon_directory="/gpfs/alpine/gen006/scratch/xf9/dicom/dcm$(printf %03d $i)/"
	mkdir ${recon_directory}

done


echo
echo " ============================ END DEMO ==========================="

