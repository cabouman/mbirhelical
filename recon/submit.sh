#!/bin/sh -l
# FILENAME: submit.sh


#SBATCH -q regular
#SBATCH -N 1
#SBATCH -C knl,quad,cache
#SBATCH -t 08:00:00
#SBATCH -J helical_plug_and_play
#SBATCH -L SCRATCH
#SBATCH -A m3128

export OMP_NUM_THREADS=68
export OMP_PROC_BIND=true   # new recommendations for hybrid MPI/OpenMP
export OMP_PLACES=threads
export NUM_PROCESSES=1


srun -n 1 -c 272 ../src/ct ../data/headscan/param/geom_recon.txt ../data/headscan/param/info_recon.txt ../data/headscan/param/prior_qggmrf.txt recon 15 






echo
echo " ============================ END DEMO ==========================="

