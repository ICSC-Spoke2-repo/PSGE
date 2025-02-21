#!/bin/bash

#SBATCH --time=240:00:00
##SBATCH --nodes=1
##SBATCH --sockets-per-node=1
##SBATCH -n 1
##SBATCH --cores-per-socket=2
#SBATCH -p moon
#SBATCH --job-name=17
#SBATCH --err=run.err
#SBATCH --out=run.out

date

singularity exec --bind=/blasco /home/ciabattoni/containers/g4_11_1_HPC.sif mpiexec -n 1 /home/ciabattoni/ICSC_G4_HPC/BoGEMMS-HPC/BoGEMMS-HPC-build/bogemms COSI.conf 0 beam.mac
