# BoGEMMS-HPC test case

## Description
This is a test case for using BoGEMMS-HPC to perform a simulation and analyze the output.
The BoGEMMS-HPC repository is available [here](https://www.ict.inaf.it/gitlab/icsc_g4_hpc/BoGEMMS-HPC).
 
The test case involves the irradiation with 122 keV photons of a cesium iodide (CsI) module, contained in an aluminum casing and coupled with a photomultiplier tube (PMT). The X photons hit perpendicularly the CsI in 40 different positions. When a 122 keV photon hit the CsI, optical photons are generated inside the scintillator and part of them is collected and absorbed by the PMT. The final output is the relative response of the PMT for each position, which is proportional to the average number of optical photons detected, compared to the corresponding experimental data obtained in laboratory.

## Usage
There are two directories: **'Simulations'** and **'Analysis'**. 'Simulations' contains the files for performing the simulation, while 'Analysis' provides the files for analyzing the simulation output and plotting the results.

### Simulation
Files provided in 'Simulations':

- **beam.mac**: macro file for the 122 keV photons
- **runCLAIRE.conf**: configuration file for the simulation
- **runCLAIRE.py**: script to start the simulation
- **runCLAIRE.txt**: input file for runCLAIRE.py
- **bogemms_container.slurm**: file to start the simulation with slurm
- **bogemms_container.sh**: file to start the simulation with nohup
- **geom/**: directory with the geometry files
- **phys/**: directory with the physics files
- **cad_files/**: directory with the CAD files (for the casing and Teflon)
- **config/**: directory containing the coordinates of the beam positions (stored in **claire_beams.dat**)

Update 'runCLAIRE.txt' to define the parameters of the simulation:

- **run_start** and **run_stop**: define the initial and final run (the total number of runs will be run_stop - run_start + 1)
- **csi_absl**: define the CsI absorption length in mm
- **refl_type**: define the reflective surface type
- **G4PATH**: define the path where to create the simulation directory
- **N_in**: define the number of events (122 keV photons) to simulate
- **num_threads** and **num_tasks**: define the number of threads and MPI tasks
- **cad_path**: define the path to the CAD file
- **job_name**: define the job name (if executed with slurm)
- **isslurm**: define how to run the simulation (0: nohup, 1: SLURM)

If you use the SLURM scheduler, update 'bogemms_container.slurm' to match the the number of threads and/or tasks used, the correct partition and the correct paths for the container and executable. 

When using MPI, the total number of events (N_in) is evenly distributed among the MPI tasks. However, N_in should be a multiple of the number of MPI tasks. If not, each task will receive N_in / num_tasks, rounded down to the nearest integer, which results in fewer total simulated events than intended.

When using multithreading, Geant4 automatically handles the distribution of events among the threads. In this case, N_in is always preserved, even if it is not a multiple of the number of threads. However, the workload distribution among the threads may be uneven.

When using both MPI and multithreading, Geant4 first distributes the events evenly among the MPI tasks. Within each MPI task, the events are then divided among the threads (whose number is specified by 'num_threads').

To start the simulation run:

```
$ python3 runCLAIRE.py runCLAIRE.txt
```

At the end of the simulation, the results are stored in the following directory:

{G4PATH}/{bgo_absl_type}/REFL_TYPE{refl_type}/{N_in}/BEAM{beam_id}/run{run_id}/

where 'beam_id' goes from 1 to 40 and 'run_id' from 'run_start' to 'run_stop'.

### Analysis
In 'Analysis' you can find the following files:

- **filterCLAIREPMT.py**: python script to filter the simulation output
- **plotCLAIREPMT.py**: python script to plot the result
- **config/**: directory containing the file (**CLAIRE_pulse.dat**) with the experimental relative response for each beam position

Use filterCLAIRE.py to filter the output files from the simulation and write the relevant information to 'PMT.dat'.
In 'Analysis/' run:

```
python3 filterCLAIREPMT.py <SIM_DIR> <run_start> <run_stop>
```

Example:

```
$ python3 filterCLAIREPMT.py ../Simulations/2000/REFL_TYPE14/1000/ 1 1
```

This will create a directory like:

{bgo_absl_type}/REFL_TYPE{refl_type}/{N_in}/

In this directory a file called 'PMT.dat' is created with the following columns:

- **EventID**: ID of the event
- **N_BEAM**: ID of the beam
- **Edep[keV]**: energy deposit [in keV] of the detected optical photon

Finally, to plot the relative response for each beam position, compared with the experimental one, run:

```
$ python3 plotCLAIREPMT.py <file1> <file2> ...
```

Example:

```
$ python3 plotCLAIREPMT.py 2000/REFL_TYPE14/1000/PMT.dat
```

## Acknowledgment
This is supported by Centro Nazionale di Ricerca in High-Performance Computing, Big Data and Quantum Computing (CN_00000013 - CUP C53C22000350006).

## License
Copyright 2025 Valentina Fioretti, Alex Ciabattoni

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
