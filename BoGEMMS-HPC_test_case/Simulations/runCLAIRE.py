"""
 runCLAIRE.py  -  description
 ---------------------------------------------------------------------------------
 running the Geant4 simulation for the CLAIRE experiment
 ---------------------------------------------------------------------------------
 copyright            : (C) 2022 Valentina Fioretti
 email                : valentina.fioretti@inaf.it
 ----------------------------------------------
// ********************************************************************
//
// Copyright 2025 Valentina Fioretti
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// **********************************************************************
"""

import sys, os
import subprocess
import shutil
import random


# Import line command parameters
arg_list = sys.argv
file_input = arg_list[1]

# Read input parameters from file
input_params = []
with open(file_input, "r") as f_in:
        for line in f_in:
                line = line.strip()
                if not line.startswith("#") and line != "":
                        columns = line.split("=")
                        input_params.append(columns[-1].strip())

run_start = int(input_params[0])
run_stop = int(input_params[1])
ABSL = input_params[2] # mm
REFL_type = int(input_params[3])
G4PATH = input_params[4]
N_in = int(input_params[5])
num_threads = int(input_params[6])
num_tasks = int(input_params[7])
cad_path = input_params[8]
job_name = input_params[9]
isslurm = int(input_params[10]) # 0: nohup, 1: slurm

REFL_dir = "/REFL_TYPE"+str(REFL_type)

# reading beam position

vecNBeam = []
vecXBeam = []
vecYBeam = []

filepath = "./config/claire_beams.dat"
f_beam = open(filepath, 'r')
for line in f_beam:
	line = line.strip()
	if not line.startswith("#"):
		columns = line.split()
		columns[0] = int(columns[0])
		columns[1] = float(columns[1])
		columns[2] = float(columns[2])

		vecNBeam.append(columns[0])
		vecXBeam.append(columns[1])
		vecYBeam.append(columns[2])
	
f_beam.close()


N_runs = run_stop - run_start + 1

path_sim_main = G4PATH+ABSL+REFL_dir+"/"+str(N_in)


if not os.path.exists(path_sim_main):
	print("... Creating "+path_sim_main)
	os.makedirs(path_sim_main)
	
# loop in the beams
for jbeam in range(len(vecNBeam)):

	print("Simulating beam n. ..."+str(jbeam+1))
	
	beam_dir = "/BEAM"+str(vecNBeam[jbeam])
	if not os.path.exists(path_sim_main+beam_dir):
		print("... Creating "+path_sim_main)
		os.makedirs(path_sim_main+beam_dir)
	
	os.chdir(path_sim_main+beam_dir)
	
	# loop in the run
	for jrun in range(run_start, run_stop + 1):
		rundir = 'run'+str(jrun)+'/'
		gorun = 0
		# creating directory
		if not os.path.exists(rundir):
			print("... Creating "+rundir)
			os.makedirs(rundir)
			gorun = 1
		else:
			if len(os.listdir('./'+rundir)) == 0:
				gorun = 1
			else:
				check_fits = 0
				check_fits_gz = 0
				for fname in os.listdir('./'+rundir):
					if fname.endswith('.fits'):
						check_fits = 1
					if fname.endswith('.fits.gz'):
						check_fits_gz = 1
				if check_fits == 1: 
					shutil.rmtree(rundir)
					os.makedirs(rundir)    
					gorun = 1
				else:
					if check_fits_gz == 0:
						shutil.rmtree(rundir)
						os.makedirs(rundir)    
						gorun = 1
					else:
						gorun = 0
	if gorun == 1:
		# copying files in directory
		subprocess.call(["pwd"])
		subprocess.call(["cp", G4PATH+'/COSI.conf', rundir])
		subprocess.call(["cp", G4PATH+'/currentEvent.rndm', rundir])
		subprocess.call(["cp", G4PATH+'/beam.mac', rundir])
		subprocess.call(["cp", G4PATH+'/'+file_input, rundir])
		if isslurm: 
			subprocess.call(["cp", G4PATH+'/bogemms_container.slurm', rundir])
		else:
			subprocess.call(["cp", G4PATH+'/bogemms_container.sh', rundir])
		# changing dir to the G4 run
		os.chdir(rundir)

		f = open('COSI.conf', "r")
		list_of_lines = f.readlines()
		l = 0
		for line in list_of_lines:
			if line.startswith('PHYS.COSI.2.CsI.ABSL'):
				list_of_lines[l] = 'PHYS.COSI.2.CsI.ABSL = '+ABSL+'\n' # mm
			if line.startswith('PHYS.COSI.2.OPTSURFACE.WRAPPER'):
				list_of_lines[l] = 'PHYS.COSI.2.OPTSURFACE.WRAPPER = '+str(REFL_type)+'\n' # CsI-Teflon
			if line.startswith('GEOM.CAD.PATH'):
				list_of_lines[l] = 'GEOM.CAD.PATH = '+cad_path+'\n'
			if (num_threads > 1):
				if line.startswith('RUN.MT.ACTIVATE'): 
					list_of_lines[l] = 'RUN.MT.ACTIVATE = 1\n'
				if line.startswith('MT.NUM.THREADS'):
					list_of_lines[l] = 'MT.NUM.THREADS = '+str(num_threads)+'\n'
			l = l + 1
		f = open('COSI.conf', "w")
		f.writelines(list_of_lines)
		f.close()

		f = open('beam.mac', "r")
		list_of_lines = f.readlines()
		N_in_task = int(N_in / num_tasks)
		l = 0
		for line in list_of_lines:
			if line.startswith('/gps/pos/centre'):
				list_of_lines[l] = '/gps/pos/centre '+str(vecXBeam[jbeam])+' '+str(vecYBeam[jbeam])+' -100 mm\n'
			if line.startswith('/run/beamOn'):
				list_of_lines[l] = '/run/beamOn '+str(N_in_task)+'\n'
			l = l + 1
		f = open('beam.mac', "w")
		f.writelines(list_of_lines)
		f.close()

		if isslurm:
			f = open('bogemms_container.slurm', "r")
			list_of_lines = f.readlines()
			l = 0
			for line in list_of_lines:
				if line.startswith('#SBATCH --job-name'):
					list_of_lines[l] = '#SBATCH --job-name='+job_name+'\n'
				l = l + 1
			f = open("bogemms_container.slurm", "w")
			f.writelines(list_of_lines)
			f.close()
		else:
			f = open('bogemms_container.sh', "r")
			list_of_lines = f.readlines()
			l = 0
			for line in list_of_lines:
				if line.startswith('bogemms'):
					if num_tasks > 1:
						list_of_lines[l] = 'mpiexec -n '+str(num_tasks)+' bogemms COSI.conf 0 beam.mac\n'
					else:
						list_of_lines[l] = 'bogemms COSI.conf 0 beam.mac\n'
				l = l + 1
			f = open("bogemms_container.sh", "w")
			f.writelines(list_of_lines)
			f.close()
		
		f = open('currentEvent.rndm', "r")
		list_of_lines = f.readlines()
		list_of_lines[3] = str(random.randrange(sys.maxsize))[:6]+'\n'
		list_of_lines[4] = str(random.randrange(sys.maxsize))[:6]+'\n'
		f = open('currentEvent.rndm', "w")
		f.writelines(list_of_lines)
		f.close()

		print("... Running "+rundir)
		if isslurm == 0: subprocess.Popen(["nohup", "bogemms", "COSI.conf", "0", "beam.mac", "&"])
		#if isslurm == 0: subprocess.Popen(["sh", "bogemms_container.sh", ">", "out.log"])	
		if isslurm == 1: subprocess.Popen(["sbatch", "bogemms_container.slurm"])	



