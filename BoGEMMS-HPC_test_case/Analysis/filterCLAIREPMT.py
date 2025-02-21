"""
 filterCLAIRE.py  -  description
 ---------------------------------------------------------------------------------
 filtering the events of the TEST
 ---------------------------------------------------------------------------------
 copyright            : (C) 2021 Valentina Fioretti
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

from astropy.io import fits
import numpy as np
import sys, os

import glob

# Import line command parameters
arg_list = sys.argv
PATH_SIM = arg_list[1]
run_start = int(arg_list[2])
run_stop = int(arg_list[3])

# Read input parameters from config file
input_params = []
file_input = PATH_SIM + "/BEAM1/run" + str(run_start) + "/runCLAIRE.txt"
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

# Volume IDs
sipad_copyno = 200	
glass_copyno = 201	
bialkali_copyno = 202	
pmt_copyno = 203	
aux_copyno = 204	
csi_copyno = 100	
teflon_copyno = 101	
almain_copyno = 102	
albot_copyno = 103	

write_to_file = False

REFL_dir = "/REFL_TYPE" + str(REFL_type)

N_beams = 40

# absorbed
vecPMTEventID = []
vecPMTEdep = []
vecPMTNbeam = []


vecGammaEent_inc = []
vecGammaEexit_inc = []

# entering
vecPMTEventID_inc = []
vecPMTEent_inc = []
vecPMTNbeam_inc = []


vecGammaEent = []
vecGammaEexit = []

N_runs = run_stop - run_start + 1

path_sim = "../Simulations/"+ABSL+REFL_dir+"/"+str(N_in)	
N_tot = N_in*N_runs
path_output = "./"+ABSL+REFL_dir+"/"+str(N_tot)+"/"	
	


for jbeam in range(N_beams):
	nbeam = jbeam + 1
	beam_dir = "/BEAM"+str(nbeam)
	
	vecPMTEventID_histo = []
	vecPMTEdep_histo = []
	vecPMTNbeam_histo = []

	for jrun in range(run_start, run_stop + 1):
		rundir = "/run"+str(jrun)
		N_tasks = len(glob.glob1(path_sim+beam_dir+rundir,"*0.task*.fits*"))
		N_in_per_task = int(N_in / N_tasks) # number of events per task

		for task in range(N_tasks):
			N_fits = len(glob.glob1(path_sim+beam_dir+rundir,"*.task"+str(task)+".fits.gz"))
	
			for jfits in range(N_fits):

				print('%%%%%%%%%%%%%% READING BoGEMMS-HPC FILE: ', path_sim+beam_dir+rundir+'/xyz.'+str(jfits)+'.task'+str(task)+'.fits.gz')
				hdulist = fits.open(path_sim+beam_dir+rundir+'/xyz.'+str(jfits)+'.task'+str(task)+'.fits.gz')

				tbdata = hdulist[1].data
				evt_id = tbdata.field('EVT_ID') + N_in_per_task * task + N_in * (jrun-1)
				trk_id = tbdata.field('TRK_ID')
				parent_trk_id = tbdata.field('PARENT_TRK_ID')
				vol_id = tbdata.field('VOLUME_ID')
				moth_id = tbdata.field('MOTHER_ID')
				ene_ent = tbdata.field('E_KIN_ENT')
				ene_exit = tbdata.field('E_KIN_EXIT')
				e_dep = tbdata.field('E_DEP')
				mdx_ent = tbdata.field('MDX_ENT')
				mdy_ent = tbdata.field('MDY_ENT')
				mdz_ent = tbdata.field('MDZ_ENT')
				mdx_exit = tbdata.field('MDX_EXIT')
				mdy_exit = tbdata.field('MDY_EXIT')
				mdz_exit = tbdata.field('MDZ_EXIT')
				x_ent = tbdata.field('X_ENT')
				y_ent = tbdata.field('Y_ENT')
				z_ent = tbdata.field('Z_ENT')
				x_exit = tbdata.field('X_EXIT')
				y_exit = tbdata.field('Y_EXIT')
				z_exit = tbdata.field('Z_EXIT')
				part_id = tbdata.field('PARTICLE_ID')
				gtime_ent = tbdata.field('GTIME_ENT')

				index = 0
				where_sameevent = [index]
				temp_index = index
				while (where_sameevent[-1] < (len(evt_id)-1)):
					while evt_id[temp_index] == evt_id[temp_index+1]:
						where_sameevent.append(temp_index+1)
						temp_index += 1
					else:
						sameev_evt_id = evt_id[where_sameevent]
						sameev_trk_id = trk_id[where_sameevent]
						sameev_vol_id = vol_id[where_sameevent]
						sameev_moth_id = moth_id[where_sameevent]
						sameev_ene_ent = ene_ent[where_sameevent]
						sameev_ene_exit = ene_exit[where_sameevent]
						sameev_ene_dep = e_dep[where_sameevent]
						sameev_mdx_ent = mdx_ent[where_sameevent]
						sameev_mdy_ent = mdy_ent[where_sameevent]
						sameev_mdz_ent = mdz_ent[where_sameevent]
						sameev_mdx_exit = mdx_exit[where_sameevent]
						sameev_mdy_exit = mdy_exit[where_sameevent]
						sameev_mdz_exit = mdz_exit[where_sameevent]
						sameev_x_ent = x_ent[where_sameevent]
						sameev_y_ent = y_ent[where_sameevent]
						sameev_z_ent = z_ent[where_sameevent]
						sameev_x_exit = x_exit[where_sameevent]
						sameev_y_exit = y_exit[where_sameevent]
						sameev_z_exit = z_exit[where_sameevent]
						sameev_part_id = part_id[where_sameevent]
						sameev_gtime_ent = gtime_ent[where_sameevent]


						where_PMT = np.where((sameev_vol_id == bialkali_copyno) & (sameev_part_id == -22) & (sameev_ene_dep > 0))
						where_Gamma = np.where((sameev_vol_id == csi_copyno) & (sameev_part_id == 22))
						where_PMT_inc = np.where((sameev_vol_id == glass_copyno) & (sameev_part_id == -22) & (sameev_mdz_ent >= 0))
						#where_CsI = np.where((sameev_vol_id == csi_copyno) & (sameev_part_id == -22) & (sameev_ene_dep > 0))
						
						#edep_CsI = sameev_ene_dep[where_CsI]
						
						# PMT 
						if (where_PMT_inc[0].size):
							evt_id_PMT_inc = sameev_evt_id[where_PMT_inc]
							e_PMT_inc = sameev_ene_ent[where_PMT_inc]
							
							if (where_PMT[0].size):
								evt_id_PMT = sameev_evt_id[where_PMT]
								edep_PMT = sameev_ene_dep[where_PMT]
							
							#edep_tot = np.sum(edep_CsI) + np.sum(edep_PMT)
							edep = np.sum(edep_PMT)
							print(f"Beam n. {nbeam}, event n. {sameev_evt_id[0]}, E_dep in PMT = {edep:.2f} keV\n")
							if write_to_file == True:	
								with open("Energy_deposited_PMTnoCsI.txt", "w") as f:
									f.write(f"{edep}\n")
							
							e_ent_Gamma = sameev_ene_ent[where_Gamma]
							e_exit_Gamma = sameev_ene_exit[where_Gamma]
							trk_id_gamma = sameev_trk_id[where_Gamma]
							
							where_xray_fluo = np.where(trk_id_gamma > 1)
							
							#if (not where_xray_fluo[0].size & len(e_ent_Gamma == 1)): 

							for jp in range(len(evt_id_PMT_inc)):

								vecPMTEventID_inc.append(evt_id_PMT_inc[jp])
								vecPMTEent_inc.append(e_PMT_inc[jp])
								vecPMTNbeam_inc.append(nbeam)

								#vecGammaEent_inc.append(e_ent_Gamma[0])
								#vecGammaEexit_inc.append(e_exit_Gamma[0])

							if (where_PMT[0].size):
								vecPMTEventID_histo.append(evt_id_PMT[0])
								vecPMTEdep_histo.append(edep)
								vecPMTNbeam_histo.append(nbeam)

							#if (where_PMT[0].size and edep_tot > 9):	
							if (where_PMT[0].size):						
								for jp in range(len(evt_id_PMT)):
									vecPMTEventID.append(evt_id_PMT[jp])
									vecPMTEdep.append(edep_PMT[jp])
									vecPMTNbeam.append(nbeam)

									#vecGammaEent.append(e_ent_Gamma[0])
									#vecGammaEexit.append(e_exit_Gamma[0])				

								

		
						N_event_eq = len(where_sameevent)
						if (evt_id[where_sameevent[-1]+1] != evt_id[-1]):
							index = where_sameevent[N_event_eq-1] + 1
							where_sameevent = [index]
							temp_index = index
						else:
							first_sameevent = where_sameevent[-1]+1
							where_sameevent = np.arange(first_sameevent, len(evt_id), 1)
			
							sameev_evt_id = evt_id[where_sameevent]
							sameev_trk_id = trk_id[where_sameevent]
							sameev_vol_id = vol_id[where_sameevent]
							sameev_moth_id = moth_id[where_sameevent]
							sameev_ene_ent = ene_ent[where_sameevent]
							sameev_ene_exit = ene_exit[where_sameevent]
							sameev_ene_dep = e_dep[where_sameevent]
							sameev_mdx_ent = mdx_ent[where_sameevent]
							sameev_mdy_ent = mdy_ent[where_sameevent]
							sameev_mdz_ent = mdz_ent[where_sameevent]
							sameev_mdx_exit = mdx_exit[where_sameevent]
							sameev_mdy_exit = mdy_exit[where_sameevent]
							sameev_mdz_exit = mdz_exit[where_sameevent]
							sameev_x_ent = x_ent[where_sameevent]
							sameev_y_ent = y_ent[where_sameevent]
							sameev_z_ent = z_ent[where_sameevent]
							sameev_x_exit = x_exit[where_sameevent]
							sameev_y_exit = y_exit[where_sameevent]
							sameev_z_exit = z_exit[where_sameevent]
							sameev_part_id = part_id[where_sameevent]
							sameev_gtime_ent = gtime_ent[where_sameevent]
			
							where_PMT = np.where((sameev_vol_id == bialkali_copyno) & (sameev_part_id == -22) & (sameev_ene_dep > 0))
							where_Gamma = np.where((sameev_vol_id == csi_copyno) & (sameev_part_id == 22))
							where_PMT_inc = np.where((sameev_vol_id == glass_copyno) & (sameev_part_id == -22) & (sameev_mdz_ent >= 0))	
							#where_CsI = np.where((sameev_vol_id == csi_copyno) & (sameev_part_id == -22) & (sameev_ene_dep > 0))
						
							#edep_CsI = sameev_ene_dep[where_CsI]
							
							# PMT 
							if (where_PMT_inc[0].size):
								evt_id_PMT_inc = sameev_evt_id[where_PMT_inc]
								e_PMT_inc = sameev_ene_ent[where_PMT_inc]
							
								if (where_PMT[0].size):
									evt_id_PMT = sameev_evt_id[where_PMT]
									edep_PMT = sameev_ene_dep[where_PMT]	
									
								#edep_tot = np.sum(edep_CsI) + np.sum(edep_PMT)
								edep = np.sum(edep_PMT)
								#print(f"Beam n. {nbeam}, event n. {sameev_evt_id[0]}, E_dep in PMT = {edep:.2f} keV\n")					
								if write_to_file == True:	
									with open("Energy_deposited_PMTnoCsI.txt", "w") as f:
										f.write(f"{edep}\n")

								e_ent_Gamma = sameev_ene_ent[where_Gamma]
								e_exit_Gamma = sameev_ene_exit[where_Gamma]
								trk_id_gamma = sameev_trk_id[where_Gamma]
							
								where_xray_fluo = np.where(trk_id_gamma > 1)
							
								#if (not where_xray_fluo[0].size & len(e_ent_Gamma == 1)): 

								for jp in range(len(evt_id_PMT_inc)):

									vecPMTEventID_inc.append(evt_id_PMT_inc[jp])
									vecPMTEent_inc.append(e_PMT_inc[jp])
									vecPMTNbeam_inc.append(nbeam)

									#vecGammaEent_inc.append(e_ent_Gamma[0])
									#vecGammaEexit_inc.append(e_exit_Gamma[0])	
								
								if (where_PMT[0].size):
									vecPMTEventID_histo.append(evt_id_PMT[0])
									vecPMTEdep_histo.append(edep)
									vecPMTNbeam_histo.append(nbeam)

								#if (where_PMT[0].size and edep_tot > 9):
								if (where_PMT[0].size):							
									for jp in range(len(evt_id_PMT)):
										vecPMTEventID.append(evt_id_PMT[jp])
										vecPMTEdep.append(edep_PMT[jp])
										vecPMTNbeam.append(nbeam)

										#vecGammaEent.append(e_ent_Gamma[0])
										#vecGammaEexit.append(e_exit_Gamma[0])	
							
							break

				hdulist.close()

	# remove bad cases (Edep in PMT 2 sigma away)
	events_outliers = []
	media = np.mean(vecPMTEdep_histo)
	std = np.std(vecPMTEdep_histo)
	for kk in range(len(vecPMTEdep_histo)):
		if vecPMTEdep_histo[kk] < (media - 2*std) or vecPMTEdep_histo[kk] > (media + 2*std):
			events_outliers.append(vecPMTEventID_histo[kk])
	
	vecPMTNbeam_old = vecPMTNbeam
	vecPMTEdep = [vecPMTEdep[i] for i in range(len(vecPMTEdep)) if (vecPMTEventID[i] not in events_outliers or vecPMTNbeam[i] != nbeam)]
	vecPMTNbeam = [vecPMTNbeam[i] for i in range(len(vecPMTNbeam)) if (vecPMTEventID[i] not in events_outliers or vecPMTNbeam[i] != nbeam)]
	vecPMTEventID = [vecPMTEventID[i] for i in range(len(vecPMTEventID)) if (vecPMTEventID[i] not in events_outliers or vecPMTNbeam_old[i] != nbeam)]
	#print(events_outliers)


vecPMTEventID = np.array(vecPMTEventID)
vecPMTEdep = np.array(vecPMTEdep)
vecPMTNbeam = np.array(vecPMTNbeam)
#vecGammaEent = np.array(vecGammaEent)
#vecGammaEexit = np.array(vecGammaEexit)

vecPMTEventID_inc = np.array(vecPMTEventID_inc)
vecPMTEent_inc = np.array(vecPMTEent_inc)
vecPMTNbeam_inc = np.array(vecPMTNbeam_inc)
#vecGammaEent_inc = np.array(vecGammaEent_inc)
#vecGammaEexit_inc = np.array(vecGammaEexit_inc)

# Write output
if not os.path.exists(path_output):
	os.makedirs(path_output)


# Open and write file 
f_out = open(path_output+"PMT.dat", 'w')
print('... Writing '+path_output+"PMT.dat")

f_out.write("# EventID N_BEAM Edep[keV] \n")
for je in range(len(vecPMTEventID)):

	f_out.write('{0:30}'.format(str(vecPMTEventID[je])))
	f_out.write('{0:30}'.format(str(vecPMTNbeam[je])))
	f_out.write('{0:30}'.format(str(vecPMTEdep[je])))
	#f_out.write('{0:30}'.format(str(vecGammaEent[je])))
	#f_out.write('{0:30}'.format(str(vecGammaEexit[je])))

	f_out.write("\n")

f_out.close()