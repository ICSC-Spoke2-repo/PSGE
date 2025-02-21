"""
 plotCLAIRE.py  -  description
 ---------------------------------------------------------------------------------
 plotting the events of COSI CLAIRE experiment
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


import numpy as np
import math
import sys
import matplotlib.pyplot as plt

from matplotlib import gridspec
from matplotlib import cm
from astropy import units as u


cmap = cm.get_cmap('jet')

# Import the input parameters
arg_list = sys.argv
filepath_input = []
for a in arg_list[1:]:
	filepath_input.append(arg_list[1])

cases_label = ["Simulation"]
cases_color = ["red"]

ncases = len(filepath_input)

N_beams = 40

# entering
vecPMTEventID = []
vecPMTEdep = []
vecPMTNbeam = []

vecGammaEent = []
vecGammaEexit = []	

minimum_list = []
maximum_list = []

fig1 = plt.figure(1, figsize=(10, 7))
gs = gridspec.GridSpec(5, 1)
gs.update(hspace = 0.)
ax1 = plt.Subplot(fig1, gs[0:4, 0:1])
fig1.add_subplot(ax1)
ax11 = plt.Subplot(fig1, gs[4:5, 0:1])
fig1.add_subplot(ax11)


vecPulse_real = []
f_pulse = open("./config/CLAIRE_pulse.dat", 'r')
for line in f_pulse:
	line = line.strip()
	if not line.startswith("#"):
		columns = line.split()
		columns[1] = float(columns[1])
	
		vecPulse_real.append(columns[1])

f_pulse.close()

for jc in range(ncases):
	print("%%%%%%%%% Processing "+cases_label[jc]+" ...")
	# entering
	vecPMTEventID = []
	vecPMTEdep = []
	vecPMTNbeam = []
	vecGammaEent = []
	vecGammaEexit = []
		
	print("Reading file ...\n")
	f_g4 = open(filepath_input[jc], 'r')
	for line in f_g4:
		line = line.strip()
		if not line.startswith("#"):
			columns = line.split()
			columns[0] = int(columns[0])
			columns[1] = int(columns[1])
			columns[2] = float(columns[2])
			
			vecPMTEventID.append(columns[0])
			vecPMTNbeam.append(columns[1])
			wave_temp = ((columns[2]*1000.) * u.eV).to(u.nm, equivalencies=u.spectral())
			vecPMTEdep.append(wave_temp.value) #nm

	f_g4.close()
	
	vecPMTEventID = np.array(vecPMTEventID)
	vecPMTEdep = np.array(vecPMTEdep)
	vecPMTNbeam = np.array(vecPMTNbeam)

	vecResponse = []
	vecErrResponse = []

	vecEdep = []
	vecErrEdep = []

	vecPulse = []
	vecErrPulse = []
		
	for jbeam in range(N_beams):
		nbeam = jbeam + 1
		##### N ABSORBED
		vecPMTEventID_sel = vecPMTEventID[np.where(vecPMTNbeam == nbeam)]
		vecPMTEdep_sel = vecPMTEdep[np.where(vecPMTNbeam == nbeam)]
		n_abs = len(np.unique(vecPMTEventID_sel))
		Eff_tot = len(vecPMTEdep_sel)
		Eff = Eff_tot/n_abs

		Edep_tot = np.sum(vecPMTEdep_sel)
		Edep = Edep_tot/n_abs
		err_Eff = np.sqrt((np.sqrt(Eff_tot)/n_abs)**2+(Eff_tot*np.sqrt(n_abs)/n_abs**2)**2)
		err_Edep = Edep_tot*np.sqrt(n_abs)/n_abs**2
		
		vecResponse.append(Eff)
		vecErrResponse.append(err_Eff)
		vecEdep.append(Edep)
		vecErrEdep.append(err_Edep)
	
	vecResponse = np.array(vecResponse)
	vecResponse_tot = np.sum(vecResponse)	
	vecErrResponse = np.array(vecErrResponse)
	vecRelResponse = vecResponse/vecResponse_tot
	vecErrRelResponse = vecErrResponse/vecResponse_tot

	vecEdep = np.array(vecEdep)
	vecErrEdep = np.array(vecErrEdep)


	ax1.errorbar(range(1, N_beams+1), vecRelResponse, yerr=vecErrRelResponse, capsize=0, fmt='-o', lw = 1, ms=5, zorder=10, color=cases_color[jc], ecolor=cases_color[jc], label = cases_label[jc])
	ax11.errorbar(range(1, N_beams+1), (vecRelResponse-vecPulse_real)*100/vecPulse_real, yerr=vecErrRelResponse*100/vecPulse_real, capsize=0, fmt='o', lw = 1, ms=5, zorder=10, color=cases_color[jc], ecolor=cases_color[jc], label = cases_label[jc])
	ax11.errorbar([0, 100], [0, 0], fmt = 'k--')

	name_file = "data.txt"
	with open(name_file, "w") as f:
		f.write("# Beam RelResponse ErrRelResponse RealResponse Edep ErrEdep\n")
		for i in range(len(vecRelResponse)):
			f.write(f"{i+1} {vecRelResponse[i]} {vecErrRelResponse[i]} {vecPulse_real[i]} {vecEdep[i]} {vecErrEdep[i]}\n")

	minim = min((vecRelResponse-vecPulse_real)*100/vecPulse_real - vecErrRelResponse*100/vecPulse_real)
	maxim = max((vecRelResponse-vecPulse_real)*100/vecPulse_real + vecErrRelResponse*100/vecPulse_real)
	minimum_list.append(math.floor(minim/10)*10)
	maximum_list.append(math.ceil(maxim/10)*10)
	
ax1.plot(range(1, N_beams+1), vecPulse_real, '-ko', label = "experiment")

ax1.set_title(r"Geant4 vs experiment (122 keV, GBP, SS reflection, $\sigma=40^\circ$)", fontsize=12)
ax1.set_ylabel(r'PMT relative response', fontsize=12)
ax11.set_ylabel(r'(sim - obs) / obs (%)', fontsize=8)
ax1.set_xticklabels([])
ax11.set_xlabel(r'Beam position', fontsize=12)
ax1.set_ylim(0.01, 0.06)
ax1.set_xlim(0, N_beams+1)
ax11.set_xlim(0, N_beams+1)
ax1.grid(True)
nstep = 5
minimum = min(minimum_list)
maximum = max(maximum_list)
stepsize = (maximum - minimum) / nstep
ax11.set_ylim(minimum, maximum)
ax11.yaxis.set_ticks(np.arange(minimum, maximum, stepsize))
ax11.grid(True)
ax1.legend(numpoints=1, fontsize = 9)

plt.savefig("plot.pdf")
