# -*- coding: utf-8 -*-
"""
Created on Wed May 16 14:52:05 2018
@author: Steven C. DeCaluwe, Daniel Korff

This code is the foundation (top level script) for a model of solid-electrolye-
interface (SEI) formation in a Li-ion battery (LIB). This is a preliminary
model that will consider one species in the SEI, one species in the
electrolyte, and track temperature, concentration, and electrical potential of
the species in a discretized grid employing a finite volume method.
"""

# %% Modules imported
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sei_1d_functions import residual
from sei_1d_functions import df_2spec2var
import cantera as ct
from assimulo.solvers import IDA
from assimulo.problem import Implicit_Problem

print('\n     Importing inputs and intializing.')
from sei_1d_init import SV_0, SV_dot_0, SVptr, times, objs, params,  \
    voltage_lookup

print('\n     Running simulation\n')
# Set up problem instance
SEI_1D = Implicit_Problem(residual, SV_0, SV_dot_0)

# Define simulation parameters
simulation = IDA(SEI_1D)                # Create simulation instance
simulation.atol = 1e-7                  # Solver absolute tolerance
simulation.rtol = 1e-4                  # Solver relative tolerance
#simulation.maxh = 0.1                   # Solver max step size

# Set simulation end time, slope flag (for anode voltage cycle), and run simulation

t_f = times[-1]

#ncp_list = np.arange(0, t_f, 0.15)
#ncp = 10000

# Run simulation
t, SV, SV_dot = simulation.simulate(t_f)


# %% Organize data and plot

#SV_df = pd.DataFrame(SV)
#SV_df = df_2spec2var(SV_df,N_x,N_y,len_sol_vec,track_vars,track_temp,num_species)

#fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 9))
sei = objs['SEI']
names = list()
for i in range(1,sei.n_species):
    names.append(sei.species_names[i])
names.append('eps_sei')
names.append('Anode potential')
names.append('SEI potential')

phi_WE = np.interp(t,voltage_lookup['time'],voltage_lookup['voltage'])
phi_SEI = phi_WE + SV[:,SVptr['phi sei'][0]]

fig, ax1 = plt.subplots(1, 1, figsize=(10, 9))
ax1.plot(t,SV[:,SVptr['Ck sei'][0,1:].astype(int)])
ax1.plot(t,SV[:,SVptr['eps sei'][0]])
ax1.plot(t,phi_WE)
ax1.plot(t,phi_SEI)
ax1.legend(names)
ax1.set_ylabel('Molar concentration (kmol/m3), Vol fraction, Electric potential (V)')
ax1.set_xlabel('time (s)')
#plt.show()
"""plt.savefig('Figure1.pdf',format='pdf',dpi=350)"""

profiles = SV[-1,SVptr['Ck sei'][:,1:]]
fig2, ax2 = plt.subplots(1, 1, figsize=(10, 9))
ax2.plot(1e9*np.arange(params['Ny'])/params['dyInv'],profiles)
ax2.plot(1e9*np.arange(params['Ny'])/params['dyInv'],SV[-1,SVptr['eps sei']])
ax2.legend(names)
ax2.set_ylabel('Molar concentration (kmol/m3)')
ax2.set_xlabel('SEI Depth (from anode, nm)')
plt.show()
"""plt.savefig('Figure2.pdf',format='pdf',dpi=350)"""



#ax3.plot(t, SV_df['c_elyte1_00'], label = elyte_species[0, 0])
#ax3.plot(t, SV_df['c_elyte2_00'], label = elyte_species[0, 1])
#ax3.plot(t, SV_df['c_elyte3_00'], label = elyte_species[0, 2])
#ax3.plot(t, SV_df['c_elyte4_00'], label = elyte_species[0, 3])
#ax3.plot(t, SV_df['c_elyte5_00'], label = elyte_species[0, 4])
#ax3.plot(t, SV_df['c_elyte6_00'], label = elyte_species[0, 5])
#ax3.plot(t, SV_df['c_elyte7_00'], label = elyte_species[0, 6])
#ax3.plot(t, SV_df['c_elyte8_00'], label = elyte_species[0, 7])
#ax3.plot(t, SV_df['c_elyte9_00'], label = elyte_species[0, 8])

"""ax1.plot(t, SV_df['c_SEI2_00'], '-b', label = SEI_species[0, 1])
#ax1.plot(t, SV_df['c_SEI3_00'], '-r', label = SEI_species[0, 2])
ax1.plot(t, SV_df['c_SEI4_00'], '-g', label = SEI_species[0, 3])
ax1.set_ylabel('Concentration ' r"[\frac{mol}{m^3}]")
ax1.set_title('SEI species concentrations over time')
ax1.legend(loc='upper left')

ax2.plot(t_vec, phi_vec, label = 'Anode Voltage')
ax2.set_ylabel('Anode Potential [V]')
ax2.set_xlabel('Time [s]')
ax2.set_title('Anode voltage over time')"""




#SV_plot = SV_df.plot()
#SV_plot.legend(loc = 'upper left')
print("\n     Goodbye")

# %% Functions !!!!!NOT CURRENTLY USED!!!!!


"""----------Time events function for IDA solver----------"""
def time_events(self, t, SV, SV_dot, sw):

    if t < t_event0:
        t_event = t_event0
    elif t >= t_event0 and t < t_event1:
        t_event = t_event1
    else:
        t_event = None

    return t_event

"""------------------------------------------------------------------------"""

"""----------Time events function for IDA solver----------"""
def handle_event(self, solver, event_info):

    event_info = event_info[0]
    while True:
        self.event_switch(solver, event_info)

        b_mode = self.time_events(solver.t, solver.SV, solver.SV_dot, solver.sw)
        self.init_mode(solver)
        a_mode = self.time_events(solver.t, solver.SV, solver.SV_dot, solver.sw)

        event_info = self.checkeIter(b_mode, a_mode)

        if not True in event_info:
            break



"""------------------------------------------------------------------------"""
