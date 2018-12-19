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
from sei_functions import df_2spec2var
from sei_functions import CV_voltage
import cantera as ct
from assimulo.solvers import IDA
from assimulo.problem import Implicit_Problem

ctifile = 'W_anode_chem.cti'

# Phase names in cti file
elyte_phase = 'electrolyte'
sei_phase = 'SEI'
sei_conductorphase = 'conductor'
WE_phase = 'tungsten'
WE_elyte_surfphase = 'tungsten_electrolyte_surf'
WE_sei_surfphase = 'tungsten_SEI_surf'
sei_elyte_surfphase = 'SEI_electrolyte_surf'
CE_phase = 'Lithium'
CE_surfphase = 'Li_surf'

# %% Define grid variables

"""----------Define grid dimensions----------"""

# A (2,10) grid will be indexed as shown below with each node having the
# number of tracked variables stored for it. The solution vector will be
# a row vector progressing through each tracked variable for each cell, then
# progressing in increasing column number for increasing row number.
#
#          |--> electrolyte
#          |
#          |-----------------------------------------------------------|
#          | 0,0 | 0,1 | 0,2 | 0,3 | 0,4 | 0,5 | 0,6 | 0,7 | 0,8 | 0,9 |
# electrode|-----------------------------------------------------------|
#          | 1,0 | 1,1 | 1,2 | 1,3 | 1,4 | 1,5 | 1,6 | 1,7 | 1,8 | 1,9 |
#          |-----------------------------------------------------------|
#
# The properties of species will be stored in the order shown below for each
# cell.
#
# [T, V_elyte, V_SEI_an, V_SEI_elyte, c_SEI_1, ... c_SEI_k, c_elyte_1, ... c_elyte_k, phi]
#

N_x = 1         # USER INPUT number of grids in plane of electrode
N_y = 50         # USER INPUT number of grids perpendicular to electrode
x = 1e-7       # USER INPUT x length of domain [m]
y = 1           # USER INPUT y length of domain [m]

# %% Define state variables

"""----------Initialize what state variables are tracked----------"""
track_temp = 0          # USER INPUT 1 turns temperature tracking on
track_vol_frac = 1      # USER INPUT 1 turns volume fraction tracking on
T_0 = np.array([300.])              # USER INPUT sets an initial temperature [K] assuming it's uniform
P_0 = 101325.

eps_0 = np.array([0.])              # initial volume fraction of SEI

# total number of state variables tracked
num_state_vars = track_temp + track_vol_frac

# %% Define phase variables

"""----------Define number of phase variables to track----------"""

# ---------------------------------------------------------------------------
#   Initially we are tracking 1 variable per phase. If more are desired change
#   track_phase_vars value accordingly and update variable list commented
#   below.
# ---------------------------------------------------------------------------

"""
    1. electric potential of electrolyte [V]
    2. electric potential of SEI/anode interface [V]
    3. electric potential of SEI/electrolyte interface [V]
"""

sweep_rate = 0.01  #...Voltage sweep rate [V/s]
sweep_dirn_0 = -1
phi_0 = 1.0;
phi_1 = 0.05;
phi_2 = 1.5;
n_cycles = 1.

check_profile = 0
#
#...For these half-cell simulations, let us assume a constant cathode
#       voltage for our energy storage and efficiency calcs:
V_CE = 0.0; #...(V)yte interface

phi_SEI_0 = 0.
phi_elyte_0 = 0.

# %% Define species and species variables


"""
=========================================================================================
==================================END USER INPUTS========================================
=========================================================================================
"""

global objs
global params
global voltage_lookup
global SVptr


"""----------Geometry calcs----------"""
dx = x/N_x      # USER INPUT length of step in x direction
dy = y/N_y      # USER INPUT length of step in y direction

"""----------Create Cantera objects----------"""
CE, elyte, sei, sei_conductor, WE = ct.import_phases(ctifile, \
    [CE_phase, elyte_phase, sei_phase, sei_conductorphase, WE_phase])
WE_elyte = ct.Interface(ctifile, WE_elyte_surfphase, [WE, elyte, sei])
WE_sei =  ct.Interface(ctifile, WE_sei_surfphase, \
    [WE, sei, sei_conductor])
sei_elyte = ct.Interface(ctifile, sei_elyte_surfphase, \
    [sei, elyte, sei_conductor])
CE_elyte = ct.Interface(ctifile, CE_surfphase, [CE, elyte])

print('\n Cantera phases created. \n')

# %% Set phase properties for cantera objects
TP_o = T_0, P_0
elyte.TP = TP_o
sei.TP = TP_o
WE.TP = TP_o
CE.TP = TP_o
WE_elyte.TP = TP_o
WE_sei.TP = TP_o
sei_elyte.TP = TP_o
CE_elyte.TP = TP_o

elyte.electric_potential = phi_elyte_0
sei.electric_potential = phi_SEI_0
WE.electric_potential = phi_0

# Store these objects in a common dictionary:
objs = {'WE':WE, 'SEI':sei, 'elyte':elyte, 'CE':CE, 'WE_SEI':WE_sei, \
    'WE_elyte':WE_elyte, 'SEI_elyte':sei_elyte, 'CE_Elyte':CE_elyte}


"""----------Define electrolyte species----------"""

# ---------------------------------------------------------------------------
# The electrolyte information is brought in from an input file using Cantera.
# ---------------------------------------------------------------------------

num_elyte_species = elyte.n_species
C_k_elyte = elyte.density_mole*elyte.X

print("The species in the electrolyte are", elyte.species_names,'\n')


"""----------Define SEI species----------"""

# ---------------------------------------------------------------------------
# The SEI information is brought in from an input file using Cantera.
# ---------------------------------------------------------------------------

num_SEI_species = sei.n_species
C_k_sei = sei.density_mole*sei.X

print("The species in the SEI are\n", sei.species_names)

# %% Set up solution vector dimensions

"""----------Set initial values for potential and concentration---------- """

# ---------------------------------------------------------------------------
# We're going to assume an initial concentration of the electrolyte and no
# concentration of the SEI species.
# ---------------------------------------------------------------------------

# %% Initializing solution vector

"""----------Set up solution vector pointer SVptr----------"""

SVptr = {}

"""     sei/elyte domain    """
# 1 variable for sei volume fraction
# 1 variable for sei electric potential
# 1 variable for elyte electric potential
# 1 variable for each species in sei and elyte
nvars_node = 3 + sei.n_species + elyte.n_species

SVptr['phi sei'] = np.arange(0,nvars_node*N_y,nvars_node)
SVptr['phi elyte'] = np.arange(1,nvars_node*N_y,nvars_node)
SVptr['eps sei'] = np.arange(2,1+nvars_node*N_y,nvars_node)
SVptr['Ck sei'] = np.ndarray(shape=(N_y,sei.n_species))
SVptr['Ck elyte'] = np.ndarray(shape=(N_y,elyte.n_species))
for i in np.arange(N_y):
    SVptr['Ck sei'][i,:] = np.arange(3+i*nvars_node,\
        3+i*nvars_node+sei.n_species)
    SVptr['Ck elyte'][i,:] = np.arange(3+i*nvars_node+sei.n_species,\
        3+i*nvars_node+sei.n_species+elyte.n_species)

# ---------------------------------------------------------------------------
# Now use the total number of tracked variables per node to calculate the
# length of the solution vector based on discretization and tracked variables
# ---------------------------------------------------------------------------
nvars_tot = N_x*(1 + N_y*nvars_node)

# ---------------------------------------------------------------------------
# Initialize solution vector using np.zeros based on number of the variables
# ---------------------------------------------------------------------------
SV_node = np.concatenate((np.array((phi_SEI_0, phi_elyte_0, eps_0)), C_k_sei, C_k_elyte))
SV_0 = np.tile(SV_node, N_y)
SV_dot_0 = np.zeros_like(SV_0)
res = np.zeros_like(SV_0)

t_0 = 0

# Set preliminary parameters for the anode voltage sweep function
phi_bounds = np.array([phi_1, phi_2])  # Upper and lower voltage bounds
R = sweep_rate                            # Sweep rate [V/s]

# Times for discontinuities (sweep sign change) in anode voltage
t_event0 = -sweep_dirn_0*(phi_0 - phi_bounds[int(0.5*(1. + sweep_dirn_0))])/(R)
dt = (phi_bounds[1] - phi_bounds[0])/R
t_events = np.arange(t_event0,t_event0+dt*(2*n_cycles+1),dt)
times = np.concatenate((np.array([0.]),t_events))

voltage_array = np.zeros_like(times)
voltage_array[0] = phi_0

direction = sweep_dirn_0
for i, t in enumerate(t_events):
    voltage_array[i+1] = voltage_array[i] + direction*(t - times[i])*R
    direction *= -1


t = np.linspace(0,times[-1],500)
v = np.interp(t,times,voltage_array)

if check_profile:
    print('Check that voltage profile is correct.  Cancel simulation, if not. \n')
    plt.figure()
    plt.plot(t,v)
    plt.xlabel('time (s)')
    plt.ylabel(r"$\Phi_{WE} - \Phi_{CE}$" ' (V)')
    plt.show()

params = {'phi bounds':phi_bounds, 'Rate':R, 'Ny':N_y, 'dyInv':1./dy, 'TP':TP_o}
voltage_lookup = {'time':times, 'voltage':voltage_array}
"""----------Run solver----------"""

"""----------Residual function for IDA solver----------"""
def residual(t, SV, SV_dot):

    # Tag global variables in residual function
    global objs
    global params
    global voltage_lookup
    global SVptr

    # Initialize residual equal to all zeros:
    res = SV_dot

    # Read out cantera objects:
    WE = objs['WE']
    sei = objs['SEI']
    elyte = objs['elyte']

    # Set temperature, density, and concentration of species
    #elyte.TPX = params['TP'], SV[]

    # This block accounts for the anode voltage based ON TIME
    """if t <= t_event0 :
        phi_anode = phi_anode_0 + R*t
    elif t > t_event0 and t <= t_event1:
        phi_anode = phi_anode_max_0 - R*t
    elif t > t_event1 and t <= t_event2:
        phi_anode = phi_anode_min_0 + R*t
    elif t > t_event2 and t <= t_event3:
        phi_anode = phi_anode_max_1 - R*t"""
    phi_WE = np.interp(t,voltage_lookup['time'],voltage_lookup['voltage'])
    #print(phi_WE)


    # Set anode and anode/electrolyte interface potentials
    WE.electric_potential = phi_WE
    #anode_elyte.electric_potential = phi_anode
    #elyte.electric_potential = 0


    # Retreive net production rates of species from cantera
#    S_dot_elyte = anode_elyte.net_production_rates[elyte_rate_ind]
    #S_dot_SEI = anode_elyte.net_production_rates[SEI_rate_ind]

    #if phi_anode >= 0.69 and phi_anode <= 0.7:
    #    print(phi_anode, S_dot_SEI)

    """res[0] = SV_dot[0]
    res[1] = SV_dot[1]
    res[2] = SV_dot[2]
    res[3] = SV_dot[3]

    if 0: #for i in Elyte_vars_range:
        # Loop for electrolyte concentrations
        res[i] = SV_dot[i]
        #elyte_ctr += 1
    if 0: #for i in SEI_vars_range:
        # Loop for SEI concentrations
        res[i] = SV_dot[i] - S_dot_SEI[SEI_ctr]/dx
        #SEI_ctr += 1

    # Volume fraction of SEI per cell
    res[-1] = SV_dot[-1]"""

    return res

"""------------------------------------------------------------------------"""

pass_in = (objs, params, voltage_lookup, SVptr)

# Set up problem instance
SEI_1D = Implicit_Problem(residual, SV_0, SV_dot_0)

# Define simulation parameters
simulation = IDA(SEI_1D)                # Create simulation instance
simulation.atol = 1e-6                  # Solver absolute tolerance
simulation.rtol = 1e-6                  # Solver relative tolerance
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

fig, ax1 = plt.subplots(1, 1, figsize=(10, 9))
ax1.plot(t,SV[:,SVptr['Ck elyte'][0].astype(int)])

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

plt.show()

#SV_plot = SV_df.plot()
#SV_plot.legend(loc = 'upper left')
print("\nGoodbye")

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
