# -*- coding: utf-8 -*-
"""
Created on Wed May 16 14:52:05 2018
@author: Daniel Korff

This code is the foundation (top level script) for a model of solid-electrolye-
interface (SEI) formation in a Li-ion battery (LIB). This is a preliminary 
model that will consider one species in the SEI, one species in the 
electrolyte, and track temperature, concentration, and electrical potential of 
the species in a discretized grid employing a finite volume method.
"""

# %% Modules imported

print("\nHello!")

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from SEI_prelim_functions import df_2spec2var
from SEI_prelim_functions import IC_fun
import cantera as ct
from assimulo.solvers import IDA
from assimulo.problem import Implicit_Problem

cantera_file = 'W_anode_chem.cti'
elyte_name = 'electrolyte'
SEI_name = 'SEI'
anode_name = 'tungsten'
anode_elyte_surf = 'tungsten_electrolyte_surf'

elyte, SEI, anode = ct.import_phases(cantera_file, \
                                     [elyte_name, SEI_name, anode_name])
anode_elyte = ct.Interface(cantera_file, anode_elyte_surf, [anode, elyte, SEI])


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
N_y = 1         # USER INPUT number of grids perpendicular to electrode
x = 10e-6       # USER INPUT x length of domain [m]
y = 1           # USER INPUT y length of domain [m]
dx = x/N_x      # USER INPUT length of step in x direction
dy = y/N_y      # USER INPUT length of step in y direction

# %% Define state variables

"""----------Initialize what state variables are tracked----------"""
track_temp = 1        # USER INPUT 1 turns temperature tracking on
track_vol_frac = 1    # USER INPUT 1 turns volume fraction tracking on
T_0 = np.array([300])   # USER INPUT sets an initial temperature [K] assuming it's uniform
P_0 = 101325.

eps_0 = np.array([0])           # initial volume fraction of SEI

# total number of state variables tracked
num_state_vars = track_temp + track_vol_frac

# %% Set phase properties for cantera objects

elyte.TP = T_0, P_0
SEI.TP = T_0, P_0
anode.TP = T_0, P_0
anode_elyte.TP = T_0, P_0

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

track_phase_vars = 3    # USER INPUT number of phase variables tracked

phi_elyte_0 = 0.0           # USER INPUT initial potential of electrolyte phase
phi_SEI_an_0 = 0.0          # USER INPUT initial potential of SEI/anode interface
phi_SEI_elyte_0 = 0.0       # USER INPUT initial potential of SEI/electrolyte interface

phi_0_vec = np.array([phi_elyte_0, phi_SEI_an_0, phi_SEI_elyte_0])

# %% Define species and species variables

"""----------Define number of species variables to track----------"""

# ---------------------------------------------------------------------------
#   Initially we are tracking 1 variable per species in the electrolyte
#   and SEI. If more are desired change track_species_vars value accordingly
#   and update variable list commented below. k subscript used for kth species
#   note that the variables will be stored in the order below in the solution
#   vector
# ---------------------------------------------------------------------------

""" 
    1. Concentration [mol_k/m^3]
"""
track_species_vars = 1  # USER INPUT number of species variables tracked

""" 
=========================================================================================
==================================END USER INPUTS========================================
=========================================================================================
"""

"""----------Define electrolyte species----------"""

# ---------------------------------------------------------------------------
# The electrolyte information is brought in from an input file using Cantera.
# ---------------------------------------------------------------------------

num_elyte_species = elyte.n_species
rho_elyte = elyte.density_mass
Y_elyte_vec = elyte.Y
c_elyte_vec = np.zeros(num_elyte_species)
c_elyte_vec = rho_elyte*Y_elyte_vec
    
print("\nThe species in the electrolyte are", elyte.species_names)

elyte_species = np.array([elyte.species_names])

"""----------Define SEI species----------"""

# ---------------------------------------------------------------------------
# The SEI information is brought in from an input file using Cantera. 
# ---------------------------------------------------------------------------

num_SEI_species = SEI.n_species
rho_SEI = SEI.density_mass
Y_SEI_vec = SEI.Y
c_SEI_vec = np.zeros(num_SEI_species)
c_SEI_vec = rho_SEI*Y_SEI_vec
    
print("The species in the SEI are\n", SEI.species_names[1:])

SEI_species = np.array([SEI.species_names])

# %% Set up solution vector dimensions

"""----------Set initial values for potential and concentration---------- """

# ---------------------------------------------------------------------------
# We're going to assume an initial concentration of the electrolyte and no
# concentration of the SEI species.
# ---------------------------------------------------------------------------

phi_0 = phi_0_vec
c_elyte_0 = c_elyte_vec   # Initial concentration of electrolyte [mol_k/m^3]
c_SEI_0   = c_SEI_vec     # Initial concentration of SEI species [mol_k/m^3]

IC_vec = np.concatenate((T_0, phi_0_vec, c_elyte_0, c_SEI_0, eps_0))

"""----------Define number of species for SEI and electrolyte----------"""

num_species       = np.array([num_elyte_species, num_SEI_species])

track_elyte_vars  = track_species_vars * num_elyte_species
track_SEI_vars    = track_species_vars * num_SEI_species

# %% Initializing solution vector

"""----------Set up solution vector----------"""

elyte_ptr = track_temp + track_phase_vars
SEI_ptr = elyte_ptr + num_elyte_species

Gen_vars_range = np.arange(0, track_temp+track_phase_vars, 1)
Elyte_vars_range = np.arange(elyte_ptr, elyte_ptr+num_elyte_species, 1)
SEI_vars_range = np.arange(SEI_ptr, SEI_ptr+num_SEI_species, 1)

# ---------------------------------------------------------------------------
# Calculate total number of tracked variables. This is how many variables
# will need to be stored for each node in the grid
# ---------------------------------------------------------------------------
track_vars = num_state_vars + track_elyte_vars + track_SEI_vars + track_phase_vars

# ---------------------------------------------------------------------------
# Now use the total number of tracked variables per node to calculate the
# length of the solution vector based on discretization and tracked variables
# ---------------------------------------------------------------------------
len_sol_vec = N_x * N_y * track_vars

# ---------------------------------------------------------------------------
# Initialize solution vector using np.zeros based on number of the variables
# ---------------------------------------------------------------------------
SV = np.zeros([1, len_sol_vec])


 
# %% Run solver and process data
"""----------Set up boundary conditions for solution vector----------"""

global phi_anode
global phi_anode_0
global R
global phi_bounds
global phi_vec
global phi_anode_vec
global ncp_list
global phi_ctr

# Initialize the solution vector, its derivative, time, and the residual vector
SV_0 = IC_fun(SV, IC_vec, len_sol_vec, track_vars, track_temp)
SV_dot_0 = SV_0
t_0 = 0
res = np.zeros([len_sol_vec])

# Set pointers for the production rates in the net production rates vector
# that Cantera hands back out based on number of species in the phases
elyte_rate_ptr = anode_elyte.n_species
SEI_rate_ptr = elyte_rate_ptr + num_elyte_species

# Create vectors of the indices that are filled by respective phase production
# rates
elyte_rate_ind = np.arange(elyte_rate_ptr, elyte_rate_ptr+num_elyte_species, 1)
SEI_rate_ind = np.arange(SEI_rate_ptr, SEI_rate_ptr+num_SEI_species, 1)

# Set preliminary parameters for the anode voltage sweep function
phi_bounds = np.array([0.05, 1.5])  # Upper and lower voltage bounds
R = 0.05                            # Sweep rate [V/s]
phi_anode_0 = np.mean(phi_bounds)   # Initial voltage of anode [V]

# Times for discontinuities (sweep sign change) in anode voltage
t_event0 = (phi_bounds[1] - phi_anode_0)/(R)
t_event1 = (phi_bounds[1] - phi_bounds[0])/R + t_event0
t_event2 = (phi_bounds[1] - phi_bounds[0])/R + t_event1

phi_anode_max_0 = phi_bounds[1] + t_event0*R
phi_anode_min_0 = phi_bounds[0] - t_event1*R

phi_ctr = 0
phi_vec = []
t_vec = []

"""----------Run solver----------"""

"""----------Residual function for IDA solver----------"""
def residual(t, SV, SV_dot):
        
    # Tag global variables in residual function
    global phi_anode_0
    global phi_anode
    global phi_vec
    global phi_ctr
    global t_vec
    
    # Counter for production rate species 
    elyte_ctr = 0
    SEI_ctr = 0
    
    # Set temperature, density, and concentration of species
    elyte.TDX = SV[0], rho_elyte, SV[elyte_ptr:elyte_ptr+num_elyte_species]
        
    # This block accounts for the anode voltage based ON TIME    
    if t <= t_event0 :
        phi_anode = phi_anode_0 + R*t
    elif t > t_event0 and t <= t_event1:
        phi_anode = phi_anode_max_0 - R*t
    elif t > t_event1 and t <= t_event2:
        phi_anode = phi_anode_min_0 + R*t
        
    phi_vec.append(phi_anode)
    t_vec.append(t)
    
    # Set anode and anode/electrolyte interface potentials
    anode.electric_potential = phi_anode
    anode_elyte.electric_potential = phi_anode
    
    # Retreive net production rates of species from cantera
#    S_dot_elyte = anode_elyte.net_production_rates[elyte_rate_ind]
    S_dot_SEI = anode_elyte.net_production_rates[SEI_rate_ind]
        
    res[0] = SV_dot[0]
    res[1] = SV_dot[1]
    res[2] = SV_dot[2]
    res[3] = SV_dot[3]
    
    for i in Elyte_vars_range:
        # Loop for electrolyte concentrations
        res[i] = SV_dot[i]
        elyte_ctr += 1
    for i in SEI_vars_range:
        # Loop for SEI concentrations
        res[i] = SV_dot[i] - S_dot_SEI[SEI_ctr]/dx
        SEI_ctr += 1
    
    elyte_ctr = 0
    SEI_ctr = 0
    phi_ctr += 1
    
    # Volume fraction of SEI per cell
    res[-1] = SV_dot[-1]  
    
    return res

"""------------------------------------------------------------------------"""

# Set up problem instance
SEI_1D = Implicit_Problem(residual, SV_0, SV_dot_0, t_0)

# Define simulation parameters
simulation = IDA(SEI_1D)                # Create simulation instance
simulation.atol = 1e-6                  # Solver absolute tolerance
simulation.rtol = 1e-6                  # Solver relative tolerance
simulation.maxh = 0.1                   # Solver max step size

# Set simulation end time, slope flag (for anode voltage cycle), and run simulation

t_f = ((phi_bounds[1] - phi_anode_0)/R)*5

#ncp_list = np.arange(0, t_f, 0.15)
#ncp = 10000
    
# Run simulation
t, SV, SV_dot = simulation.simulate(t_f)
  
    
# %% Organize data and plot

SV_df = pd.DataFrame(SV)
SV_df = df_2spec2var(SV_df,N_x,N_y,len_sol_vec,track_vars,track_temp,num_species)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 9))

#ax3.plot(t, SV_df['c_elyte1_00'], label = elyte_species[0, 0])
#ax3.plot(t, SV_df['c_elyte2_00'], label = elyte_species[0, 1])
#ax3.plot(t, SV_df['c_elyte3_00'], label = elyte_species[0, 2])
#ax3.plot(t, SV_df['c_elyte4_00'], label = elyte_species[0, 3])
#ax3.plot(t, SV_df['c_elyte5_00'], label = elyte_species[0, 4])
#ax3.plot(t, SV_df['c_elyte6_00'], label = elyte_species[0, 5])
#ax3.plot(t, SV_df['c_elyte7_00'], label = elyte_species[0, 6])
#ax3.plot(t, SV_df['c_elyte8_00'], label = elyte_species[0, 7])
#ax3.plot(t, SV_df['c_elyte9_00'], label = elyte_species[0, 8])

ax1.plot(t, SV_df['c_SEI2_00'], '-b', label = SEI_species[0, 1])
#ax1.plot(t, SV_df['c_SEI3_00'], '-r', label = SEI_species[0, 2])
ax1.plot(t, SV_df['c_SEI4_00'], '-g', label = SEI_species[0, 3])
ax1.set_ylabel('Concentration ' r"$[\frac{kmol}{m^3}]$")
ax1.set_title('SEI species concentrations over time')
ax1.legend(loc='upper left')

ax2.plot(t_vec, phi_vec, label = 'Anode Voltage')
ax2.set_ylabel('Anode Potential [V]')
ax2.set_xlabel('Time [s]')
ax2.set_title('Anode voltage over time')

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