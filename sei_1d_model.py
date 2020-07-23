# -*- coding: utf-8 -*-
"""
Created on Wed May 16 14:52:05 2018
@author: Steven C. DeCaluwe, Daniel Korff

This code is the foundation (top level script) for a model of solid-electrolye-
interface (SEI) formation in a Li-ion battery (LIB). The model is set to track
the electric potential, cheical composition, and volume fraction of an SEI as
it grows in 1 dimension through a series of discretized finite volumes.

While the model's solution vector also stores the electrolyte composition, at
present (July 24, 2019), there is no functionality in the model which considers
eveolution of the electrolyte state variables throughout the simulation domain.
"""

#  Import all necessary modules
import numpy as np
import cantera as ct
from assimulo.solvers import IDA
from assimulo.problem import Implicit_Problem
import os

print('\n     Importing inputs and intializing.')

# This function imports all user inputs, creates all needed
#    Cantera objects, dictionaries for parameters, pointers, etc., and
#    initializes the solution vector.
from sei_1d_init import SV_0, SV_dot_0, SVptr, times, objs, params,  \
    save_name, ctifile, mode

if mode == 'detailed':
    from sei_1d_functions import residual_detailed as residual
    from sei_1d_outputs_detailed import *
elif mode == 'homogeneous':
    from sei_1d_functions import residual_homogeneous as residual
    from sei_1d_outputs_homogeneous import *
elif mode == 'reduced':
    from sei_1d_functions import residual_reduced as residual
    from sei_1d_outputs_reduced import *
else:
    raise Exception('Please select a valid simulation mode." + \
        "  Available options are "detailed", "homogeneous", and "reduced"')

print('\n     Running simulation\n')
# Set up problem instance
SEI_1D = Implicit_Problem(residual, SV_0, SV_dot_0)

# Define simulation parameters
simulation = IDA(SEI_1D)                # Create simulation instance
simulation.atol = 1e-8                  # Solver absolute tolerance
simulation.rtol = 1e-6                  # Solver relative tolerance
#simulation.maxh = 55                   # Solver max step size

# Set simulation end time, slope flag (for anode voltage cycle), and run simulation

t_f = times[-1]

#ncp_list = np.arange(0, t_f, 0.15)
#ncp = 10000

" Run the simulation "
"----------------------------------------------------------------------------------"
t, SV, SV_dot = simulation.simulate(t_f)

" Organize, plot, and save the data:"
"----------------------------------------------------------------------------------"
# This function creates a dataframe and adds labels to each element in the SV:
# CURRENTLY NOT UP TO DATE - SD, 7/18/19
#   SV_df = output_names(SV,N_x,N_y,len_sol_vec,track_vars,track_temp,num_species)

# Put t and SV into a format for saving, and create array of SV variable names
#  TEMPORARY - EVENTUALLY REPLACE WITH output_names DATAFRAME APPROACH.
data, data_names = prepare_data(SV, t, objs, params)

# Save the data:
save_files(save_name, ctifile, data.T, data_names)

# Plot the data
plot_data(t, SV, SVptr, objs, params)

print("\n     Goodbye")

# Functions !!!!!NOT CURRENTLY USED!!!!!
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
