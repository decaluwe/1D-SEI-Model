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
import numpy as np

"Identify cti file:"
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
N_x = 1         # USER INPUT number of grids in plane of electrode
N_y = 5        # USER INPUT number of grids perpendicular to electrode
x = 1           # USER INPUT x length of domain [m]
y = 5e-8        # USER INPUT y length of domain [m]
d_sei = 5e-9    # USER INPUT d_SEI representative diameter of SEI grain [m]
# %% Define state variables

"""----------Initial state variables----------"""
T_0 = np.array([300.])      # USER INPUT initial temperature [K], assumed uniform
P_0 = 101325.               # USER INPUT defines initial pressure.
eps_0 = np.array([0.])      # initial volume fraction of SEI

"""
Define CV parameters:
"""
sweep_rate = 0.01  #...Voltage sweep rate [V/s]
sweep_dirn_0 = -1
phi_0 = 1.0;
phi_1 = 0.05;
phi_2 = 1.5;
n_cycles = 0.5

# If you want to verify that the electric potential input looks correct before
#     running the simulation, switch this to '1'
check_profile = 0
#
#...For these half-cell simulations, let us assume a constant cathode
#       voltage for our energy storage and efficiency calcs:
phi_CE = 0.0; #...(V)yte interface

#...Initial electric potential of SEI, relative to WE:
phi_SEI_dl_0 = 0.
#...Initial electric potential of electrolyte:
phi_elyte_0 = 0.

"""
SEI properties:
"""
# Mass density [kg of k per m3 of k] requires knowing order of species in cti
#   file.
rho_k_SEI = [1e-20, 2110, 2013, 1321]
# Electrical Resistivity [S/m]
sigma_el = 1e-6

"""
Double layer capacitances
"""
C_dl_WE_SEI = 2e-4  # F/m2
