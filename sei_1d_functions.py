# -*- coding: utf-8 -*-
"""
Created on Thu May 17 18:49:06 2018

@author: Daniel Korff

This code takes in an empty solution vector for tracking temperature,
concentration, and potential of 2 species based on the grid used in the top
level script and updates the column headers to reflect the data stored in the
respective columns.
"""
import numpy as np
import cantera as ct

"""----------Residual function for IDA solver----------"""
def residual(t, SV, SV_dot):
    from sei_1d_init import objs, params, voltage_lookup, SVptr
    # Tag global variables in residual function
    print(t)
    # Initialize residual equal to all zeros:
    res = SV_dot

    # Read out cantera objects:
    WE = objs['WE']
    sei = objs['SEI']
    elyte = objs['elyte']
    sei_elyte = objs['SEI_elyte']
    sei_conductor = objs['conductor']
    WE_sei = objs['WE_SEI']

    # Set temperature, density, and concentration of species
    #elyte.TPX = params['TP'], SV[]

    phi_WE = np.interp(t,voltage_lookup['time'],voltage_lookup['voltage'])
    #print(phi_WE)

    i_sei = np.zeros(params['Ny']+1,)
    i_Far = np.zeros(params['Ny'],)

    # Set anode and anode/electrolyte interface potentials
    WE.electric_potential = phi_WE

    "TEMPORARY:"

    sei.electric_potential = phi_WE + SV[SVptr['phi sei'][0]]
    sei_conductor.electric_potential = phi_WE + SV[SVptr['phi sei'][0]]

    i_sei[0] = WE_sei.get_net_production_rates(WE)*ct.faraday


    """Ck_sei_loc = SV[SVptr['Ck sei'][j]]
    #print(Ck_sei_loc)
    rho_sei_loc = abs(np.dot(Ck_sei_loc,sei.molecular_weights))
    Xk_sei_loc = Ck_sei_loc/sum(Ck_sei_loc)
    sei.TDX = None, rho_sei_loc, Xk_sei_loc

    elyte.electric_potential = 0.
    Ck_elyte_loc = SV[SVptr['Ck elyte'][j]]

    eps_sei_loc = SV[SVptr['eps sei'][j]]
    # SEI surface Area Per unit Volume (APV)
    A_sei_APV = 4.*eps_sei_loc*(1-eps_sei_loc)/params['d_sei']

    #anode_elyte.electric_potential = phi_anode
    #elyte.electric_potential = 0

    # Retreive net production rates of species from cantera
#    S_dot_elyte = anode_elyte.net_production_rates[elyte_rate_ind]
    #S_dot_SEI = anode_elyte.net_production_rates[SEI_rate_ind]
    #print(sei_elyte.get_net_production_rates(sei))
    dSVdt_ck_sei = sei_elyte.get_net_production_rates(sei)*A_sei_APV
    res[SVptr['Ck sei'][j]] = SV_dot[SVptr['Ck sei'][j]] - dSVdt_ck_sei

    dSVdt_eps_sei = np.dot(dSVdt_ck_sei, params['vol_k sei'])
    res[SVptr['eps sei'][j]] = SV_dot[SVptr['eps sei'][j]] - dSVdt_eps_sei"""
    rxn_scale = 1
    for j in range(params['Ny']):

        Ck_sei_loc = SV[SVptr['Ck sei'][j]]
        #print(Ck_sei_loc)
        rho_sei_loc = abs(np.dot(Ck_sei_loc,sei.molecular_weights))
        Xk_sei_loc = Ck_sei_loc/sum(Ck_sei_loc)
        sei.TDX = None, rho_sei_loc, Xk_sei_loc

        elyte.electric_potential = 0.
        Ck_elyte_loc = SV[SVptr['Ck elyte'][j]]

        eps_sei_loc = SV[SVptr['eps sei'][j]]
        # SEI surface Area Per unit Volume (APV)
        A_sei_APV = 4.*eps_sei_loc*(1-eps_sei_loc)/params['d_sei']

        #anode_elyte.electric_potential = phi_anode
        #elyte.electric_potential = 0

        # Retreive net production rates of species from cantera
    #    S_dot_elyte = anode_elyte.net_production_rates[elyte_rate_ind]
        #S_dot_SEI = anode_elyte.net_production_rates[SEI_rate_ind]
        #print(sei_elyte.get_net_production_rates(sei))
        Rates_sei_elyte = sei_elyte.get_net_production_rates(sei)*A_sei_APV
        dSVdt_ck_sei = rxn_scale*Rates_sei_elyte
        res[SVptr['Ck sei'][j]] = SV_dot[SVptr['Ck sei'][j]] - dSVdt_ck_sei

        dSVdt_eps_sei = np.dot(dSVdt_ck_sei, params['vol_k sei'])
        res[SVptr['eps sei'][j]] = SV_dot[SVptr['eps sei'][j]] - dSVdt_eps_sei

        # Calculate faradaic current density due to charge transfer at SEI-elyte
        #   interface, in A/m2 total.
        Rates_sei_conductor = sei_elyte.get_net_production_rates(sei_conductor)
        i_Far[j] = Rates_sei_conductor*A_sei_APV/params['dyInv']

        # Rates for next node are scaled by sei volume fraction in this node:
        rxn_scale = eps_sei_loc

    i_dl = i_sei[0] - i_sei[1] - i_Far[0]
    dSVdt_phi_dl = -i_dl/params['C_dl WE_sei']
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













# %% The function below will apply column labels to the data frame sol_vec
def df_2spec2var(sol_vec,N_x,N_y,len_sol_vec,track_vars,track_temp,num_species):
    # Initialize base strings

    temp_tag        = ["T_" for i in range(N_x*N_y)]
    V_elyte_tag     = ["V_elyte_" for i in range(N_x*N_y)]
    V_sei_an_tag    = ["V_sei_an_" for i in range(N_x*N_y)]
    V_sei_elyte_tag = ["V_sei_elyte_" for i in range(N_x*N_y)]
    phi_tag         = ["phi_" for i in range(N_x*N_y)]
    c_elyte_tag = {}
    c_SEI_tag = {}
    for i in range(num_species[0]):
        c_elyte_tag[str(i+1)] = ["c_elyte"+str(i+1)+"_" for j in range(N_x*N_y)]
    for i in range(num_species[1]):
        c_SEI_tag[str(i+1)] = ["c_SEI"+str(i+1)+"_" for j in range(N_x*N_y)]

    # concatenate appropriate indices onto base strings

    for j in range(N_y):
        for i in range(N_x):
            temp_tag[i+j*N_x]        = temp_tag[i+j*N_x]+str(j)
            V_elyte_tag[i+j*N_x]     = V_elyte_tag[i+j*N_x]+str(j)
            V_sei_an_tag[i+j*N_x]    = V_sei_an_tag[i+j*N_x]+str(j)
            V_sei_elyte_tag[i+j*N_x] = V_sei_elyte_tag[i+j*N_x]+str(j)
            for key in c_elyte_tag:
                c_elyte_tag[key][i+j*N_x] = c_elyte_tag[key][i+j*N_x]+str(j)
                c_elyte_tag[key][i+j*N_x] = c_elyte_tag[key][i+j*N_x]+str(i)
            for key in c_SEI_tag:
                c_SEI_tag[key][i+j*N_x] = c_SEI_tag[key][i+j*N_x]+str(j)
                c_SEI_tag[key][i+j*N_x] = c_SEI_tag[key][i+j*N_x]+str(i)
            phi_tag[i+j*N_x]         = phi_tag[i+j*N_x]+str(j)

            temp_tag[i+j*N_x]        = temp_tag[i+j*N_x]+str(i)
            V_elyte_tag[i+j*N_x]     = V_elyte_tag[i+j*N_x]+str(i)
            V_sei_an_tag[i+j*N_x]    = V_sei_an_tag[i+j*N_x]+str(i)
            V_sei_elyte_tag[i+j*N_x] = V_sei_elyte_tag[i+j*N_x]+str(i)
            phi_tag[i+j*N_x]         = phi_tag[i+j*N_x]+str(i)

    # Determine columns where each variable is being stored

    V_ctr = track_temp + 0
    elyte_species_ctr = track_temp + 3
    SEI_species_ctr = elyte_species_ctr + num_species[0]
    phi_ctr = SEI_species_ctr + num_species[1]
    c_elyte_cols = [0 for i in range(num_species[0])]
    c_SEI_cols = [0 for i in range(num_species[1])]

    temp_cols = list(range(0, len_sol_vec, track_vars))
    V_elyte_cols = list(range(V_ctr, len_sol_vec, track_vars))
    V_sei_an_cols = list(range(V_ctr+1, len_sol_vec, track_vars))
    V_sei_elyte_cols = list(range(V_ctr+2, len_sol_vec, track_vars))
    for i in range(num_species[0]):
        c_elyte_cols[i] = list(range(elyte_species_ctr+i, len_sol_vec, track_vars))
    for i in range(num_species[1]):
        c_SEI_cols[i] = list(range(SEI_species_ctr+i, len_sol_vec, track_vars))
    phi_cols = list(range(phi_ctr, len_sol_vec, track_vars))

    # Loop over all columns updating column names to appropriate var label

    count = 0
    for i in range(N_x*N_y):
        if track_temp == 1:
            sol_vec = sol_vec.rename(columns = \
                                     {temp_cols[i]:temp_tag[count]})
        sol_vec = sol_vec.rename(columns = \
                                 {V_elyte_cols[i]:V_elyte_tag[count]})
        sol_vec = sol_vec.rename(columns = \
                                 {V_sei_an_cols[i]:V_sei_an_tag[count]})
        sol_vec = sol_vec.rename(columns = \
                                 {V_sei_elyte_cols[i]:V_sei_elyte_tag[count]})
        j = 0
        for key in c_elyte_tag:
            sol_vec = sol_vec.rename(columns = \
                                     {c_elyte_cols[j][i]:c_elyte_tag[key][count]})
            j += 1
        j = 0
        for key in c_SEI_tag:
            sol_vec = sol_vec.rename(columns = \
                                     {c_SEI_cols[j][i]:c_SEI_tag[key][count]})
            j += 1
        sol_vec = sol_vec.rename(columns = \
                                 {phi_cols[i]:phi_tag[count]})
        count += 1

    # Return cleaned up solution vector to top level script

    return sol_vec

# %% The following function applies boundary conditions to sol_vec


#-----------------------------------------------------------------------------
