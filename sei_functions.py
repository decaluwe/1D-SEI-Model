# -*- coding: utf-8 -*-
"""
Created on Thu May 17 18:49:06 2018

@author: Daniel Korff

This code takes in an empty solution vector for tracking temperature,
concentration, and potential of 2 species based on the grid used in the top
level script and updates the column headers to reflect the data stored in the
respective columns.
"""
def CV_voltage():

    sweep_rate = 0.01  #...Voltage sweep rate [V/s]
    sweep_dirn_0 = -1
    phi_0 = 1.0;
    phi_1 = 0.05;
    phi_2 = 1.5;
    n_cycles = 1.


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
