"""
The following functions were created to read/write variable values from/to .csv
files. Use of these functions allows simple saving and reading of these
variables regardless of their storage method.

A SaveFiles function was also added to easily create copies of files used to
run the model. This allows the user to go back and check how the solution was
calculated at that time even if the current version of the model has been
updated to fix bugs or incorporate additional physics.
"""



""" Read and Write w.r.t. Modules """
"-----------------------------------------------------------------------------"
def ModuleWriter(file, module):
    import types, csv

    f = open(file, 'w')
    w = csv.writer(f, lineterminator='\n')

    for item in dir(module):
        if not item.startswith("__"):
            if type(vars(module)[item]) != types.ModuleType:
                w.writerow([item, vars(module)[item]])

    f.close()

def ModuleReader(file):
    import csv, numpy

    f = open(file, 'r')
    reader = csv.reader(f)

    d = {}
    for row in reader:
        k, v = row

        if '[' not in v:
            try:
                d[k] = eval(v)
            except:
                d[k] = v
        else:
            d[k] = " ".join(v.split()).replace(' ',', ')
            d[k] = numpy.asarray(eval(d[k]))

    f.close()

    return d



""" Read and Write w.r.t. Dictionaries """
"-----------------------------------------------------------------------------"
def DictWriter(file, dictionary):
    import csv

    f = open(file, 'w')
    w = csv.writer(f, lineterminator='\n')

    for k,v in dictionary.items():
        w.writerow([k,v])

    f.close()

def DictReader(file):
    import csv, numpy

    f = open(file, 'r')
    reader = csv.reader(f)

    p = {}
    for row in reader:
        k, v = row

        if '[' not in v:
            p[k] = eval(v)
        else:
            p[k] = " ".join(v.split()).replace(' ',', ')
            p[k] = numpy.asarray(eval(p[k]))

    f.close()

    return p

""" Create output arrays for the data and the variable names """
"-----------------------------------------------------------------------------"
def prepare_data(SV, t, objs, params):
    import numpy as np

    sei = objs['SEI']
    elyte = objs['elyte']

    t = np.asarray(t)
    t.shape = (t.shape[0],1)
    data = np.concatenate((np.array(t),SV),1)

    SVnames = list()
    #SVnames.append('time')
    SVnames.append('phi SEI')
    SVnames.append('phi elyte')
    SVnames.append('eps SEI')

    for i in range(1,sei.n_species):
        SVnames.append(sei.species_names[i])

    for i in range(1,elyte.n_species):
        SVnames.append(elyte.species_names[i])


    data_names = np.tile(SVnames, params['Ny'])

    return data, data_names


""" Save File Copies """
"-----------------------------------------------------------------------------"
def save_files(save_name, ctifile, data, names):
    import os, sys
    import numpy as np
    import cantera as ct
    from shutil import copy2, rmtree
    import datetime

    """ Set up saving location """
    "-------------------------------------------------------------------------"
     # Read out current working directory:
    cwd = os.getcwd()
    try:
        os.chdir(cwd + '/output')
    except:
        os.mkdir(cwd + '/output')
        os.chdir(cwd + '/output')

    folder_name = os.getcwd()+'/'+save_name+datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S")

     # Create folder for any files/outputs to be saved:
    os.makedirs(folder_name)


    """ Copy and save files to this location """
    "-------------------------------------------------------------------------"
    copy2(cwd + '/sei_1d_inputs.py', folder_name)
    copy2(cwd + '/sei_1d_functions.py', folder_name)
    copy2(cwd + '/sei_1d_init.py', folder_name)
    copy2(cwd + '/sei_1d_model.py', folder_name)
    #ModuleWriter(cwd + '/' + folder_name + '/user_inputs.csv', user_inputs)

    # Save the current cti files into new folder:
    cti_path = ct.__path__[0]
    if os.path.exists(cwd + '/' + ctifile):
        copy2(cwd + '/' + ctifile, folder_name)
    else:
        copy2(cti_path + '/data/' + ctifile, folder_name)

    # Save the current parameters for post processing:
    #DictWriter(cwd + '/' + folder_name + '/params.csv', p)

    # Save a copy of the full solution matrix:
    np.savetxt(folder_name +'/solution.csv', data, delimiter=',')
    np.savetxt(folder_name+'/names.csv',names,delimiter=",", fmt="%s")


# The function below will apply column labels to the data frame sol_vec
def output_names(sol_vec,N_x,N_y,len_sol_vec,track_vars,track_temp,num_species):
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

def plot_data(t, SV, SVptr, objs, params):
    from sei_1d_init import voltage_lookup
    from matplotlib import pyplot as plt
    import numpy as np

    sei = objs['SEI']
    elyte = objs['elyte']
    names = list()
    for i in range(sei.n_species):
        names.append(sei.species_names[i])

    names.append('Anode potential')
    for i in range(params['Ny']):
        names.append('SEI potential_'+str(i))

    phi_WE = np.interp(t,voltage_lookup['time'],voltage_lookup['voltage'])


    C_k = SV[-1,SVptr['Ck sei']]
    eps_k_sei = np.zeros_like(C_k)
    vol_k = C_k*sei.partial_molar_volumes
    v_tot = np.dot(C_k,sei.partial_molar_volumes)
    eps_k_sei = vol_k/v_tot


    names = list()
    for i in range(sei.n_species):
        names.append(sei.species_names[i])

    """fig5, ax5 = plt.subplots(1, 1, figsize=(8, 7.2))
    ax5.plot(1e9*np.arange(params['Ny'])/params['dyInv'],eps_k_sei)
    ax5.plot(1e9*np.arange(params['Ny'])/params['dyInv'],1.-SV[-1,SVptr['eps sei']])
    ax5.legend(names)
    ax5.set_ylabel('Species volume fraction')
    ax5.set_xlabel('SEI Depth (from anode, nm)')
    plt.show()"""
    """plt.savefig('Figure2.pdf',format='pdf',dpi=350)"""
