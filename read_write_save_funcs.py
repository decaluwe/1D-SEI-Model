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



""" Save File Copies """
"-----------------------------------------------------------------------------"
def SaveFiles(folder_name, ctifile,sv_save, names):
    import os, sys
    import numpy as np
    import cantera as ct
    from shutil import copy2, rmtree
    
    """ Set up saving location """
    "-------------------------------------------------------------------------"
    cwd = os.getcwd()
    
    # Create folder for any files/outputs to be saved:

    os.makedirs(folder_name)
    copy2(cwd + '/../' + 'sei_1d_inputs.py', folder_name)
    copy2(cwd + '/../' + 'sei_1d_functions.py', folder_name)
    copy2(cwd + '/../' + 'sei_1d_init.py', folder_name)
    copy2(cwd + '/../' + 'sei_1d_model.py', folder_name)
    #ModuleWriter(cwd + '/' + folder_name + '/user_inputs.csv', user_inputs)
    
    # Save the current cti files into new folder:
    cti_path = ct.__path__[0]
    if os.path.exists(cwd + '/../' + ctifile):
        copy2(cwd + '/../' + ctifile, folder_name)
    else:
        copy2(cti_path + '/data/' + ctifile, folder_name)
        
    # Save the current parameters for post processing:
    #DictWriter(cwd + '/' + folder_name + '/params.csv', p)
    
    # Save a copy of the full solution matrix:
    np.savetxt(folder_name +'/solution.csv', sv_save, delimiter=',')
    np.savetxt(folder_name+'/names.csv',names,delimiter=",", fmt="%s")
    