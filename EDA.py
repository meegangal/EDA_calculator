#!/usr/bin/python
from __future__ import print_function, absolute_import
# ^ use this if using on remote computer via ssh

###########################################################################################################
#                                           EDA.py                                                        #  
# Python script written by Meegan L. Galante to calculate energy decompositions                           #
# of a transition state from a Gaussian transition state search. Inorder to calculate                     #
# EDA using this script, a directory must be created containing only required EDA files.                  #
#                                                                                                         #
# The EDA files required are single point energy calculations for the TS,                                 #  
# for fragment structures from the TS structure for each molecule involved,                               #  
# for optimized structures for each molecule and the fort.99 output                                       #  
# file from an electrostatic potential calculation. '.log' files with multiple linked jobs                #
# can be used as long as the final calculation is the single point energy calculation                     #                                  
#                                                                                                         #
# #########################################################################################################
#                                                                                                         #  
# Run this script from the command line by writing:                                                       #  
#    $EDA.py getEDA "/Users/directory/of/your/files1"                                                     #  
# This will give you EDA values for the TS in that directory                                              #  
# If no directory is given, will give you EDA of current working directory                                #  
#     $EDA.py finalEDA "/Users/directory/of/your/files1" "/Users/directory/of/your/files2"                #
# This will give you EDA values for each TS and also give you the difference between them                 #  
# If no directory is given, will give you EDA and difference of directories in current working directory  #  
#       - best way to use is by making a directory for the EDA that contains                              #  
#         directories for the two TS' you are comparing and running the 'EDA.py finalEDA' command         #
#       - if using on local computer 'python' required prior to command                                   #
#       - if using on local computer, must be in same directory as script to run                          #    
#                                                                                                         #  
# Last modified April 23 2025                                                                             #  
###########################################################################################################

import os, sys

def getData(directory):
    ''' only '.log','.99' or'.txt' files will be accepted
    '.99' files will be converted to '.txt' files '''

    files = os.listdir(directory)
    working_files = []
    data = {}
    os.chdir(directory)
    for file in files:
        if file[-2:] == '99':
            new_file_name = file[:-2] + 'txt'
            os.rename(file, new_file_name)
            working_files.append(new_file_name)
        elif file[-3:] == 'log' or file[-3:] == 'txt':
            working_files.append(file)
    for file in working_files:
        if file[-3:] == 'log':
            with open(file) as f:
                lines = f.readlines()
            indices = []
            Normal = []
            for index, line in enumerate(lines):
                if 'Normal' in line:
                    indices.append(index)
                    Normal.append(line) 
            if len(indices) > 1:
                lines = lines[indices[-2]:]
            for line in lines:
                if 'NAtoms='in line:
                    row = line.split()
                    NAtoms = int(row[1])
            data.setdefault(file, [])
            data[file].extend(['NAtoms', NAtoms])
            sp, disp = getEnergy(file)
            data[file].extend(['Single Point Energy', sp, 'Dispersion', disp])
        if file[-3:] =='txt':
            with open(file) as f:
                lines = f.readlines()
            NAtoms = int(len(lines)/2)
            data.setdefault(file, [])
            data[file].extend(['NAtoms', NAtoms])
    return data

def getEnergy(filename):
    ''' this will parse log file for single point energy value and dispersion energy
        INPUT: File name
        OUTPUT: Tuple containing single point energy and disperion energy'''
    singlepoint = 'ERROR'
    if filename[-3:] == 'log':
        with open(filename) as f:
            lines = f.readlines()
        for line in lines:
            row = line.split()
            if 'Done:' in row:
                singlepoint = float(row[4])
            if 'Dispersion' in row:
                dispersion = float(row[4])
    
    return singlepoint, dispersion

def getEDA(current_dir=os.getcwd()):
    '''this will go through the data dict and calculate distortion and interaction energy
    INPUT: directory path, data dictionary, final_esp 
    OUTPUT: Dictionary called EDA that will contain dE, distortion for molecule 1,
    distortion for molecule 2, total distortion, interaction energy, dispersion, ESP and ERCT'''
    
    data = getData(current_dir)
    espinfo = []
    NAtoms = []
    for key in data:
        if len(data[key]) < 3:
            espinfo.append(key)
            espinfo.extend(data[key])
        NAtoms.append(data[key][1])
    del data[espinfo[0]]
    NAtoms = list(set(NAtoms))
    TS_n = max(NAtoms)
    del NAtoms[NAtoms.index(TS_n)]
    molecules = {}
    for i in range(len(NAtoms)):
        molecules.setdefault('mol{0}'.format(i+1), [])
        molecules['mol{0}'.format(i+1)].append(NAtoms[i])
    TS = []
    TS_name = ""
    #TS is list with first element as single point and second element as dispersion
    frag_dispersion = []
    if TS_n == sum(NAtoms):
       for key in data.keys():
            for natom in NAtoms:
                if natom in data[key]:
                    mol_key = [key for key, value in molecules.items() if value[0] == natom]
                    molecules[mol_key[0]].append(data[key][3])
            if TS_n in data[key]:
                TS.append(data[key][3])
                TS.append(data[key][5])
                TS_name = key
    else:
        print('Number of atoms do not add up to TS structure!')
    for key, value in molecules.items():
        del molecules[key][0]
        molecules[key] = sorted(molecules[key]) 
    frag_files = []
    for mol in molecules.keys():
        file = [key for key, value in data.items() if molecules[mol][1] in value]
        frag_dispersion.append(data[file[0]][5])
        frag_files.append(file)
    #this sorts so that first term is the optimized struct and second is fragment
    #distortion will be frag - opt, or index1 - index0
    EDA = {}
    dist = []
    optimized = []
    fragments = []
    # this is to get dE
    for key, value in molecules.items():
        optimized.append(value[0])
        fragments.append(value[1])
    EDA['dE:                '] = (TS[0] - sum(optimized))*627.5095
    # this is to get distortion energy 
    for key, value in molecules.items():
        distortion = (value[1] - value[0])*627.5095
        dist.append(distortion)
        EDA[key + ' distortion:   '] = distortion
    EDA['Total Distortion:  '] = sum(dist)
    #this is to get interaction energy
    interaction = (TS[0] - sum(fragments))*627.5095
    EDA['Interaction Energy:'] = interaction
    #this is to get dispersion energy
    dispersion = (TS[1] - sum(frag_dispersion))*627.5095
    EDA['Dispersion Energy: '] = dispersion
    #this is to get ESP 
    espfile = espinfo[0]
    chargefile = []
    for file in frag_files:
        if data[file[0]][1] == espinfo[2]:
            chargefile.append(file[0])
    chargefile.append(espinfo[2])
    dESP = calculateESP(chargefile, espfile)*627.5095
    EDA['dESP:              '] = dESP
    ERCT = interaction - (dispersion + dESP)
    EDA['ERCT:              '] = ERCT
    print(f'------------------------------------------------------\nEDA of {TS_name}')
    for key, value in EDA.items():
        print(f'{key}   {value}')
    return EDA

def OrganizeCharge(fname):
    ''' This will parse the charge information for each atom in the reactant or whichever molecule matches the fort.99 file
    INPUT: file name (fragment structure)
    OUTPUT: list of charges, len(charges) should be same length of lines in fort.99 file'''
    atom_type = []
    charge = []
    with open(fname[0]) as f:
        lines = f.readlines()
    index = []
    for line in lines:
        row = line.split()
        if 'Summary' in row and 'Natural' in row:
            index = lines.index(line)
    lines = lines[index+6:index+fname[1]+6]
    for line in lines:
        row = line.split()
        atom_type.append(row[0])
        charge.append(float(row[2]))
    return atom_type, charge

def OrganizeESP(fname):
    ''' This gets the esp values from the fort.99
    INPUT: file name (fort.txt)
    OUTPUT: list containing all esp values, len(charges) == len(esp)'''
    esp = []
    with open(fname) as f:
        lines = f.readlines()
        for line in lines:
            if line[0] != '}':
                row = line.split()
                if len(row) == 4:
                    esp.append(float(row[-1].rstrip('\\')))    
    return esp

def calculateESP(chargefile, espfile):
    ''' Calculates the energy contribution of electrostatic potential to the interaction energy
    INPUT: charges and esp
    OUTPUT: float called final_esp'''
    global atoms
    atoms, charges = OrganizeCharge(chargefile)
    esp_values = OrganizeESP(espfile)
    if len(charges) == len(esp_values):
        product = [x * y for x, y in zip(charges,esp_values)]
    final_esp = sum(product)
    return final_esp

def finalEDA(path_1, path_2):
    ''' will calculate the difference between two EDA
    INPUT: directory paths of both EDA files
    OUTPUT: dictionary with the difference between each value in both EDAs
    note: EDA_relative = EDA_minor - EDA_major'''
    EDA1 = getEDA(path_1)
    EDA2 = getEDA(path_2)
    if len(EDA1) == len(EDA2):
        difference = {key: EDA1[key] - EDA2.get(key, 0) for key in EDA1}
    else:
        print('ERROR: Mismatched number of values')
        difference = {}
    print('------------------------------------------------------\nDifference:')
    for key, value in difference.items():
        print(f'{key}   {value}')
    print('------------------------------------------------------')
    return difference 


#this will make it so you can run script from the command line
if __name__ == '__main__':
    if len(sys.argv) > 2:
        args = sys.argv[2:]
        globals()[sys.argv[1]](*args)
    elif sys.argv[1] == 'finalEDA':
        args = [os.getcwd() + '/' + file for file in os.listdir(os.getcwd())]
        globals()[sys.argv[1]](*args)
    else:
        globals()[sys.argv[1]]()


