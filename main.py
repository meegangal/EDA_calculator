import os

''' This script will be used to run energy decomposition analysis (EDA) for a transition state using the gibb's free energy obtained 
from a single point energy calculation from the TS structure. 
EDA requires 6 files for a bimolecular TS (every additional molecule involved in TS require 2 files):
(1) The TS file, (2) molecule 1 optimized structure, (3) molecule 1 TS fragment structure, 
(4) molecule 2 optimized structure, (5) molecule 2 TS fragment structure, (6) fort.99 file containing electrostatic potential information
All structure files will have '.log' extension'''

    ## each TS will have 6 files required for EDA
        # TS -- this will have the highest number of atoms (Natoms)!
        # cat frag  -- this will have x number of atoms
        # cat opt   -- this will have x numbet of atoms and be lower in energy than frag
        # rct frag  -- this will have Natoms - x number of atoms
        # rct opt   -- this will have Natoms - x number of atoms be lower in energy than frag 
        # esp   -- this will be a text file, that will be originally a .99 file

#py installer to convert to executable

def getData(directory):
    ''' only '.log','.99' or'.txt' files will be accepted
    '.99' files will be converted to '.txt' files
    Default will be current directory but can input a specific directory as an argument
    This will omit all files that do not have '.log', or '.txt' extensions and parse required data for EDA.
    Required data will be NAtoms from all files; single point energy and dispersion energy from each '.log' file
    INPUT: Directory path containing all required files for EDA
    OUTPUT: Dictionary with each file name as keys. Values will be a list containing number of atoms, single point energy and dispersion energy
    it will be organized in a way in which each even index will be the label of the numerical value and each odd index will be the numerical value
    '''
    #path = Path(directory)
    #os.chdir(path)
    #working_directory = os.getcwd()
    files = os.listdir(directory) #why wont this work!!!
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
                if 'Normal' in line: #this truncates the file to only the last calculation in the .log file
                    indices.append(index) # this could be unnecessary 
                    Normal.append(line) # so if there are multiple calculations this will take the single poiint of the last single point
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
    EDA['dE'] = (TS[0] - sum(optimized))*627.51
    # this is to get distortion energy 
    for key, value in molecules.items():
        distortion = (value[1] - value[0])*627.51
        dist.append(distortion)
        EDA[key + ' distortion'] = distortion
    EDA['Total Distortion'] = sum(dist)
    #this is to get interaction energy
    interaction = (TS[0] - sum(fragments))*627.51
    EDA['Interaction Energy'] = interaction
    #this is to get dispersion energy
    dispersion = (TS[1] - sum(frag_dispersion))*627.51
    EDA['Dispersion Energy'] = dispersion
    #this is to get ESP 
    espfile = espinfo[0]
    chargefile = []
    for file in frag_files:
        if data[file[0]][1] == espinfo[2]:
            chargefile.append(file[0])
    chargefile.append(espinfo[2])
    dESP = calculateESP(chargefile, espfile)*627.51
    EDA['dESP'] = dESP
    ERCT = interaction - (dispersion + dESP)
    EDA['ERCT'] = ERCT

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

def getDifference(path1, path2):
    ''' will calculate the difference between two EDA
    INPUT: directory paths of both EDA files
    OUTPUT: dictionary with the difference between each value in both EDAs
    note: EDA_relative = EDA_minor - EDA_major'''
    EDA1 = getEDA(path1)
    EDA2 = getEDA(path2)
    if len(EDA1) == len(EDA2):
        difference = {key: EDA1[key] - EDA2.get(key, 0) for key in EDA1}
    else:
        print('ERROR: Mismatched number of values')
        difference = {}
    return difference
    
def makeDatasheet():
    #pandas??
    pass




#spath = os.path.join("/Users","meegangalante","Downloads","python","EDA_calculator","S_TS")
#rpath = os.path.join("/Users","meegangalante","Downloads","python","EDA_calculator","R_TS")

#R_EDA = getEDA(rpath)
#S_EDA = getEDA(spath)
#print(getDifference(rpath,spath))
