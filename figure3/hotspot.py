#!/usr/bin/env python3
#
# Protein Analysis and Classifier Toolkit
# Author: Justin R. Klesmith
# Copyright (C) 2018-2019 by Justin R. Klesmith
#
# This software is released under GNU General Public License 3
# Additional license options available at http://license.umn.edu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""epitope analysis - load the epitope features and work with them"""

from sys import version_info, exit

# python version check
if version_info < (3,4):
    exit("[PACT Error] Your Python interpreter is too old.\n"
          "Minimium version required is Python 3.4 (64 bit recommended).") 
    
# library import
try:
    import numpy as np
except ImportError:
    exit("[Error] Numpy is not installed and is required.")

try:
    import scipy as sp
except ImportError:
    exit("[Error] SciPy not installed and is required.")

try:
    import matplotlib.pyplot as plt
except ImportError:
    exit("[Error] matplotlib not installed and is required.")
    
try:
    import seaborn as sns
    sns.set(rc={'figure.figsize':(6,6),
                'axes.facecolor':'white',
                'figure.facecolor':'white'
                })
except ImportError:
    exit("[Error] Seaborn not installed and is required.")
    
try:
    import pandas as pd
except ImportError:
    exit("[Error] pandas not installed and is required.")   

# python std lib imports 
from os.path import isfile, dirname, realpath
from sys import argv

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2019, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2019.3"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

"""
Function: PACT Import
"""
# set our path
pact_path = dirname(realpath(argv[0]))

# set our filename
import_features = 'CD19_epitope_features_dataset.csv'

# read in our pact file    
if isfile(pact_path + '/' + import_features):
        
    # load the csv of features into pandas
    df_pact = pd.read_csv(pact_path + '/' + import_features) 

"""
Function: PDB Import
"""
# PDB
import_pdb = '6al5_0001f_0001_0001.pdb'

# read in our pact file    
if isfile(pact_path + '/' + import_pdb):
    
    # open the file
    with open(pact_path + '/' + import_pdb, 'r') as infile:
        list_lines = infile.read().splitlines()

"""
Function: Parse PDB
"""

#Create the dicts for the atom and hetatms
dict_atom = {}

# iterate through each line
for line in list_lines:

    # get the chain
    chain = line[21:22].strip()

    # check if the chain is in the dict
    if chain not in dict_atom:
        dict_atom[chain] = {}

    # get the residue number
    resi_num = int(line[22:26].strip())

    # add the residue to the dict
    if resi_num not in dict_atom[chain]:
        dict_atom[chain][resi_num] = {}
        
    # get the atom name
    atom_name = line[12:16].strip()
    
    # add the atom to the dict
    if atom_name not in dict_atom[chain][resi_num]:
        dict_atom[chain][resi_num][atom_name] = []    

    # add the information to the residue
    dict_atom[chain][resi_num][atom_name].append({
    'atom_name':line[12:16].strip(),
    'res_name':line[17:20].strip(),
    'chain':line[21:22].strip(),
    'resi_num':resi_num,
    'x_coor':float(line[30:38].strip()),
    'y_coor':float(line[38:46].strip()),
    'z_coor':float(line[46:54].strip()),
    'occupancy':float(line[54:60].strip()),
    'temp':float(line[60:66].strip()),
    'element':line[76:78].strip()})
   
#dict_atom[chain][resi_num][atom_name]    
    
"""
Function: Pre-process dataset
"""
# strip out WT, UNK, STOPs, and low read counts
mask_drop = (
        (df_pact['FMC63_Epitope_Drop_classified'] != 'WT') &      
        (df_pact['FMC63_Epitope_Drop_classified'] != 'STOP') &  
        (df_pact['FMC63_Epitope_Drop_ref_counts'] >= 12)
        )

# get the cleaned up drop dataframe
df_masked = df_pact[mask_drop]

# return the searched locations
np_loc = np.unique(df_masked['Location'].astype(int))

# create lists to work into
list_fmc = []
list_4g7 = []


# bin each site and calculate the mean z-score
for loc in np_loc:
    
    # return a df at the location
    df_loc = df_masked[df_masked['Location'] == loc]
    
    # calculate the mean z-score
    mean_fmc = np.mean(df_loc['4G7_Epitope_Drop_sd_from_wt'])
    mean_4g7 = np.mean(df_loc['4G7_20nm_Epitope_Drop_ref_counts'])
    

