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
    exit("[Error] Seaborni not installed and is required.")
    
try:
    import pandas as pd
except ImportError:
    exit("[Error] pandas not installed and is required.")   

"""
try:
    #Modify the pymol to launch quietly
    import __main__
    __main__.pymol_argv = [ 'pymol', '-q' ]
    
    #Import PyMOL
    import pymol
    
    #Load the module
    pymol.finish_launching()

    #Load the cmd module
    from pymol import cmd
except ImportError:
    exit("[Error] Pymol is not installed and is required.")
"""

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
Function: PACT Import and Pre-processing
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
Function: Pre-process dataset
"""
# strip out WT, UNK, STOPs, and low read counts
mask_fmc = (
        (df_pact['FMC63_Epitope_Drop_classified'] != 'WT') &      
        (df_pact['FMC63_Epitope_Drop_classified'] != 'STOP') &  
        (df_pact['FMC63_Epitope_Drop_ref_counts'] >= 12)
        )

mask_4g7 = (
        (df_pact['4G7_Epitope_Drop_classified'] != 'WT') &      
        (df_pact['4G7_Epitope_Drop_classified'] != 'STOP') &
        (df_pact['4G7_Epitope_Drop_ref_counts'] >= 12)
        )

# get the cleaned up drop dataframe
df_fmc_mask = df_pact[mask_fmc]
df_4g7_mask = df_pact[mask_4g7]

# return the searched locations
np_fmc_loc = np.unique(df_fmc_mask['Location'].astype(int))
np_4g7_loc = np.unique(df_4g7_mask['Location'].astype(int))


# create lists to work into
list_fmc = []
list_4g7 = []
list_color = []
list_rasa = []

# bin each site and calculate the mean z-score
for loc in np_4g7_loc:
    
    # return a df at the location
    df_loc_f = df_fmc_mask[df_fmc_mask['Location'] == loc]
    df_loc_4 = df_4g7_mask[df_4g7_mask['Location'] == loc]
    
    # calculate the mean z-score
    mean_fmc = np.mean(df_loc_f['FMC63_Epitope_Drop_sd_from_wt'])
    mean_4g7 = np.mean(df_loc_4['4G7_Epitope_Drop_sd_from_wt'])
    
    # append the mean
    list_fmc.append(mean_fmc)
    list_4g7.append(mean_4g7)
           
    rasa = 1 - df_loc_f['frac_burial'].values[0]
    list_rasa.append(rasa)
    
    if rasa > 0.15:
        list_color.append('blue')
    else:
        list_color.append('red')

# plot the estimate
plt.figure(figsize=(8,8))
plt.scatter(list_fmc, list_4g7, s=4, color=list_color)
plt.ylabel('4G7')
plt.xlabel('FMC')
plt.show()



np_fmc = np.array(list_fmc)
np_4g7 = np.array(list_4g7)
np_rasa = np.array(list_rasa)