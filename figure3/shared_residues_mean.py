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

# return the searched locations
#np_loc = np.unique(df_pact['Location'].astype(int))

np_loc = np.array([48,75,77,82,83,84,85,86,88,89,99,100,
           101,105,136,137,138,140,142,143,144,145,
           146,147,149,196,197,198,199,200,201,202,
           203,204,205,206])

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

mask_3b10 = (
        (df_pact['3B10_Epitope_Drop_classified'] != 'WT') &      
        (df_pact['3B10_Epitope_Drop_classified'] != 'STOP') &
        (df_pact['3B10_Epitope_Drop_ref_counts'] >= 12)
        )

# get the cleaned up drop dataframe
df_fmc = df_pact[mask_fmc]
df_4g7 = df_pact[mask_4g7]
df_3b10 = df_pact[mask_3b10]


# create lists to work into
list_fmc_mean = []
list_4g7_mean = []
list_3b10_mean = []

list_fmc_intol = []
list_4g7_intol = []
list_3b10_intol = []

# bin each site and calculate the mean z-score
for loc in np_loc:
    
    # return a df at the location
    df_loc_f = df_fmc[df_fmc['Location'] == loc]
    df_loc_4 = df_4g7[df_4g7['Location'] == loc]
    df_loc_3 = df_3b10[df_3b10['Location'] == loc]
      
    # append the mean
    list_fmc_mean.append(np.mean(df_loc_f['FMC63_Epitope_Drop_sd_from_wt']))
    list_4g7_mean.append(np.mean(df_loc_4['4G7_Epitope_Drop_sd_from_wt']))
    list_3b10_mean.append(np.mean(df_loc_3['3B10_Epitope_Drop_sd_from_wt']))
    
    # get the number of winners
    win_f = np.count_nonzero(df_loc_f['FMC63_Epitope_Drop_sd_from_wt'] > 2)
    win_4 = np.count_nonzero(df_loc_4['4G7_Epitope_Drop_sd_from_wt'] > 2)
    win_3 = np.count_nonzero(df_loc_3['3B10_Epitope_Drop_sd_from_wt'] > 2)
    
    # get the number of obs
    total_f = df_loc_f['FMC63_Epitope_Drop_sd_from_wt'].shape[0]
    total_4 = df_loc_4['4G7_Epitope_Drop_sd_from_wt'].shape[0]
    total_3 = df_loc_3['3B10_Epitope_Drop_sd_from_wt'].shape[0]
    
    list_fmc_intol.append(win_f / total_f)
    list_4g7_intol.append(win_4 / total_4)
    list_3b10_intol.append(win_3 / total_3)


# from pymol on 6al5 (my numbering)
"""
select c. A within 7 of c. L

48, 82-86, 88, 204-206

select c. A within 7 of c. H

48, 75, 77, 83, 86, 88, 89, 99, 100,
 101, 105, 136, 137, 138, 140, 142, 143, 144
145, 146, 147, 149, 196-206

"""

list_7a = [48,75,77,82,83,84,85,86,88,89,99,100,
           101,105,136,137,138,140,142,143,144,145,
           146,147,149,196,197,198,199,200,201,202,
           203,204,205,206]


np_fmc_mean = np.array(list_fmc_mean)
np_4g7_mean = np.array(list_4g7_mean)
np_3b10_mean = np.array(list_3b10_mean)

np_fmc_intol = np.array(list_fmc_intol)
np_4g7_intol = np.array(list_4g7_intol)
np_3b10_intol = np.array(list_3b10_intol)