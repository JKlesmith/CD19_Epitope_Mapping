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

"""
FMC63
"""
# strip out WT, UNK, STOPs, and low read counts
mask_fmc = (
        (df_pact['FMC63_Epitope_Drop_classified'] != 'WT') &      
        (df_pact['FMC63_Epitope_Drop_classified'] != 'STOP') &  
        (df_pact['FMC63_Epitope_Drop_ref_counts'] >= 12) & 
        (df_pact['FMC63_Epitope_Drop_sd_from_wt'] > 2) &
        (df_pact['frac_burial'] < 0.85)
        )

core_fmc = np.count_nonzero(mask_fmc)
core_fmc_sites = np.unique(df_pact[mask_fmc]['Location'].astype(int))

list_fmc_frac = []

for site in core_fmc_sites:
    
    # get the number of BEN at this site
    count_ben = df_pact[mask_fmc][df_pact[mask_fmc]['Location'] == site].shape[0]
    
    # get the total at this site
    count_total = np.count_nonzero(
            (df_pact['Location'] == site) &
            (df_pact['FMC63_Epitope_Drop_classified'] != 'WT') &      
            (df_pact['FMC63_Epitope_Drop_classified'] != 'STOP') &  
            (df_pact['FMC63_Epitope_Drop_ref_counts'] >= 12)
            )
    
    list_fmc_frac.append(count_ben/count_total)

np_fmc_frac = np.array(list_fmc_frac)
print(np.count_nonzero(np_fmc_frac >= 0.25))
 


   
"""
4G7
"""
mask_4g7 = (
        (df_pact['4G7_Epitope_Drop_classified'] != 'WT') &      
        (df_pact['4G7_Epitope_Drop_classified'] != 'STOP') &  
        (df_pact['4G7_Epitope_Drop_ref_counts'] >= 12) & 
        (df_pact['4G7_Epitope_Drop_sd_from_wt'] > 2) &
        (df_pact['frac_burial'] < 0.85)
        )

core_4g7 = np.count_nonzero(mask_4g7)
core_4g7_sites = np.unique(df_pact[mask_4g7]['Location'].astype(int))

list_4g7_frac = []

for site in core_4g7_sites:
    
    # get the number of BEN at this site
    count_ben = df_pact[mask_4g7][df_pact[mask_4g7]['Location'] == site].shape[0]
    
    # get the total at this site
    count_total = np.count_nonzero(
            (df_pact['Location'] == site) &
            (df_pact['4G7_Epitope_Drop_classified'] != 'WT') &      
            (df_pact['4G7_Epitope_Drop_classified'] != 'STOP') &  
            (df_pact['4G7_Epitope_Drop_ref_counts'] >= 12)
            )
    
    list_4g7_frac.append(count_ben/count_total)

np_4g7_frac = np.array(list_4g7_frac)
print(np.count_nonzero(np_4g7_frac >= 0.25))

"""
3B10
"""
mask_3b10 = (
        (df_pact['3B10_Epitope_Drop_classified'] != 'WT') &      
        (df_pact['3B10_Epitope_Drop_classified'] != 'STOP') &  
        (df_pact['3B10_Epitope_Drop_ref_counts'] >= 12) & 
        (df_pact['3B10_Epitope_Drop_sd_from_wt'] > 2) &
        (df_pact['frac_burial'] < 0.85)
        )

core_3b10 = np.count_nonzero(mask_3b10)
core_3b10_sites = np.unique(df_pact[mask_3b10]['Location'].astype(int))

list_3b10_frac = []

for site in core_3b10_sites:
    
    # get the number of BEN at this site
    count_ben = df_pact[mask_3b10][df_pact[mask_3b10]['Location'] == site].shape[0]
    
    # get the total at this site
    count_total = np.count_nonzero(
            (df_pact['Location'] == site) &
            (df_pact['3B10_Epitope_Drop_classified'] != 'WT') &      
            (df_pact['3B10_Epitope_Drop_classified'] != 'STOP') &  
            (df_pact['3B10_Epitope_Drop_ref_counts'] >= 12)
            )
    
    list_3b10_frac.append(count_ben/count_total)

np_3b10_frac = np.array(list_3b10_frac)
print(np.count_nonzero(np_3b10_frac >= 0.25))





"""
FMC vs 4G7
"""
mask_f_site = (np_fmc_frac >= 0.25)
mask_4_site = (np_4g7_frac >= 0.25)
mask_3_site = (np_3b10_frac >= 0.25)

shared_f = core_fmc_sites[mask_f_site]
shared_4 = core_4g7_sites[mask_4_site]
shared_3 = core_3b10_sites[mask_3_site]

f_to_4 = np.count_nonzero(np.isin(shared_f, shared_4))
f_to_3 = np.count_nonzero(np.isin(shared_f, shared_3))