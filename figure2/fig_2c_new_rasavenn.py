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

# find the unique locations
np_locations = np.unique(df_pact['Location'].astype(int))

"""
Function: basal wt rasa
"""
# find wt
mask_wt = (df_pact['FMC63_Epitope_Drop_classified'] == 'WT')

# make a new df
df_wt = df_pact[mask_wt]

# return the basal dataset
np_basal = np.array(1 - df_wt['frac_burial'].values)


"""
Function: per-site rasa
"""
# set dataset name
dset_name = 'FMC63'

# mask out WT, UNK, STOPs
mask_drop = (
        (df_pact[dset_name + '_Epitope_Drop_classified'] != 'WT') &
        (df_pact[dset_name + '_Epitope_Drop_classified'] != 'STOP') &
        (df_pact[dset_name + '_Epitope_Drop_classified'] != 'UNK') &
        (df_pact[dset_name + '_Epitope_Drop_ref_counts'] >= 12)
        )

# make a new df
df_drop = df_pact[mask_drop]

# define lists to work into
list_24 = []
list_49 = []
list_74 = []
list_99 = []
list_100 = []

list_all_rasa = []
list_all_frac = []

# loop each location
for loc in np_locations:
    
    # mask a df of the location
    df_loc = df_drop[df_drop['Location'] == loc]

    # get a list of z-scores
    np_drp = df_loc[dset_name + '_Epitope_Drop_sd_from_wt'].values
      
    # winners
    numdrp = np.count_nonzero(np_drp > 2)
    numtot = np_drp.shape[0]
        
    # fraction won
    win_frac = numdrp / numtot
    
    # rasa
    rasa = 1 - df_loc['frac_burial'].values[0]

    # bin the results
    if win_frac > 0 and win_frac < 0.25:
        list_24.append(rasa)
    elif win_frac >= 0.25 and win_frac < 0.5:    
        list_49.append(rasa)
    elif win_frac >= 0.5 and win_frac < 0.75:    
        list_74.append(rasa)
    elif win_frac >= 0.75 and win_frac < 1:    
        list_99.append(rasa)
    elif win_frac == 1:
        list_100.append(rasa)
        
    list_all_rasa.append(rasa)
    list_all_frac.append(win_frac)

    
    
np_24 = np.array(list_24)
np_49 = np.array(list_49)
np_74 = np.array(list_74)
np_99 = np.array(list_99)
np_100 = np.array(list_100)

np_all_rasa = np.array(list_all_rasa)
np_all_frac = np.array(list_all_frac)