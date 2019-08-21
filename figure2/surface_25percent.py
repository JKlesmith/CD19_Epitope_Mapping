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
from itertools import product 
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

# constant defines
mutation_types = '*FWYPMILVAGCSTNQDEHKR'

# find the unique locations
np_locations = np.sort(np.unique(df_pact['Location'].astype(int)))


"""
Function: Pre-process dataset
"""
# strip out WT, UNK, STOPs, and low read counts
# RASA ≤ 0.15 ==> buried >= 0.85 frac_burial

# set dataset name
dset_name = 'FMC63'

# process the dataset
mask_drop = (
        (df_pact[dset_name + '_Epitope_Drop_classified'] != 'WT') &
        (df_pact[dset_name + '_Epitope_Drop_classified'] != 'UNK') &        
        (df_pact[dset_name + '_Epitope_Drop_classified'] != 'STOP') &        
        (df_pact[dset_name + '_Epitope_Drop_ref_counts'] >= 12) & 
        (df_pact['frac_burial'] < 0.85)
        )

df_drop = df_pact[mask_drop]

# find the unique locations
np_locations = np.unique(df_drop['Location'].astype(int))

# define three lists to work into
list_frac = []


# loop each location and calculate the frac of BEN at site vs max freq
for loc in np_locations:
    
    # return a df of the location
    df_loc = df_drop[df_drop['Location'] == loc]
    
    # return np arrs of BEN, non WT,UNK muts, and max site
    count_ben = np.count_nonzero(df_loc[dset_name + '_Epitope_Drop_classified'] == 'BEN')
    
    # total at site
    count_all = df_loc.values.shape[0]
    
    # add to list
    list_frac.append(count_ben/count_all)        

np_frac = np.array(list_frac)