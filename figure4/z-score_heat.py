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

# constant defines
mutation_types = '*FWYPMILVAGCSTNQDEHKR'

# find the unique locations
np_locations = np.sort(np.unique(df_pact['Location'].astype(int)))

# make lists to work into
str_top = ''
str_drp = ''
str_wt = ''

# loop the mutations in my types
for mut in mutation_types:
     
    # return a mask for that mutation, and sort on location
    df_mut = df_pact[(df_pact['Mutation'] == mut)].sort_values(by='Location')
       
    # check to see if all locations are represented
    if not np.array_equal(np_locations, df_mut['Location'].values):
        print("[Error] Expected locations mis-match.")
        print(np_locations)
        print(df_mut['Location'].values)
    
    # normalize the top
    np_top = (df_mut['3B10_Epitope_Top_sd_from_wt'].values - 
              df_mut['3B10_Epitope_Display_sd_from_wt'].values)
    
    # get the drop
    np_drp = df_mut['3B10_Epitope_Drop_sd_from_wt'].values
    
    # return a mask of read counts
    mask_counts = df_mut['3B10_Epitope_Drop_ref_counts'] < 12
    
    # set nan at sites that fail count filter
    np_top[mask_counts] = np.nan
    np_drp[mask_counts] = np.nan
    
    # append to our output string
    str_top += mut + ',' + ','.join(map(str, np_top.tolist())) + '\n'
    str_drp += mut + ',' + ','.join(map(str, np_drp.tolist())) + '\n'
    
# append to our w/t string (use the last mutation type df)
str_head_wt = ',' + ','.join(map(str, df_mut['wt_resi'].tolist())) + '\n'

# create the header location
str_head_loc = ',' + ','.join(map(str, np_locations.tolist())) + '\n'

# write our csv
with open(pact_path + '/z-score_top3B10.csv', 'w') as outfile:
    outfile.write(str_head_loc)
    outfile.write(str_head_wt)
    outfile.write(str_top)
    
# write our csv
with open(pact_path + '/z-score_drp3B10.csv', 'w') as outfile:
    outfile.write(str_head_loc)
    outfile.write(str_head_wt)
    outfile.write(str_drp)
