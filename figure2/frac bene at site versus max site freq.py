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

# define three lists to work into
list_50 = []
list_51t74 = []
list_75t100 = []

# loop each location and calculate the frac of BEN at site vs max freq
for loc in np_locations:
    
    # return a df of the location
    df_loc = df_pact[df_pact['Location'] == loc]
    
    # return np arrs of BEN, non WT,UNK muts, and max site
    count_ben = np.count_nonzero(df_loc['fmc_clin_classified'] == 'BEN')
    
    count_all = np.count_nonzero(
            (df_loc['fmc_clin_classified'] != 'WT') &
            (df_loc['fmc_clin_classified'] != 'STOP') &
            (df_loc['fmc_clin_classified'] != 'UNK'))

    max_per = df_loc['max_percent'].values[0]
    
    # add to list
    if max_per < 50:
        list_50.append(count_ben/count_all)
    elif (max_per >= 50) and (max_per < 75):
        list_51t74.append(count_ben/count_all)        
    else: 
        list_75t100.append(count_ben/count_all)    

# make a boxplot of enriched, not-enriched, and basal burials (remove WT)
np_50 = np.array(list_50)
np_51t74 = np.array(list_51t74)
np_75t100 = np.array(list_75t100)

# make a dataframe that assigns groups
a = pd.DataFrame({'group':np.repeat('A',np_50.shape[0]), 'value':np_50})
b = pd.DataFrame({'group':np.repeat('B',np_51t74.shape[0]), 'value':np_51t74})
c = pd.DataFrame({'group':np.repeat('C',np_75t100.shape[0]), 'value':np_75t100})
df=a.append(b).append(c)
 
# Usual boxplot
sns.boxplot(x='group', y='value', data=df)
ax = sns.boxplot(x='group', y='value', data=df, color='white')
ax = sns.stripplot(x='group', y='value', data=df,
                   color="red", jitter=0.2, size=2.5)

# iterate over boxes
for i,box in enumerate(ax.artists):
    box.set_edgecolor('black')
    box.set_facecolor('white')

    # iterate over whiskers and median lines
    for j in range(6*i,6*(i+1)):
         ax.lines[j].set_color('black')