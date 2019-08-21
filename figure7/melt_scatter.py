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

"""
Function: per-mutation scatter
"""
# set the tolerance
read_tol = 12
dset = 'FMC63'

# strip out WT, UNK, stop, 12 count mutations
mut_filter = (
        (df_pact[dset + '_Epitope_Drop_classified'] != 'WT') &
        (df_pact[dset + '_Epitope_Drop_classified'] != 'STOP') &
        (df_pact[dset + '_Epitope_Drop_classified'] != 'UNK') &
        (df_pact['fmc_melt_classified'] != 'UNK') &
        (df_pact['fmc_melt_ref_counts'] >= read_tol) &
        (df_pact[dset + '_Epitope_Drop_ref_counts'] >= read_tol) &
        (df_pact['Location'] < 259)
        )

# return a df that is filtered for high quality mutations
df_filtered = df_pact[mut_filter]

# return np_arrays
np_drp = df_filtered[dset + '_Epitope_Drop_sd_from_wt'].values
np_melt = df_filtered['fmc_melt_sd_from_wt'].values

# return the searched locations
np_loc = np.unique(df_filtered['Location'].astype(int))

# plot
plt.figure(figsize=(6,6))
plt.scatter(np_drp, np_melt, s=3, c='black', alpha=0.5)
plt.axhline(y=2, color='b', linestyle='--')
plt.axvline(x=2, color='b', linestyle='--')
plt.show()

"""
Function:
"""

list_pm_mean = []

# bin each site and calculate the mean z-score
for loc in np_loc:
        
    # mask the loc
    mask_loc = (df_filtered['Location'] == loc)
        
    # return a df at the site
    df_site = df_filtered[mask_loc]
    
    # get the site mean
    site_mean = np.mean(df_site[dset + '_Epitope_Drop_sd_from_wt'])
    #list_pm_mean.append([site_mean] * df_site.shape[0])
    
    # is there melt muts?
    count_wins = np.count_nonzero(df_site['fmc_melt_sd_from_wt'] > 2)
    
       
        
        
    # get the number of winners
    count_win = np.count_nonzero(df_site[dset + '_Epitope_Drop_sd_from_wt'] > 2)
    
    # get the number of obs
    count_obs = df_site[dset + '_Epitope_Drop_sd_from_wt'].shape[0]
    
    # calculate fraction
    if count_obs == 0:
        continue
    else:
        frac = count_win/count_obs
        list_pm_mean.append([frac] * count_obs)
        
        
    # output mutations
    if count_win > 0 and count_wins > 0:
        print(str(loc) + " " + str(round(frac,3)))
    
np_wins = np.hstack(list_pm_mean)
 
# plot
plt.figure(figsize=(6,6))
plt.scatter(np_wins, np_melt, s=3, c='black', alpha=0.5)
plt.axhline(y=2, color='b', linestyle='--')
plt.axvline(x=.25, color='b', linestyle='--')
plt.show()


