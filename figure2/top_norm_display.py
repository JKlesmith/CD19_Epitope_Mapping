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


# set the tolerance
tol = 2

# strip out WT and UNK mutations
strip1 = df_pact[df_pact['FMC63_Epitope_Display_classified'] != 'WT']
strip2 = strip1[~np.isnan(strip1['FMC63_Epitope_Display_sd_from_wt'])]    
df_search = strip2[~np.isnan(strip2['FMC63_Epitope_Top_sd_from_wt'])]    


# get fmc vs 4g7 of display
np_dis_d_n = df_search['FMC63_Epitope_Display_sd_from_wt'].values
np_dis_t_n = df_search['FMC63_Epitope_Top_sd_from_wt'].values



list_color = []

for i in range(np_dis_d_n.shape[0]):
    # y = x - tol

    if np_dis_t_n[i] < np_dis_d_n[i] - tol:
        list_color.append('red')
    else:
        list_color.append('black')
    

# plot
plt.figure(figsize=(9,9))
plt.scatter(np_dis_d_n, np_dis_t_n, s=0.5, c=list_color, alpha=0.5)
plt.show()








# get fmc vs 4g7 of display
np_dis_d_n = df_search['FMC63_Epitope_Display_sd_from_wt'].values
np_dis_t_n = df_search['FMC63_Epitope_Top_sd_from_wt'].values
np_top = np_dis_t_n - np_dis_d_n


list_color = []

for i in range(np_dis_d_n.shape[0]):
    # y = x - tol

    if np_top[i] < -1 * tol:
        list_color.append('red')
    else:
        list_color.append('black')
 

# plot
plt.figure(figsize=(9,9))
plt.scatter(np_dis_d_n, np_top, s=0.5, c=list_color, alpha=0.5)
plt.show()






# remove stops
mask_stops = ~(df_search['Mutation'] == '*')
np_top = np_top[mask_stops]



# get fmc vs 4g7 of display
np_drop = df_search['FMC63_Epitope_Drop_sd_from_wt'][mask_stops].values

# get the fraction burial
np_burial = df_search['frac_burial'][mask_stops].values

"""
# plot
plt.figure(figsize=(9,9))
plt.scatter(np_top, np_drop, s=0.5, c='black', alpha=0.5)
plt.show()
"""
mask = ~np.isnan(np_drop) & ~np.isnan(np_top)
"""
print(np.count_nonzero(
        (np_drop[mask] >= 2) &
        (np_top[mask] <= -2)
        ))

print(np.count_nonzero(
        (np_drop[mask] >= 2) &
        (np_top[mask] > -2)
        ))

print(np.count_nonzero(
        (np_drop[mask] < 2) &
        (np_top[mask] <= -2)
        ))
"""



# get a mask of high count display variants
mask_counts = (
        (df_search['FMC63_Epitope_Display_ref_counts'][mask_stops] >= 12) &
        (df_search['FMC63_Epitope_Display_sel_counts'][mask_stops] >= 12)
        )

np_drop = np_drop[mask_counts]
np_top = np_top[mask_counts]
np_burial = np_burial[mask_counts]

# plot
plt.figure(figsize=(9,9))
plt.scatter(np_top, np_drop, s=0.5, c='black', alpha=0.5)
plt.show()

mask = ~np.isnan(np_drop) & ~np.isnan(np_top)

# top left true hits
print(np.count_nonzero(
        (np_drop[mask] > 2) &
        (np_top[mask] < -2)
        ))

# top right, false hits drop
print(np.count_nonzero(
        (np_drop[mask] > 2) &
        (np_top[mask] >= -2)
        ))

# bottom left, false hits top
print(np.count_nonzero(
        (np_drop[mask] <= 2) &
        (np_top[mask] < -2)
        ))

# bottom right, true negative both
print(np.count_nonzero(
        (np_drop[mask] <= 2) &
        (np_top[mask] >= -2)
        ))

# 3797 mutations that pass the filters and had matching data


# report our mutations
out_top_core = np_top[(np_burial >= 0.9)]
out_drp_core = np_drop[(np_burial >= 0.9)]

out_top_surf = np_top[~(np_burial >= 0.9)]
out_drp_surf = np_drop[~(np_burial >= 0.9)]