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
Function: pre-process the dataframes
"""
# mask out WT, UNK, STOPs, and enforce a count threshold in display
mask_top = (
        (df_pact['FMC63_Epitope_Display_classified'] != 'WT') &
        (df_pact['FMC63_Epitope_Display_classified'] != 'STOP') &
        (df_pact['FMC63_Epitope_Display_classified'] != 'UNK') &
        (df_pact['FMC63_Epitope_Top_classified'] != 'UNK') &
        (df_pact['FMC63_Epitope_Display_ref_counts'] >= 12) &
        (df_pact['FMC63_Epitope_Display_sel_counts'] >= 12)
        )

# make a new df
df_top = df_pact[mask_top]

# mask out WT, UNK, STOPs
mask_drop = (
        (df_pact['FMC63_Epitope_Drop_classified'] != 'WT') &
        (df_pact['FMC63_Epitope_Drop_classified'] != 'STOP') &
        (df_pact['FMC63_Epitope_Drop_classified'] != 'UNK')
        )

# make a new df
df_drop = df_pact[mask_drop]

"""
Function: Let's find out the number of mutations at the site
Drop <= 2, Top < -2

We want a Venn Diagram
    Overlaps of total sites ID'ed by the Q3 mutations
    Sites ID'ed by all BEN drps
    Sites ID'ed by double drp/top
    Three circle Venn Diagram like the fusion paper
"""
# define lists to work into
list_numq3 = [] # number of q3 mutations (from top df, strict counts)
list_totdrp = [] # number of BEN drp at that site (seperate df, not as strict)
list_totdbl = [] # number of q1 true hits (from top df, strict counts)

list_tot_muts_top = []
list_tot_muts_drp = []

# loop each location
for loc in np_locations:
    
    # mask a df of the location
    df_top_loc = df_top[df_top['Location'] == loc]
    df_drp_loc = df_drop[df_drop['Location'] == loc]
    
    # get a list of z-scores for top (normalized) and drop from top df
    np_top = (df_top_loc['FMC63_Epitope_Top_sd_from_wt'].values - 
              df_top_loc['FMC63_Epitope_Display_sd_from_wt'].values)
    np_drp = df_top_loc['FMC63_Epitope_Drop_sd_from_wt'].values
      
    # return a count of the number of q3 mutations
    numq3 = np.count_nonzero((np_top < -2) & (np_drp <= 2))
    if numq3 > 1:
        list_numq3.append(numq3)
    else:
        list_numq3.append(0)
    
    # return the total number of BEN drps
    #np_drp_all = df_drp_loc['FMC63_Epitope_Drop_sd_from_wt'].values
    list_tot_muts_drp.append(np_drp.shape[0])
    list_totdrp.append(np.count_nonzero((np_drp > 2)))
    
    # return the total number of hits
    list_totdbl.append(np.count_nonzero((np_top < -2) & (np_drp > 2)))
    
    # total number of obs mutations at site
    list_tot_muts_top.append(np_top.shape[0])

# change lists to np arrays
np_a = np.array(list_numq3)
np_b = np.array(list_totdrp)
np_e = np.array(list_tot_muts_drp)

# calcuate our Venn diagram
mask_both = np.count_nonzero((np_a > 0) & (np_b > 0))
mask_q3only = np.count_nonzero((np_a > 0) & (np_b == 0))
mask_benonly = np.count_nonzero((np_a == 0) & (np_b > 0))
mask_none = np.count_nonzero((np_a == 0) & (np_b == 0))


print(str(mask_both) + "," + str(mask_q3only) + "," +
      str(mask_benonly) + "," +str(mask_none))






"""
# make a dataframe that assigns groups
a = pd.DataFrame({'group':np.repeat('A',np_a.shape[0]), 'value':np_a})
b = pd.DataFrame({'group':np.repeat('B',np_b.shape[0]), 'value':np_b})
c = pd.DataFrame({'group':np.repeat('C',np_c.shape[0]), 'value':np_c})
d = pd.DataFrame({'group':np.repeat('D',np_d.shape[0]), 'value':np_d})
df=a.append(b).append(c).append(d)
 
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

plt.show()
"""