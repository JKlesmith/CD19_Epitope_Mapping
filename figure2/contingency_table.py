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

# constant defines
mutation_types = '*FWYPMILVAGCSTNQDEHKR'

# read in our pact file    
if isfile(pact_path + '/' + import_features):
        
    # load the csv of features into pandas
    df_pact = pd.read_csv(pact_path + '/' + import_features)     


"""
pM Top vs pM drop
"""
# get a mask to remove wild-type
mask_wt = (df_pact['FMC63_Epitope_Drop_classified'] != "WT")

# define seperate np arrays
np_drp_f = df_pact['FMC63_Epitope_Drop_sd_from_wt'][mask_wt].values.T
np_top_f = df_pact['FMC63_Epitope_Top_sd_from_wt'][mask_wt].values.T

# create a mask where it is not nan
mask = ~np.isnan(np_top_f) & ~np.isnan(np_drp_f)

# calculate the regression
_, _, r_val, _, _ = sp.stats.linregress(np_top_f[mask], np_drp_f[mask])
print(round(r_val**2, 3))
print(np_top_f[mask].shape[0])

# plot
plt.figure(figsize=(6,6))
plt.scatter(np_top_f, np_drp_f, s=0.5, c='black', alpha=0.7)
plt.show()


"""
Make a contingency table

y-axis: drop
x-axis: top

mask_xy (where xy is dnb for del, neu, ben)

"""    
# what are our possible states
list_states = ['DEL', 'NEU', 'BEN']

# get all of the combinations
list_prod = list(product(list_states, list_states))

# loop the combos
for combo in list_prod:
    
    # create a mask
    mask = ((df_pact['FMC63_Epitope_Drop_classified'] == combo[0]) & 
            (df_pact['FMC63_Epitope_Top_classified'] == combo[1]))
    
    # count trues
    print(' '.join(['Drop', combo[0], '-',
                    'Top', combo[1], '-',            
                    str(np.count_nonzero(mask))
                    ]))

    
"""
Look at the Drop:DEL and Top:DEL mutations
258 of them exist in this cross-classification
"""
# create a mask
mask = ((df_pact['FMC63_Epitope_Drop_classified'] == 'DEL') & 
        (df_pact['FMC63_Epitope_Top_classified'] == 'DEL'))

# return a table with the search
df_search = df_pact[mask]

# of this subset, how many have DEL display
# for FMC of 258, 50 have bad display (~19%)
print(np.count_nonzero(df_search['FMC63_Epitope_Display_classified'] == 'DEL'))


# exclude poor displayers
# create a mask
mask = ((df_pact['FMC63_Epitope_Drop_classified'] == 'DEL') & 
        (df_pact['FMC63_Epitope_Top_classified'] == 'DEL') & 
        (df_pact['FMC63_Epitope_Display_classified'] != 'DEL'))

# return a table with the search
df_search = df_pact[mask]



"""
The hypothesis is that the remaining muts are a sensitivity issue

questions:
    what is the dist. of mutations at the sites
    i.e.
    num of muts classified as D in top and drp
    what is the per-site 'missing' top mut count
    
"""

# find the unique locations
np_locations = np.unique(df_search['Location'].astype(int))

list_count = []
list_fract = []

list_count_drp = []
list_fract_drp = []

# loop the search locations
for l in np_locations:
    
    # get the total number of muts at that site
    msk_tot = ((df_pact['FMC63_Epitope_Top_classified'] != 'STOP') & 
            (df_pact['FMC63_Epitope_Top_classified'] != 'WT') & 
            (df_pact['FMC63_Epitope_Top_classified'] != 'UNK') &
            (df_pact['Location'] == l)
            )
    count_tot = np.count_nonzero(msk_tot)
    
    # get the total number of epitopes at that site
    msk_tot_top = ((df_pact['FMC63_Epitope_Top_classified'] == 'DEL') &
            (df_pact['Location'] == l)
            )
    count_tot_top = np.count_nonzero(msk_tot_top)
    
    # append to our list
    list_count.append(count_tot_top)
    list_fract.append(count_tot_top/count_tot)
    
    
    
    # get the total number of muts at that site
    msk_tot = ((df_pact['FMC63_Epitope_Drop_classified'] != 'STOP') & 
            (df_pact['FMC63_Epitope_Drop_classified'] != 'WT') & 
            (df_pact['FMC63_Epitope_Drop_classified'] != 'UNK') &
            (df_pact['Location'] == l)
            )
    count_tot = np.count_nonzero(msk_tot)
    
    # get the total number of epitopes at that site
    msk_tot_drp = ((df_pact['FMC63_Epitope_Drop_classified'] == 'BEN') &
            (df_pact['Location'] == l)
            )
    count_tot_drp = np.count_nonzero(msk_tot_drp)
    
    # append to our list
    list_count_drp.append(count_tot_drp)
    list_fract_drp.append(count_tot_drp/count_tot)
    
    






# make a dataframe that assigns groups
a = pd.DataFrame({'group':np.repeat('A',len(list_count)), 'value':list_count})
b = pd.DataFrame({'group':np.repeat('B',len(list_count_drp)), 'value':list_count_drp})
df=a.append(b)
 
# Usual boxplot
sns.boxplot(x='group', y='value', data=df)
ax = sns.boxplot(x='group', y='value', data=df, color='white')
ax = sns.stripplot(x='group', y='value', data=df, color="red", jitter=0.2, size=2.5)

# iterate over boxes
for i,box in enumerate(ax.artists):
    box.set_edgecolor('black')
    box.set_facecolor('white')

    # iterate over whiskers and median lines
    for j in range(6*i,6*(i+1)):
         ax.lines[j].set_color('black')







