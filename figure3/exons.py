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


"""
Function: Pre-process dataset
"""
# strip out WT, UNK, STOPs, and low read counts
# RASA ≤ 0.15 ==> buried >= 0.85 frac_burial

# set dataset name
dset_name = 'FMC63'
#dset_name = '4G7'
#dset_name = '3B10'

show_surface = True


# get a list of frac burial
mask_wt = (df_pact[dset_name + '_Epitope_Drop_classified'] == 'WT')
df_wt = df_pact[mask_wt]

# preprocess the rest of the list
mask_drop = (
        (df_pact[dset_name + '_Epitope_Drop_classified'] != 'WT') &
        (df_pact[dset_name + '_Epitope_Drop_classified'] != 'UNK') &        
        (df_pact[dset_name + '_Epitope_Drop_classified'] != 'STOP') &        
        (df_pact[dset_name + '_Epitope_Drop_ref_counts'] >= 12)
        )

# get the cleaned up drop dataframe
df_drp = df_pact[mask_drop]


"""
Function: Find winners based on fraction of obs
"""
# setup lists of 0.XX fraction of observed
list_25 = []

# loop the locations
for loc in np_locations:
    
    # return a df for that location
    df_loc = df_drp[df_drp['Location'] == loc]
       
    # get the number of winners
    count_win = np.count_nonzero(df_loc[dset_name + '_Epitope_Drop_sd_from_wt'] > 2)
    
    # get the number of obs
    count_obs = df_loc[dset_name + '_Epitope_Drop_sd_from_wt'].shape[0]
    
    # calculate fraction
    if count_obs == 0:
        continue
    else:
        frac = count_win/count_obs
    
    # assign into buckets
    if frac >= 0.25:
        list_25.append(loc)

"""
Function: format lists to pymol selection api
""" 
def process_selection(list_resis):
    return '+'.join(map(str, set(list_resis)))

def pymol_resi_list_gen(start, end):
    return '+'.join(map(str, set([x for x in range(start, end + 1)])))

"""
Function: __main__
"""
#Open the pdb file
cmd.load(pact_path + '/' + '6al5_0001f_0001_0001.pdb')

# set a white background
cmd.bg_color(color="white")

# set a base color of gray
cmd.color('gray90', selection = 'all')

# set the colors to match the exons
cmd.color('orange', selection = 'resi ' + pymol_resi_list_gen(1,10))
cmd.color('green', selection = 'resi ' + pymol_resi_list_gen(11,99))
cmd.color('magenta', selection = 'resi ' + pymol_resi_list_gen(100,167))
cmd.color('gray90', selection = 'resi ' + pymol_resi_list_gen(168,259))
cmd.color('cyan', selection = 'resi ' + pymol_resi_list_gen(260,272))

# show as cartoon
cmd.show_as("cartoon")

# make sheets flat
cmd.set("cartoon_flat_sheets", value=1)

# color the winners
try:
    if len(list_25) > 0:
        cmd.show(representation='sticks', selection = 'resi ' + process_selection(list_25))

except pymol.parsing.QuietException:
    print("selection error, empty list?")
    quit()

# hide c-term
cmd.hide(selection = 'resi 259+260+261+262+263+264+265+266+267+268+269+270+271+272')

# set ray shadows to Zero
cmd.set("ray_shadows", value=0)

# hide hydrogens
cmd.remove('hydrogens')