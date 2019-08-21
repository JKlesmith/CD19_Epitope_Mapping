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

"""
Function: 
"""
# set the tolerance
z_tol = 2
read_tol = 12
dset = 'FMC63'

# strip out WT, UNK, stop, 12 count mutations
mut_filter = (
        (df_pact[dset + '_Epitope_Drop_classified'] != 'WT') &
        (df_pact[dset + '_Epitope_Drop_classified'] != 'STOP') &
        (df_pact['fmc_melt_classified'] != 'UNK') &
        (df_pact['fmc_melt_ref_counts'] >= read_tol)
        )

# return a df that is filtered for high quality mutations
df_filtered = df_pact[mut_filter]

# return np_arrays
np_drp = df_filtered[dset + '_20nm_Epitope_Drop_sd_from_wt'].values
np_melt = df_filtered['fmc_melt_sd_from_wt'].values

# return the searched locations
np_loc = np.unique(df_filtered['Location'].astype(int))

# plot
#x = np.linspace(-5, 5, 20)
plt.figure(figsize=(6,6))
plt.scatter(np_drp, np_melt, s=3, c='black', alpha=0.5)
plt.axhline(y=2, color='b', linestyle='--')
plt.axvline(x=2, color='b', linestyle='--')
#plt.plot(x, x+1, '-r')
#plt.plot(x, x, '-r')
#plt.plot(x, x-1, '-r')
plt.show()






list_any = []

# bin each site and calculate the mean z-score
for loc in np_loc:
        
    # mask the loc
    mask_loc = (df_filtered['Location'] == loc)
        
    # return a df at the site
    df_site = df_filtered[mask_loc]
    
    # find mutations
    melt = df_site['fmc_melt_sd_from_wt']

    # count muts
    if np.count_nonzero(melt > 2) > 0:
        list_any.append(loc)








"""
Function: format lists to pymol selection api
""" 
def process_selection(list_resis):
    return '+'.join(map(str, set(list_resis)))    

"""
Function: Pymol
"""


#Open the pdb file
cmd.load(pact_path + '/' + '6al5_0001f_0001_0001.pdb')

# set a white background
cmd.bg_color(color="white")

# set a base color of gray
cmd.color('gray90', selection = 'all')

# show as cartoon
cmd.show_as("cartoon")

# make sheets flat
cmd.set("cartoon_flat_sheets", value=1)

# color the winners
try:
   
    #if len(np_sites_top.tolist()) > 0:
        #cmd.color('red', selection = 'resi ' + process_selection(np_sites_top.tolist()))
        #cmd.show(representation='sticks', selection = 'resi ' + process_selection(np_sites_top.tolist()))
    """
    if len(list_100) > 0:
        cmd.color('red', selection = 'resi ' + process_selection(list_100))
        cmd.show(representation='sticks', selection = 'resi ' + process_selection(list_100))
        
    if len(list_75) > 0:
        cmd.color('magenta', selection = 'resi ' + process_selection(list_75))
        cmd.show(representation='sticks', selection = 'resi ' + process_selection(list_75))
        
    if len(list_50) > 0:
        cmd.color('cyan', selection = 'resi ' + process_selection(list_50))
        cmd.show(representation='sticks', selection = 'resi ' + process_selection(list_50))
        
    if len(list_25) > 0:
        cmd.color('blue', selection = 'resi ' + process_selection(list_25))
        cmd.show(representation='sticks', selection = 'resi ' + process_selection(list_25))
    """
    if len(list_any) > 0:
        cmd.color('red', selection = 'resi ' + process_selection(list_any))
        cmd.show(representation='sticks', selection = 'resi ' + process_selection(list_any))

        
except pymol.parsing.QuietException:
    print("selection error, empty list?")
    quit()

# hide c-term
cmd.hide(selection = 'resi 259+260+261+262+263+264+265+266+267+268+269+270+271+272')

# set ray shadows to Zero
cmd.set("ray_shadows", value=0)

# hide hydrogens
cmd.remove('hydrogens')
