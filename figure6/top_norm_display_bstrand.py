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
    
"""
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
"""

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
        
        (df_pact[dset + '_Epitope_Drop_classified'] != 'UNK') &
        (df_pact[dset + '_Epitope_Top_classified'] != 'UNK') &
        (df_pact[dset + '_Epitope_Display_classified'] != 'UNK') &

        (df_pact[dset + '_Epitope_Display_ref_counts'] >= read_tol) &
        (df_pact[dset + '_Epitope_Display_sel_counts'] >= read_tol)
        )

# return a df that is filtered for high quality mutations
df_filtered = df_pact[mut_filter]

# return np_arrays
np_drp = df_filtered[dset + '_Epitope_Drop_sd_from_wt'].values
np_top = (df_filtered[dset + '_Epitope_Top_sd_from_wt'].values - 
          df_filtered[dset + '_Epitope_Display_sd_from_wt'].values)

# return the searched locations
np_loc = np.unique(df_filtered['Location'].astype(int))


list_color = []
list_drp_mean = []
list_top_mean = []

# bin each site and calculate the mean z-score
for loc in np_loc:
    
    # mask the loc
    mask_loc = (df_filtered['Location'] == loc)
      
    # append the mean
    #list_drp_mean.append(np.mean(np_drp[mask_loc]))
    #list_top_mean.append(np.mean(np_top[mask_loc]))
    
    drp_win = np.count_nonzero(np_drp[mask_loc] > 2)
    top_win = np.count_nonzero(np_top[mask_loc] < -2)
    
    if drp_win/np_drp[mask_loc].shape[0] >= 0.25:
        list_color.append(['red'] * np_drp[mask_loc].shape[0])
    else:
        list_color.append(['black'] * np_drp[mask_loc].shape[0])
    #top_win/np_top[mask_loc].shape[0]
    
    
#np_drp = np.array(list_drp_mean)
#np_top = np.array(list_top_mean)
    

list_c = np.hstack(list_color).tolist()

# plot
plt.figure(figsize=(6,6))
plt.scatter(np_top, np_drp, s=1, c=list_c, alpha=0.5)
plt.axhline(y=z_tol, color='b', linestyle='--')
plt.axvline(x=-z_tol, color='b', linestyle='--')
plt.show()


"""
What additional sites are implicated by the 51% of mutations in Q3?
Also, which sites were missed by the top gate (the balance of the 85% TPR)? 
I’d focus the analysis on these sites; 
i.e. would there have been a functionally different conclusion 
if we had used top instead of drop?


Definition: >= 25%
"""


site_drp = []
site_top = []
site_q3_win = []
site_q3_dblwin = []
site_q3_lose = []
site_both = []

# bin each site and calculate the mean z-score
for loc in np_loc:
    
    # mask the loc
    mask_loc = (df_filtered['Location'] == loc)
    
    # get the wins
    drp_win = np.count_nonzero(np_drp[mask_loc] > 2)
    top_win = np.count_nonzero(np_top[mask_loc] < -2)
    
    # what is the winning fraction
    frac_drp = drp_win/np_drp[mask_loc].shape[0]
    frac_top = top_win/np_top[mask_loc].shape[0]
    
    # do we want to force a min number of observed
    if np_drp[mask_loc].shape[0] < 10 or np_top[mask_loc].shape[0] < 10:
        continue
    
    # ID sites in Q3 not ID'ed in top left
    df_site = df_filtered[mask_loc]
    search = ((df_site[dset + '_Epitope_Drop_sd_from_wt'] > 2) & 
              (df_site[dset + '_Epitope_Top_sd_from_wt'] < -2))
    
    # do we have mutations in Q3?
    if np.count_nonzero(search) > 0:
        
        # are these mutations at sites not ID'ed in drp but in top?
        if frac_drp < 0.25 and frac_top >= 0.25:
            site_q3_win.append(loc)
            
        # are these mutations at sites not ID'ed in drp but in top?
        if frac_drp < 0.25 and frac_top < 0.25:
            site_q3_lose.append(loc)
            
        if frac_drp >= 0.25 and frac_top >= 0.25:
            site_q3_dblwin.append(loc)

    # find sites missed by top gate      
    if frac_drp >= 0.25 and frac_top < 0.25:
        print('fail drp, win top' + str(loc))
        
    # get a list of sites IDed as important
    if frac_drp >= 0.25:
        site_drp.append(loc)
        
    if frac_top >= 0.25:
        site_top.append(loc)
    
    if frac_drp >= 0.25 and frac_top >= 0.25:
        site_both.append(loc)
        
    #if frac_drp < 0.25 and frac_top >= 0.25:
    #    site_both.append(loc)


# ID sites in/not in either
np_comp = np.isin(site_top, site_drp)

# new sites
np_sites_top = np.array(site_top)[~np_comp]

"""
Import beta-strand information
"""

# set our filename
import_strand = 'res_bstrandcount.csv'

# read in our pact file    
if isfile(pact_path + '/' + import_strand):
        
    # load the csv of features into pandas
    df_strand = pd.read_csv(pact_path + '/' + import_strand)  




# match the new sites to the index in strand
dist_index_new = []

# interate each location
for loc in np_sites_top:
    
    # exclude locs above 258
    if loc > 258:
        continue
    
    # return a mask for the location in the df
    mask_loc = df_strand['Residue'] == loc

    # append the b-strand location at the new location
    dist_index_new.append(df_strand[mask_loc]['# from end'].values[0])
           
                          
                          
# match the new sites to the index in strand
dist_index_drp = []

# interate each location
for loc in site_drp:
    
    # exclude locs above 258
    if loc > 258:
        continue
    
    # return a mask for the location in the df
    mask_loc = df_strand['Residue'] == loc

    # append the b-strand location at the new location
    dist_index_drp.append(df_strand[mask_loc]['# from end'].values[0])



# match the new sites to the index in strand
dist_index_top = []

# interate each location
for loc in site_top:
    
    # exclude locs above 258
    if loc > 258:
        continue
    
    # return a mask for the location in the df
    mask_loc = df_strand['Residue'] == loc

    # append the b-strand location at the new location
    dist_index_top.append(df_strand[mask_loc]['# from end'].values[0])





