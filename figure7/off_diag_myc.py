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
dset = '4G7'

# strip out WT, UNK, stop, 12 count mutations
mut_filter = (
        (df_pact[dset + '_Epitope_Drop_classified'] != 'WT') &
        (df_pact[dset + '_Epitope_Drop_classified'] != 'STOP') &
        
        (df_pact[dset + '_Epitope_Drop_classified'] != 'UNK') &
        (df_pact[dset + '_20nm_Epitope_Myc_classified'] != 'UNK') &

        (df_pact[dset + '_Epitope_Drop_ref_counts'] >= read_tol) &
        (df_pact[dset + '_20nm_Epitope_Myc_ref_counts'] >= read_tol)
        )

# return a df that is filtered for high quality mutations
df_filtered = df_pact[mut_filter]

# return np_arrays
np_drp = df_filtered[dset + '_Epitope_Drop_sd_from_wt'].values
np_top = df_filtered[dset + '_20nm_Epitope_Myc_sd_from_wt'].values

# return the searched locations
np_loc = np.unique(df_filtered['Location'].astype(int))

# plot
x = np.linspace(-5, 5, 20)
plt.figure(figsize=(6,6))
plt.scatter(np_top, np_drp, s=3, c='black', alpha=0.5)
plt.axhline(y=2, color='b', linestyle='--')
plt.axvline(x=2, color='b', linestyle='--')
plt.plot(x, x+1, '-r')
plt.plot(x, x, '-r')
plt.plot(x, x-1, '-r')
plt.show()


list_any = []

# bin each site and calculate the mean z-score
for loc in np_loc:
        
    # mask the loc
    mask_loc = (df_filtered['Location'] == loc)
    
    # are we at the surface?
    #rasa = 1 - df_filtered[mask_loc]['frac_burial'].values[0]
    
    # return a df at the site
    df_site = df_filtered[mask_loc]
    
    # find mutations above 2z and off diag
    pm = df_site[dset + '_Epitope_Drop_sd_from_wt']
    nm = df_site[dset + '_20nm_Epitope_Myc_sd_from_wt']

    # find the abs change of mutations
    dset_abs = np.subtract(pm, nm)

    # find a mask for off diag
    off_diag_mask = dset_abs <= 1
    
    # for mutations off diag are any above 2z
    pm_od = pm[off_diag_mask]
    nm_od = nm[off_diag_mask]
    
    pm_od_count = np.count_nonzero(pm_od[pm_od > 2])
    nm_od_count = np.count_nonzero(nm_od[nm_od > 2])


    if pm_od_count > 0 and nm_od_count > 0:
        list_any.append(loc)


"""
Function: PDB Import
"""
# PDB
import_pdb = '6al5_0001f_0001_0001.pdb'

# read in our pact file    
if isfile(pact_path + '/' + import_pdb):
    
    # open the file
    with open(pact_path + '/' + import_pdb, 'r') as infile:
        list_lines = infile.read().splitlines()

#Create the dicts for the atom and hetatms
list_chain = []
list_loc = []
list_atom = []
list_resi = []
list_xcor = []
list_ycor = []
list_zcor = []

# iterate through each line
for line in list_lines:
    
    # break at TER
    if line.startswith("TER"):
        break

    # get the chain
    list_chain.append(line[21:22].strip())

    # get the residue number
    list_loc.append(int(line[22:26].strip()))

    # get the atom name
    list_atom.append(line[12:16].strip())
    
    # get the residue name
    list_resi.append(line[17:20].strip())

    # get the xyz coords
    list_xcor.append(float(line[30:38].strip()))
    list_ycor.append(float(line[38:46].strip()))
    list_zcor.append(float(line[46:54].strip()))

# create a dataframe in one shot
df_pdb = pd.DataFrame(
    {'chain': list_chain,
     'loc': list_loc,
     'resi': list_resi,
     'atom': list_atom,
     'xcoor': list_xcor,
     'ycoor': list_ycor,
     'zcoor': list_zcor
    })
    
# defines
chain = 'A'
epi_resi = [144, 140, 146, 202]

# use CB or CA in case of GLY
search_resi = ((df_pdb['atom'] == 'CB') |
              ((df_pdb['atom'] == 'CA') &
              (df_pdb['resi'] == 'GLY')))

# return the xyz coors of the pdb
np_xzy = df_pdb.loc[search_resi, ['xcoor', 'ycoor', 'zcoor']].values

# list of np arrays with dists
list_dists = []

# return the xyz coor of the epitope
for loc in epi_resi:
    search_epi = ((df_pdb['loc'] == loc) & 
                  (df_pdb['atom'] == 'CB'))
    
    # check for gly
    if np.count_nonzero(search_epi) == 0:
        search_epi = ((df_pdb['loc'] == loc) & 
                      (df_pdb['atom'] == 'CA'))        
    
    np_epi = df_pdb.loc[search_epi, ['xcoor', 'ycoor', 'zcoor']].values
    
    
    # calculate the distance, append up to resi 258
    list_dists.append(sp.spatial.distance.cdist(np_xzy, 
                                                np_epi, metric='euclidean')[:258])

# stack the lists and get the min
np_dist_min = np.min(np.hstack(list_dists), axis=1)


"""
Look at myc winners
"""
df_myc_only = df_filtered[df_filtered[dset + '_20nm_Epitope_Myc_sd_from_wt'] > 2]


list_rasa = 1 - df_myc_only['frac_burial'].values
list_hi = df_myc_only['hydropathy'].values
list_dist = [np_dist_min[x - 1] for x in df_myc_only['Location'].values]



