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
    sns.set(rc={'figure.figsize':(9,9),
                'axes.facecolor':'white',
                'figure.facecolor':'white'
                })
except ImportError:
    exit("[Error] Seaborn not installed and is required.")
    
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

# constant defines
mutation_types = '*FWYPMILVAGCSTNQDEHKR'

# read in our pact file    
if isfile(pact_path + '/' + import_features):
    
    # load the csv of features into pandas
    df_pact = pd.read_csv(pact_path + '/' + import_features)

# get a mask to remove wild-type and stops
mask = ((df_pact['FMC63_20nm_Epitope_Drop_classified'] != "WT") &
        (df_pact['FMC63_20nm_Epitope_Drop_classified'] != "STOP"))

# define seperate np array
np_20drp_f = df_pact['FMC63_20nm_Epitope_Drop_sd_from_wt'][mask].values.T
np_20drp_4 = df_pact['4G7_20nm_Epitope_Drop_sd_from_wt'][mask].values.T

np_20myc_f = df_pact['FMC63_20nm_Epitope_Myc_sd_from_wt'][mask].values.T
np_20myc_4 = df_pact['4G7_20nm_Epitope_Myc_sd_from_wt'][mask].values.T

np_location = df_pact['Location'][mask].values.T




"""
nM Drop
"""
# create a mask where it is not nan
mask = ~np.isnan(np_20drp_f) & ~np.isnan(np_20drp_4)

# calculate the regression
_, _, r_val, _, _ = sp.stats.linregress(np_20drp_f[mask], np_20drp_4[mask])
print(round(r_val**2, 3))
print(np_20drp_f[mask].shape[0])

# plot
plt.figure(figsize=(6,6))
plt.scatter(np_20drp_f, np_20drp_4, s=0.5, c='black', alpha=0.7)
plt.show()

"""
pM Myc
"""
# create a mask where it is not nan
mask = ~np.isnan(np_20myc_f) & ~np.isnan(np_20myc_4)

# calculate the regression
_, _, r_val, _, _ = sp.stats.linregress(np_20myc_f[mask], np_20myc_4[mask])
print(round(r_val**2, 3))
print(np_20myc_f[mask].shape[0])

# plot
plt.figure(figsize=(6,6))
plt.scatter(np_20myc_f, np_20myc_4, s=0.5, c='black', alpha=0.7)
plt.show()



"""
remove tile 3 181 to 272
"""
# create a mask where it is not nan
mask = (~np.isnan(np_20myc_f) & ~np.isnan(np_20myc_4) & (np_location < 181))

# calculate the regression
_, _, r_val, _, _ = sp.stats.linregress(np_20myc_f[mask], np_20myc_4[mask])
print(round(r_val**2, 3))
print(np_20myc_f[mask].shape[0])

# plot
plt.figure(figsize=(6,6))
plt.scatter(np_20myc_f[mask], np_20myc_4[mask], s=0.5, c='black', alpha=0.7)
plt.show()