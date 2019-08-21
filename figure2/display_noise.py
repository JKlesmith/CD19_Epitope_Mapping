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
Fig SI X
A) FMC vs 4G7 display z-score (ID red scatter)
B) FMC sel vs ref counts (same muts from A in red == low read counts)
C) FMC vs 4G7 top z-score (red mutns removed, 
     mutations that are expected to be in epitope in purple)
D) FMC Top vs Display (red muts removed, 
     muts that are expected to be epitope in bright blue)
E) Table cross comparing the mutations in (C) and (D)

    
    
Do D next as I want to see how filtering display changes that coorlelation    
    
"""



# set the tolerance
tol = 2.5

# strip out WT and UNK mutations
strip1 = df_pact[df_pact['FMC63_Epitope_Display_classified'] != 'WT']
strip2 = strip1[~np.isnan(strip1['FMC63_Epitope_Display_sd_from_wt'])]    
df_search = strip2[~np.isnan(strip2['4G7_Epitope_Display_sd_from_wt'])]    

# isolate non-noise
nonnoise = ((df_search['FMC63_Epitope_Display_sd_from_wt'] < tol) &
                  (df_search['FMC63_Epitope_Display_sd_from_wt'] > -1 * tol) &
                  (df_search['4G7_Epitope_Display_sd_from_wt'] < tol) &
                  (df_search['4G7_Epitope_Display_sd_from_wt'] > -1 * tol))

# isolate noise
noise = ~((df_search['FMC63_Epitope_Display_sd_from_wt'] < tol) &
              (df_search['FMC63_Epitope_Display_sd_from_wt'] > -1 * tol) &
              (df_search['4G7_Epitope_Display_sd_from_wt'] < tol) &
              (df_search['4G7_Epitope_Display_sd_from_wt'] > -1 * tol))

"""
A) FMC vs 4G7 display z-score (ID red scatter)
"""
# get fmc vs 4g7 of display
np_dis_f_n = df_search['FMC63_Epitope_Display_sd_from_wt'][noise].values
np_dis_4_n = df_search['4G7_Epitope_Display_sd_from_wt'][noise].values
np_dis_f_nn = df_search['FMC63_Epitope_Display_sd_from_wt'][nonnoise].values
np_dis_4_nn = df_search['4G7_Epitope_Display_sd_from_wt'][nonnoise].values


dis_f_combined = np.concatenate((np_dis_f_nn, np_dis_f_n))
dis_4_combined = np.concatenate((np_dis_4_nn, np_dis_4_n))

list_color = np.concatenate((
                np.repeat('black', np_dis_f_nn.shape[0]),
                np.repeat('red', np_dis_f_n.shape[0])
                ))

# plot
plt.figure(figsize=(9,9))
plt.scatter(dis_f_combined, dis_4_combined, s=0.5, c=list_color, alpha=0.5)
plt.show()


"""
B) FMC sel vs ref counts (same muts from A in red == low read counts)
"""

# counts
np_fr_n = df_search['FMC63_Epitope_Display_ref_counts'][noise].values.astype(float)
np_fs_n = df_search['FMC63_Epitope_Display_sel_counts'][noise].values.astype(float)
#np_4r_n = df_search['4G7_Epitope_Display_ref_counts'][noise].values
#np_4s_n = df_search['4G7_Epitope_Display_sel_counts'][noise].values

# counts
np_fr_nn = df_search['FMC63_Epitope_Display_ref_counts'][nonnoise].values.astype(float)
np_fs_nn = df_search['FMC63_Epitope_Display_sel_counts'][nonnoise].values.astype(float)
#np_4r_nn = df_search['4G7_Epitope_Display_ref_counts'][nonnoise].values
#np_4s_nn = df_search['4G7_Epitope_Display_sel_counts'][nonnoise].values

# replace 0 values with 0.9 for log plotting
np_fr_n[np_fr_n == 0] = 0.9
np_fs_n[np_fs_n == 0] = 0.9

np_fr_nn[np_fr_nn == 0] = 0.9
np_fs_nn[np_fs_nn == 0] = 0.9

ref_combined = np.concatenate((np_fr_nn, np_fr_n))
sel_combined = np.concatenate((np_fs_nn, np_fs_n))

list_color = np.concatenate((
                np.repeat('black', np_fr_nn.shape[0]),
                np.repeat('red', np_fr_n.shape[0])
                ))

# plot
plt.figure(figsize=(9,9))
plt.scatter(ref_combined, sel_combined, s=0.5, c=list_color, alpha=0.5)
plt.yscale('log')
plt.xscale('log')
plt.show()



"""
D) FMC Top vs Display (red muts removed, 
     muts that are expected to be epitope in bright blue)
"""

# get fmc vs 4g7 of display
np_dis_f_n = df_search['FMC63_Epitope_Display_sd_from_wt'][noise].values
np_dis_4_n = df_search['FMC63_Epitope_Top_sd_from_wt'][noise].values
np_dis_f_nn = df_search['FMC63_Epitope_Display_sd_from_wt'][nonnoise].values
np_dis_4_nn = df_search['FMC63_Epitope_Top_sd_from_wt'][nonnoise].values


dis_f_combined = np.concatenate((np_dis_f_nn, np_dis_f_n))
dis_4_combined = np.concatenate((np_dis_4_nn, np_dis_4_n))

list_color = np.concatenate((
                np.repeat('black', np_dis_f_nn.shape[0]),
                np.repeat('red', np_dis_f_n.shape[0])
                ))

# plot
plt.figure(figsize=(9,9))
plt.scatter(dis_f_combined, dis_4_combined, s=0.5, c=list_color, alpha=0.5)
plt.show()