# -*- coding: utf-8 -*-
"""
Created on Sat Jul 19 12:46:57 2025

@author: utente
"""

#%% Import the files
### - "Scenario_generation_function": 
###    - data processing and scenario generation by clustering in the training dataset (2005-2017)
###    - postprocessing
### - "Model_REC_function":
###    - clustering-based design-operation optimization in the training dataset (2005-2017)
###    - full timeseries operation-only optimization in the testing dataset (2018-2023)

import Scenario_generation_function as sg
import Model_REC_function as mr

### CHOICE
loc = "Padova"  # location
plot_ts = "no" # plot full timeseries
opt_K = 10      # "best" number of clusters and representative days

#%% ------------------- Scenario generation by clustering ---------------------

sg.scenario_generation(loc=loc,
                       year_start=2005,
                       year_end=2023,
                       seasons=0,
                       extreme_method='A', 
                       plot_ts=plot_ts) 

#%% ------------------- Post-processing clustering ----------------------------

index_new, index_representative_r_new, number_elements_each_cluster_new, \
TS, d_rts = sg.post_clustering(loc=loc,
                               seasons=0,
                               extreme_method='A', 
                               opt_K=opt_K)

#%% ------------------------- REC optimization --------------------------------
### Clustering-based design-operation optimization model (training dataset)

mr.clustering_design_operation(loc=loc,
                               seasons=0,
                               extreme_method='A') 

#%% ------------------------- REC optimization --------------------------------
### Full timeseries operation-only optimization model (testing dataset)

mr.fullts_operation(loc=loc,
                    seasons=0,
                    extreme_method='A',
                    opt_K=opt_K) 