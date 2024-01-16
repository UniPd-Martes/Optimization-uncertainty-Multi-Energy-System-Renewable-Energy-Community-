# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 09:37:45 2024

@author: utente
"""

#%% MODEL

# Mono-objective optimization of a local Multi Energy System (MES), 
# representing an Energy Community (EC) for the "building" sector
# Members with electricity and heating demands:
# - residential, house (res)
# - tertiary, small office (ter)
# - commercial, small restaurant (com)
# - public, primary school (pub)

### Renewable Energy Community (REC) model
# Individual self-consumption for each member
# Virtual energy sharing (net electricity) among members
# Heating demand satisfied by each member independently

### Energy conversion and storage units
# Energy conversion systems:
    # - Photovoltaic (PV)
    # - Boiler (BOIL), fuelled by natural gas
    # - Air-water heat pump (HP)
# Storage systems:
    # - Thermal energy storage (TES): hot water tank
    # - Electrical energy storage (EES): lithium battery
### How are these units distributed among the members?
# - res (prosumer): PV, EES, HP, TES
# - ter (consumer): boiler
# - com (consumer): boiler
# - pub (prosumer): PV, EES, HP, TES

### The model is solved in each year of the training dataset (2005-2014)
### - Without uncertainty: 4 average seasonal days in each year
### - With uncertainty of solar irradiance and ambient temperature: from 2 to 7 stochastic scenarios for each season in each year
#### Example with uncertainty:
#### First case: 2 scenarios for each season --> 8 scenarios in one year
#### Last case: 7 scenarios for each season --> 28 scenarios in one year 

#%% PYTHON PACKAGES

import gurobipy as grb            # Gurobi, for optimization
from gurobipy import *            # Gurobi, for optimization
import numpy as np                # Numpy, for mathematics
import pandas as pd               # Pandas, for data analysis
import matplotlib.pyplot as plt   # Pyplot, for plots, graphs, diagrams
import time                       # time required to solve the optimization
import os                         # to select the folder for saving

#%% TECHNO-ECONOMIC PARAMETERS 

PATH_data = "C:/Users/utente/Desktop/Università Gabriele/Assegno_biennale/Models/Model_ML2/Input_data/TechnoEconomicData.xlsx"

# TECHNO-ECONOMIC PARAMETERS OF EACH TECHNOLOGY
TEparam = pd.read_excel(PATH_data,
                        sheet_name='Conversion&StorageSystems',
                        index_col=[0,1])

# PV

# Investment cost in [€/kWp]
invPV = TEparam.loc["PV","InvestmentCost"]["Value"]
# Fixed O&M cost in [%investment/year]
OnM_fix_PV = TEparam.loc["PV","O&Mcost,fix"]["Value"]
OnM_fix_PV = OnM_fix_PV/100 # form [%] to []
# Lifetime in [years]
nPV = TEparam.loc["PV","Lifetime"]["Value"]

# BOIL

# Investment cost in [€/kWth]
invBOIL = TEparam.loc["BOIL","InvestmentCost"]["Value"]
# Fixed O&M cost in [%investment/year]
OnM_fix_BOIL = TEparam.loc["BOIL","O&Mcost,fix"]["Value"]
OnM_fix_BOIL = OnM_fix_BOIL/100 # form [%] to []
# Lifetime in [years]
nBOIL = TEparam.loc["BOIL","Lifetime"]["Value"]
# Thermal efficiency in []
eff_th_BOIL = TEparam.loc["BOIL","EfficiencyTh"]["Value"]

# HP

# Investment cost in [€/kWth]
invHP = TEparam.loc["HP","InvestmentCost"]["Value"]
# Fixed O&M cost in [%investment/year]
OnM_fix_HP = TEparam.loc["HP","O&Mcost,fix"]["Value"]
OnM_fix_HP = OnM_fix_HP/100 # form [%] to []
# Lifetime in [years]
nHP = TEparam.loc["HP","Lifetime"]["Value"]
# Supply temperature in [°C]
Tsupply = 70 
# Power consumption linearization: P(Q) = (e*Q + f)/COPideal [kW]
ee = TEparam.loc["HP","COPid*P(Q)_var"]["Value"]
ff = TEparam.loc["HP","COPid*P(Q)_fix"]["Value"]
# Minimum load in part load operation in [%] of the nominal capacity
min_load_HP = TEparam.loc["HP","MinLoad"]["Value"]
min_load_HP = min_load_HP/100 # from [%] to []

# TES

# Investment cost in [€/kWh]
invTES = TEparam.loc["TES","InvestmentCost"]["Value"]
# Fixed O&M cost in [%investment/year]
OnM_fix_TES = TEparam.loc["TES","O&Mcost,fix"]["Value"]
OnM_fix_TES = OnM_fix_TES/100 # form [%] to []
# Lifetime in [years]
nTES = TEparam.loc["TES","Lifetime"]["Value"]
# Round trip efficiency in []
rteTES = TEparam.loc["TES","RoundTripEfficiency"]["Value"]
# Self-discharge rate in [%SOC/hour], with SOC = state of charge
sdTES = TEparam.loc["TES","SelfDischarge"]["Value"]
sdTES = sdTES/100 # from [%] to []
# Output capacity in [kW/kWh]
coutTES = TEparam.loc["TES","OutputCapacity"]["Value"]
# Input capacity in [kW/kWh]
cinpTES = TEparam.loc["TES","InputCapacity"]["Value"]

# EES

# Investment cost in [€/kWh]
invEES = TEparam.loc["EES","InvestmentCost"]["Value"]
# Fixed O&M cost in [%investment/year]
OnM_fix_EES = TEparam.loc["EES","O&Mcost,fix"]["Value"]
OnM_fix_EES = OnM_fix_EES/100 # form [%] to []
# Lifetime in [years]
nEES = TEparam.loc["EES","Lifetime"]["Value"]
# Round trip efficiency in []
rteEES = TEparam.loc["EES","RoundTripEfficiency"]["Value"]
# Self-discharge rate in [%SOC/hour], with SOC = state of charge
sdEES = TEparam.loc["EES","SelfDischarge"]["Value"]
sdEES = sdEES/100 # from [%] to []
# Output capacity in [kW/kWh]
coutEES = TEparam.loc["EES","OutputCapacity"]["Value"]
# Input capacity in [kW/kWh]
cinpEES = TEparam.loc["EES","InputCapacity"]["Value"]

# PRICES OF ENERGY CARRIERS
ECparam = pd.read_excel(PATH_data,
                        sheet_name='Carriers',
                        index_col=[0,1])
# Natural gas
# Purchasing price of natural gas in [€/Mwh]
cbuy_ng = ECparam.loc["NatGas","CostHouseholds"]["Value"]
cbuy_ng = cbuy_ng/1000 # from [€/MWh] to [€/kWh]

# FINANCIAL PARAMETERS (INTEREST RATE, INCENTIVE FOR THE SHARED ENERGY IN THE REC)
Fparam = pd.read_excel(PATH_data,
                       sheet_name='Parameters',
                       index_col=[0])
r = Fparam.loc["InterestRate"]["Value"]
inc = Fparam.loc["REC_incentive"]["Value"]

# MAXIMUM HOURLY VARIATION OF FLEXIBLE DEMANDS 
D_var = 0.1

#%% DESIGN-OPERATION OPTIMIZATION MODEL

array=[] # to save the optimal values of the decision variables
    
def REC_MES_typical_day(K): # K is the number of days considered in one year

    # Initialize the Gurobi model
    m = grb.Model()
    
    # Set extra termination criteria
    m.Params.MIPGap = 2e-2        # tolerance of the optimal solution in []
    
    ### DECISION VARIABLES for the community
    # Total net electrical power exchanged with the grid
    # Imported power
    P_imp = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="P_imp")
    # Exported power
    P_exp = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="P_exp")
    # Electricity shared (that is incentivised), minimum between "P_imp" and "P_exp"
    # Shared energy according to the italian legislation
    Econd = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS,name='Econd') 
    
    ### DECISION VARIABLES for each member
    
    # res (PV, EES, HP, TES)
               
    # PV 
    CapPV_res = m.addVar(vtype=grb.GRB.CONTINUOUS, name ="CapPV_res")     
     
    # EES
    CapEES_res = m.addVar(vtype=grb.GRB.CONTINUOUS, name="CapEES_res")              # capacity
    SocEES_res = m.addVars(H+1, K,    # the SOC requires one more time step (H+1)   # State of Charge (SOC)
                           vtype=grb.GRB.CONTINUOUS, name="SocEES_res")
    PcEES_res = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="PcEES_res")         # charging power
    PdEES_res = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="PdEES_res")         # discharging power
    deltaEES_res = m.addVars(H, K, vtype=grb.GRB.BINARY, name="deltaEES_res")       # binary variable for charge/discharge
    thetacEES_res = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="thetacEES_res") # auxiliary variable for charge
    thetadEES_res = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="thetadEES_res") # auxiliary variable for discharge
     
    # HP
    CapHP_res = m.addVar(vtype=grb.GRB.CONTINUOUS, name="CapHP_res")                # capacity
    QHP_res = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="QHP_res")             # output heat
    deltaHP_res = m.addVars(H, K, vtype=grb.GRB.BINARY, name="deltaHP_res")         # binary for on-off status
    thetaHP_res = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="thetaHP_res")     # auxiliary variable
     
    # TES
    CapTES_res = m.addVar(vtype=grb.GRB.CONTINUOUS, name="CapTES_res")              # capacity
    SocTES_res = m.addVars(H+1, K,    # the SOC requires one more time step (H+1)   # SOC
                           vtype=grb.GRB.CONTINUOUS, name="SocTES_res")  
    QcTES_res = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="QcTES_res")         # charging heat
    QdTES_res = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="QdTES_res")         # discharging heat
    deltaTES_res = m.addVars(H, K, vtype=grb.GRB.BINARY, name="deltaTES_res")       # binary variable for charge/discharge
    thetacTES_res = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="thetacTES_res") # auxiliary variable for charge
    thetadTES_res = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="thetadTES_res") # auxiliary variable for discharge
     
    # Electrical power exchanged with the grid for member "res"
    P_imp_res = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="P_imp_res")
    P_exp_res = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="P_exp_res")
    
    # Flexible demand
    DemEl_res_op = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="DemEl_res_op")
    
    # ter (BOIL)
               
    # BOIL
    CapBOIL_ter = m.addVar(vtype=grb.GRB.CONTINUOUS, name="CapBOIL_ter")            # capacity
    QBOIL_ter = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="QBOIL_ter")         # output heat
             
    # Electrical power exchanged with the grid for member "ter" (only imported)
    P_imp_ter = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="P_imp_ter")
    
    # Flexible demand
    DemEl_ter_op = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="DemEl_ter_op")
        
    # com (BOIL)
               
    # BOIL
    CapBOIL_com = m.addVar(vtype=grb.GRB.CONTINUOUS, name="CapBOIL_com")            # capacity
    QBOIL_com = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="QBOIL_com")         # output heat
              
    # Electrical power exchanged with the grid for member "com" (only imported)
    P_imp_com = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="P_imp_com")
    
    # Flexible demand
    DemEl_com_op = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="DemEl_com_op")
               
    # pub (PV, EES, HP, TES)
               
    # PV 
    CapPV_pub = m.addVar(vtype=grb.GRB.CONTINUOUS, name="CapPV_pub")   
     
    # EES
    CapEES_pub = m.addVar(vtype=grb.GRB.CONTINUOUS, name="CapEES_pub")              # capacity
    SocEES_pub = m.addVars(H+1, K,    # the SOC requires one more time step (H+1)   # SOC
                           vtype=grb.GRB.CONTINUOUS, name="SocEES_pub")
    PcEES_pub = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="PcEES_pub")         # charging power
    PdEES_pub = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="PdEES_pub")         # discharging power
    deltaEES_pub = m.addVars(H, K, vtype=grb.GRB.BINARY, name="deltaEES_pub")       # binary variable for charge/discharge
    thetacEES_pub = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="thetacEES_pub") # auxiliary variable for charge
    thetadEES_pub = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="thetadEES_pub") # auxiliary variable for discharge
     
    # HP
    CapHP_pub = m.addVar(vtype=grb.GRB.CONTINUOUS, name="CapHP_pub")                # capacity
    QHP_pub = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="QHP_pub")             # output heat
    deltaHP_pub = m.addVars(H, K, vtype=grb.GRB.BINARY, name="deltaHP_pub")         # binary for on-off status
    thetaHP_pub = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="thetaHP_pub")     # auxiliary variable
     
    # TES
    CapTES_pub = m.addVar(vtype=grb.GRB.CONTINUOUS, name="CapTES_pub")              # capacity
    SocTES_pub = m.addVars(H+1, K,    # the SOC requires one more time step (H+1)   # SOC
                           vtype=grb.GRB.CONTINUOUS, name="SocTES_pub")  
    QcTES_pub = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="QcTES_pub")         # charging heat
    QdTES_pub = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="QdTES_pub")         # discharging heat
    deltaTES_pub = m.addVars(H, K, vtype=grb.GRB.BINARY, name="deltaTES_pub")       # binary variable for charge/discharge
    thetacTES_pub = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="thetacTES_pub") # auxiliary variable for charge
    thetadTES_pub = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="thetadTES_pub") # auxiliary variable for discharge
     
    # Electrical power exchanged with the grid for member "pub"
    P_imp_pub = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="P_imp_pub")
    P_exp_pub = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="P_exp_pub")
    
    # Flexible demand
    DemEl_pub_op = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="DemEl_pub_op")
    
    m.update()
    
    ### CONSTRAINTS
    
    # Ideal COP for each operating condition
    COPid = 1/(1-(Tamb+273.15)/(Tsupply+273.15)) # Carnot
    # Big M parameter (big enough threshold) to linearize constraints
    M = 1e4
    
    ### Constraints for each member
        
    # res
              
    ### PV 
    PPV_res = [[CapPV_res*SolRad[h,k] for k in range(K)] for h in range(H)] # SolRad in [kW/m^2]
    
    ### EES
    m.addConstrs(SocEES_res[h+1,k] == SocEES_res[h,k]*(1-sdEES) \
                  + PcEES_res[h,k]*np.sqrt(rteEES) \
                  - PdEES_res[h,k]/np.sqrt(rteEES)
                  for h in range(H) for k in range(K))
    # Auxiliary variable charge
    m.addConstrs(thetacEES_res[h,k] <= deltaEES_res[h,k]*M
                  for h in range(H) for k in range(K))
    m.addConstrs(CapEES_res - thetacEES_res[h,k] <= (1 - deltaEES_res[h,k])*M
                  for h in range(H) for k in range(K))
    m.addConstrs(CapEES_res - thetacEES_res[h,k] >= (1 - deltaEES_res[h,k])*0
                  for h in range(H) for k in range(K))
    # Auxiliary variable discharge
    m.addConstrs(thetadEES_res[h,k] <= (1-deltaEES_res[h,k])*M
                  for h in range(H) for k in range(K))
    m.addConstrs(CapEES_res - thetadEES_res[h,k] <= deltaEES_res[h,k]*M
                  for h in range(H) for k in range(K))
    m.addConstrs(CapEES_res - thetadEES_res[h,k] >= deltaEES_res[h,k]*0
                  for h in range(H) for k in range(K))
    # SOC upper bound
    m.addConstrs(SocEES_res[h,k] <= CapEES_res
                  for h in range(H+1) for k in range(K))
    # Input and output capaciities
    m.addConstrs(PcEES_res[h,k] <= cinpEES*thetacEES_res[h,k]
                  for h in range(H) for k in range(K))
    m.addConstrs(PdEES_res[h,k] <= coutEES*thetadEES_res[h,k]
                  for h in range(H) for k in range(K))
    # Limit operation to daily storage
    m.addConstrs(SocEES_res[0,k] == SocEES_res[H,k] for k in range(K))
    
    ### HP
    m.addConstrs(thetaHP_res[h,k] <= deltaHP_res[h,k]*M
                  for h in range(H) for k in range(K))
    m.addConstrs(CapHP_res - thetaHP_res[h,k] <= (1 - deltaHP_res[h,k])*M
                  for h in range(H) for k in range(K))
    m.addConstrs(CapHP_res - thetaHP_res[h,k] >= (1 - deltaHP_res[h,k])*0
                  for h in range(H) for k in range(K))
    m.addConstrs(QHP_res[h,k] <= thetaHP_res[h,k]
                  for h in range(H) for k in range(K))
    m.addConstrs(QHP_res[h,k] >= min_load_HP*thetaHP_res[h,k]
                  for h in range(H) for k in range(K))
    PHP_res = [[(QHP_res[h,k]*ee + deltaHP_res[h,k]*ff)/COPid[h,k] 
                 for k in range(K)] for h in range(H)]
    
    ### TES
    # SOC equation
    m.addConstrs(SocTES_res[h+1,k] == SocTES_res[h,k]*(1-sdTES) \
                  + QcTES_res[h,k]*np.sqrt(rteTES) \
                  - QdTES_res[h,k]/np.sqrt(rteTES)
                  for h in range(H) for k in range(K))
    # Auxiliary variable charge
    m.addConstrs(thetacTES_res[h,k] <= deltaTES_res[h,k]*M
                  for h in range(H) for k in range(K))
    m.addConstrs(CapTES_res - thetacTES_res[h,k] <= (1 - deltaTES_res[h,k])*M
                  for h in range(H) for k in range(K))
    m.addConstrs(CapTES_res - thetacTES_res[h,k] >= (1 - deltaTES_res[h,k])*0
                  for h in range(H) for k in range(K))
    # Auxiliary variable discharge
    m.addConstrs(thetadTES_res[h,k] <= (1-deltaTES_res[h,k])*M
                  for h in range(H) for k in range(K))
    m.addConstrs(CapTES_res - thetadTES_res[h,k] <= deltaTES_res[h,k]*M
                  for h in range(H) for k in range(K))
    m.addConstrs(CapTES_res - thetadTES_res[h,k] >= deltaTES_res[h,k]*0
                  for h in range(H) for k in range(K))
    # SOC upper bound
    m.addConstrs(SocTES_res[h,k] <= CapTES_res
                  for h in range(H+1) for k in range(K))  
    # Input and output capaciities
    m.addConstrs(QcTES_res[h,k] <= cinpTES*thetacTES_res[h,k]
                  for h in range(H) for k in range(K))
    m.addConstrs(QdTES_res[h,k] <= coutTES*thetadTES_res[h,k]
                  for h in range(H) for k in range(K))
    # Limit operation to daily storage
    m.addConstrs(SocTES_res[0,k] == SocTES_res[H,k] for k in range(K))  
    
    ### Electrical balance (P_imp_res[h,k]-P_exp_res[h,k] contributes to the shared energy)
    m.addConstrs(P_imp_res[h,k]-P_exp_res[h,k]+PPV_res[h][k] \
                  +PdEES_res[h,k]-PcEES_res[h,k]-DemEl_res_op[h,k]-PHP_res[h][k] == 0
                  for h in range(H) for k in range(K))
    ### Heating balance
    m.addConstrs(QHP_res[h,k]+QdTES_res[h,k]-QcTES_res[h,k]-DemTh_res[h,k] == 0
                  for h in range(H) for k in range(K))
    
    # Flexible demand
    m.addConstrs(sum(DemEl_res_op[h,k] for h in range(H)) == sum(DemEl_res[h,k] for h in range(H)) for k in range(K))
    m.addConstrs(DemEl_res_op[h,k] <= np.max(DemEl_res,axis=0)[k] for h in range(H) for k in range(K))           
    m.addConstrs(DemEl_res_op[h,k] >= np.min(DemEl_res,axis=0)[k] for h in range(H) for k in range(K))
    m.addConstrs(DemEl_res_op[h,k] <= (1+D_var)*DemEl_res[h,k] for h in range(H) for k in range(K))  
    m.addConstrs(DemEl_res_op[h,k] >= (1-D_var)*DemEl_res[h,k] for h in range(H) for k in range(K))
        
    # ter
              
    ### BOIL
    # Heat output upper bound
    m.addConstrs(QBOIL_ter[h,k] <= CapBOIL_ter
                  for h in range(H) for k in range(K))
    # Fuel consumption (linear expression)
    FBOIL_ter = [[QBOIL_ter[h,k]/eff_th_BOIL for k in range(K)] for h in range(H)] 
    
    ### Electrical balance (P_imp_ter[h,k] contributes to the shared energy)
    m.addConstrs(P_imp_ter[h,k]-DemEl_ter_op[h,k] == 0
                  for h in range(H) for k in range(K))
    ### Heating balance
    m.addConstrs(QBOIL_ter[h,k]-DemTh_ter[h,k] == 0
                  for h in range(H) for k in range(K))
    
    # Flexible demand
    m.addConstrs(sum(DemEl_ter_op[h,k] for h in range(H)) == sum(DemEl_ter[h,k] for h in range(H)) for k in range(K))
    m.addConstrs(DemEl_ter_op[h,k] <= np.max(DemEl_ter,axis=0)[k] for h in range(H) for k in range(K))           
    m.addConstrs(DemEl_ter_op[h,k] >= np.min(DemEl_ter,axis=0)[k] for h in range(H) for k in range(K))
    m.addConstrs(DemEl_ter_op[h,k] <= (1+D_var)*DemEl_ter[h,k] for h in range(H) for k in range(K))  
    m.addConstrs(DemEl_ter_op[h,k] >= (1-D_var)*DemEl_ter[h,k] for h in range(H) for k in range(K))
        
    # com
              
    ### BOIL
    # Heat output upper bound
    m.addConstrs(QBOIL_com[h,k] <= CapBOIL_com
                  for h in range(H) for k in range(K))
    # Fuel consumption (linear expression)
    FBOIL_com = [[QBOIL_com[h,k]/eff_th_BOIL for k in range(K)] for h in range(H)]
    
    ### Electrical balance (P_imp_com[h,k] contributes to the shared energy)
    m.addConstrs(P_imp_com[h,k]-DemEl_com_op[h,k] == 0
                  for h in range(H) for k in range(K))
    ### Heating balance
    m.addConstrs(QBOIL_com[h,k]-DemTh_com[h,k] == 0
                  for h in range(H) for k in range(K))
    
    # Flexible demand
    m.addConstrs(sum(DemEl_com_op[h,k] for h in range(H)) == sum(DemEl_com[h,k] for h in range(H)) for k in range(K))
    m.addConstrs(DemEl_com_op[h,k] <= np.max(DemEl_com,axis=0)[k] for h in range(H) for k in range(K))           
    m.addConstrs(DemEl_com_op[h,k] >= np.min(DemEl_com,axis=0)[k] for h in range(H) for k in range(K))
    m.addConstrs(DemEl_com_op[h,k] <= (1+D_var)*DemEl_com[h,k] for h in range(H) for k in range(K))  
    m.addConstrs(DemEl_com_op[h,k] >= (1-D_var)*DemEl_com[h,k] for h in range(H) for k in range(K))
        
    # pub
              
    ### PV 
    PPV_pub = [[CapPV_pub*SolRad[h,k] for k in range(K)] for h in range(H)] # SolRad in [kW/m^2]
        
    ### EES
    m.addConstrs(SocEES_pub[h+1,k] == SocEES_pub[h,k]*(1-sdEES) \
                  + PcEES_pub[h,k]*np.sqrt(rteEES) \
                  - PdEES_pub[h,k]/np.sqrt(rteEES)
                  for h in range(H) for k in range(K))
    # Auxiliary variable charge
    m.addConstrs(thetacEES_pub[h,k] <= deltaEES_pub[h,k]*M
                  for h in range(H) for k in range(K))
    m.addConstrs(CapEES_pub - thetacEES_pub[h,k] <= (1 - deltaEES_pub[h,k])*M
                  for h in range(H) for k in range(K))
    m.addConstrs(CapEES_pub - thetacEES_pub[h,k] >= (1 - deltaEES_pub[h,k])*0
                  for h in range(H) for k in range(K))
    # Auxiliary variable discharge
    m.addConstrs(thetadEES_pub[h,k] <= (1-deltaEES_pub[h,k])*M
                  for h in range(H) for k in range(K))
    m.addConstrs(CapEES_pub - thetadEES_pub[h,k] <= deltaEES_pub[h,k]*M
                  for h in range(H) for k in range(K))
    m.addConstrs(CapEES_pub - thetadEES_pub[h,k] >= deltaEES_pub[h,k]*0
                  for h in range(H) for k in range(K))
    # SOC upper bound
    m.addConstrs(SocEES_pub[h,k] <= CapEES_pub
                  for h in range(H+1) for k in range(K))
    # Input and output capaciities
    m.addConstrs(PcEES_pub[h,k] <= cinpEES*thetacEES_pub[h,k]
                  for h in range(H) for k in range(K))
    m.addConstrs(PdEES_pub[h,k] <= coutEES*thetadEES_pub[h,k]
                  for h in range(H) for k in range(K))
    # Limit operation to daily storage
    m.addConstrs(SocEES_pub[0,k] == SocEES_pub[H,k] for k in range(K))
    
    ### HP
    m.addConstrs(thetaHP_pub[h,k] <= deltaHP_pub[h,k]*M
                  for h in range(H) for k in range(K))
    m.addConstrs(CapHP_pub - thetaHP_pub[h,k] <= (1 - deltaHP_pub[h,k])*M
                  for h in range(H) for k in range(K))
    m.addConstrs(CapHP_pub - thetaHP_pub[h,k] >= (1 - deltaHP_pub[h,k])*0
                  for h in range(H) for k in range(K))
    m.addConstrs(QHP_pub[h,k] <= thetaHP_pub[h,k]
                  for h in range(H) for k in range(K))
    m.addConstrs(QHP_pub[h,k] >= min_load_HP*thetaHP_pub[h,k]
                  for h in range(H) for k in range(K))
    PHP_pub = [[(QHP_pub[h,k]*ee + deltaHP_pub[h,k]*ff)/COPid[h,k] 
                 for k in range(K)] for h in range(H)]
    
    ### TES
    # SOC equation
    m.addConstrs(SocTES_pub[h+1,k] == SocTES_pub[h,k]*(1-sdTES) \
                  + QcTES_pub[h,k]*np.sqrt(rteTES) \
                  - QdTES_pub[h,k]/np.sqrt(rteTES)
                  for h in range(H) for k in range(K))
    # Auxiliary variable charge
    m.addConstrs(thetacTES_pub[h,k] <= deltaTES_pub[h,k]*M
                  for h in range(H) for k in range(K))
    m.addConstrs(CapTES_pub - thetacTES_pub[h,k] <= (1 - deltaTES_pub[h,k])*M
                  for h in range(H) for k in range(K))
    m.addConstrs(CapTES_pub - thetacTES_pub[h,k] >= (1 - deltaTES_pub[h,k])*0
                  for h in range(H) for k in range(K))
    # Auxiliary variable discharge
    m.addConstrs(thetadTES_pub[h,k] <= (1-deltaTES_pub[h,k])*M
                  for h in range(H) for k in range(K))
    m.addConstrs(CapTES_pub - thetadTES_pub[h,k] <= deltaTES_pub[h,k]*M
                  for h in range(H) for k in range(K))
    m.addConstrs(CapTES_pub - thetadTES_pub[h,k] >= deltaTES_pub[h,k]*0
                  for h in range(H) for k in range(K))
    # SOC upper bound
    m.addConstrs(SocTES_pub[h,k] <= CapTES_pub
                  for h in range(H+1) for k in range(K))  
    # Input and output capaciities
    m.addConstrs(QcTES_pub[h,k] <= cinpTES*thetacTES_pub[h,k]
                  for h in range(H) for k in range(K))
    m.addConstrs(QdTES_pub[h,k] <= coutTES*thetadTES_pub[h,k]
                  for h in range(H) for k in range(K))
    # Limit operation to daily storage
    m.addConstrs(SocTES_pub[0,k] == SocTES_pub[H,k] for k in range(K))  
    
    ### Electrical balance (P_imp_pub[h,k]-P_exp_pub[h,k] contributes to the shared energy)
    m.addConstrs(P_imp_pub[h,k]-P_exp_pub[h,k]+PPV_pub[h][k] \
                  +PdEES_pub[h,k]-PcEES_pub[h,k]-DemEl_pub_op[h,k]-PHP_pub[h][k] == 0
                  for h in range(H) for k in range(K))
        
    ### Heating balance
    m.addConstrs(QHP_pub[h,k]+QdTES_pub[h,k]-QcTES_pub[h,k]-DemTh_pub[h,k] == 0
                  for h in range(H) for k in range(K))
    
    # Flexible demand
    m.addConstrs(sum(DemEl_pub_op[h,k] for h in range(H)) == sum(DemEl_pub[h,k] for h in range(H)) for k in range(K))
    m.addConstrs(DemEl_pub_op[h,k] <= np.max(DemEl_pub,axis=0)[k] for h in range(H) for k in range(K))           
    m.addConstrs(DemEl_pub_op[h,k] >= np.min(DemEl_pub,axis=0)[k] for h in range(H) for k in range(K))
    m.addConstrs(DemEl_pub_op[h,k] <= (1+D_var)*DemEl_pub[h,k] for h in range(H) for k in range(K))  
    m.addConstrs(DemEl_pub_op[h,k] >= (1-D_var)*DemEl_pub[h,k] for h in range(H) for k in range(K))
    
    # Total net electrical power imported from (P_imp) and exported to (P_exp) the grid 
    m.addConstrs(P_imp[h,k] == P_imp_res[h,k]+P_imp_ter[h,k]+P_imp_com[h,k]+P_imp_pub[h,k] for h in range(H) for k in range(K))
    m.addConstrs(P_exp[h,k] == P_exp_res[h,k]+P_exp_pub[h,k] for h in range(H) for k in range(K))
    
    ### Calculation of electricity shared
    m.addConstrs(Econd[h,k] == min_(P_imp[h,k],P_exp[h,k]) for h in range(H) for k in range(K))
    
    m.update()
    
    ### DEFINE COSTS
    
    ### The solution is based on 4 average seasonal days (without uncertainty)
    ### or on 2 to 7 typical days for each season (with uncertainty)
    ### Time horizon of 1 year    
    ### The investment cost is actualized on 1 year 
    ### The yearly operational cost is weighted by the scenario probabilities
    
    # Investment cost [€]
    CostINV = (r*np.power(1+r,nPV)/(np.power(1+r,nPV)-1)+OnM_fix_PV)*invPV*(CapPV_res+CapPV_pub) +\
              (r*np.power(1+r,nBOIL)/(np.power(1+r,nBOIL)-1)+OnM_fix_BOIL)*invBOIL*(CapBOIL_ter+CapBOIL_com) +\
              (r*np.power(1+r,nHP)/(np.power(1+r,nHP)-1)+OnM_fix_HP)*invHP*(CapHP_res+CapHP_pub) +\
              (r*np.power(1+r,nTES)/(np.power(1+r,nTES)-1)+OnM_fix_TES)*invTES*(CapTES_res+CapTES_pub) +\
              (r*np.power(1+r,nEES)/(np.power(1+r,nEES)-1)+OnM_fix_EES)*invEES*(CapEES_res+CapEES_pub)
    
    # Operational costs [€], including the REC incentive
    CostOP = sum(w_scenarios[k]*sum((FBOIL_ter[h][k]+FBOIL_com[h][k])*cbuy_ng 
                                         + P_imp[h,k]*cbuy_el[h,k] - P_exp[h,k]*csell_el[h,k] - Econd[h,k]*inc
                                         for h in range(H)) for k in range(K))
    
    m.update()
        
    ### SET THE OBJECTIVE FUNCTION
    
    m.setObjective(CostINV+CostOP, grb.GRB.MINIMIZE)    
    
    def data_cb(m, where):
        if where == grb.GRB.Callback.MIPSOL:
           cur_obj = m.cbGet(grb.GRB.Callback.MIPSOL_OBJBST)  # best incumbent solution (upper bound, primal solution)
           cur_bd = m.cbGet(grb.GRB.Callback.MIPSOL_OBJBND)   # best bound solution (lower bound, dual solution)
           opt_gap.append(((cur_obj-cur_bd)/cur_obj)*100)     # all gaps during iterations
    
    opt_gap = []
    m.optimize(callback=data_cb)   # optimization with callback to retrieve optimal gaps
    opt_gap = np.array(opt_gap)
    opt_gap = opt_gap[-1]          # optimal gap
    print(opt_gap)

    # Assess the objectives
    C_TOT = m.objVal  # optimal value of the objective function [€]
    design_cost = CostINV.getValue()
    operation_cost = CostOP.getValue()
    
    for v in m.getVars():
        array.append(v.x)

    ### RESULTS ###
    
    ### Capacities
    # res
    cap_PV_res = CapPV_res.X
    cap_EES_res = CapEES_res.X
    cap_HP_res = CapHP_res.X
    cap_TES_res = CapTES_res.X
    # ter
    cap_BOIL_ter = CapBOIL_ter.X
    # com
    cap_BOIL_com = CapBOIL_com.X
    # pub
    cap_PV_pub = CapPV_pub.X
    cap_EES_pub = CapEES_pub.X
    cap_HP_pub = CapHP_pub.X
    cap_TES_pub = CapTES_pub.X
    
    ### Operational variables
    array_P_imp=[]
    array_P_exp=[]
    array_Econd=[]
    # res
    array_SocEES_res=[]
    array_PcEES_res=[]
    array_PdEES_res=[]
    array_deltaEES_res=[]
    array_QHP_res=[]
    array_deltaHP_res=[]
    array_SocTES_res=[]
    array_QcTES_res=[]
    array_QdTES_res=[]
    array_deltaTES_res=[]
    array_P_imp_res=[]
    array_P_exp_res=[]
    # ter
    array_QBOIL_ter=[]
    array_P_imp_ter=[]
    # com
    array_QBOIL_com=[]
    array_P_imp_com=[]
    # pub
    array_SocEES_pub=[]
    array_PcEES_pub=[]
    array_PdEES_pub=[]
    array_deltaEES_pub=[]
    array_QHP_pub=[]
    array_deltaHP_pub=[]
    array_SocTES_pub=[]
    array_QcTES_pub=[]
    array_QdTES_pub=[]
    array_deltaTES_pub=[]
    array_P_imp_pub=[]
    array_P_exp_pub=[]
    # DR
    array_DemEl_res_op=[]
    array_DemEl_ter_op=[]
    array_DemEl_com_op=[]
    array_DemEl_pub_op=[]
    
    
    for h in range(H): 
        for k in range(K):
            array_P_imp.append(P_imp[h,k].X)
            array_P_exp.append(P_exp[h,k].X)
            array_Econd.append(Econd[h,k].X)
    
    for h in range(H+1): 
        for k in range(K):
            array_SocEES_res.append(SocEES_res[h,k].X)
            array_SocTES_res.append(SocTES_res[h,k].X)
            array_SocEES_pub.append(SocEES_pub[h,k].X)
            array_SocTES_pub.append(SocTES_pub[h,k].X)
    
    for h in range(H): 
        for k in range(K):
            # res 
            array_PcEES_res.append(PcEES_res[h,k].X)
            array_PdEES_res.append(PdEES_res[h,k].X)
            array_deltaEES_res.append(deltaEES_res[h,k].X)
            array_QHP_res.append(QHP_res[h,k].X)
            array_deltaHP_res.append(deltaHP_res[h,k].X)
            array_QcTES_res.append(QcTES_res[h,k].X)
            array_QdTES_res.append(QdTES_res[h,k].X)
            array_deltaTES_res.append(deltaTES_res[h,k].X)
            array_P_imp_res.append(P_imp_res[h,k].X)
            array_P_exp_res.append(P_exp_res[h,k].X)
            # ter
            array_QBOIL_ter.append(QBOIL_ter[h,k].X)
            array_P_imp_ter.append(P_imp_ter[h,k].X)
            # com
            array_QBOIL_com.append(QBOIL_com[h,k].X)
            array_P_imp_com.append(P_imp_com[h,k].X)
            # pub
            array_PcEES_pub.append(PcEES_pub[h,k].X)
            array_PdEES_pub.append(PdEES_pub[h,k].X)
            array_deltaEES_pub.append(deltaEES_pub[h,k].X)
            array_QHP_pub.append(QHP_pub[h,k].X)
            array_deltaHP_pub.append(deltaHP_pub[h,k].X)
            array_QcTES_pub.append(QcTES_pub[h,k].X)
            array_QdTES_pub.append(QdTES_pub[h,k].X)
            array_deltaTES_pub.append(deltaTES_pub[h,k].X)
            array_P_imp_pub.append(P_imp_pub[h,k].X)
            array_P_exp_pub.append(P_exp_pub[h,k].X)
            # DR
            array_DemEl_res_op.append(DemEl_res_op[h,k].X)
            array_DemEl_ter_op.append(DemEl_ter_op[h,k].X)
            array_DemEl_com_op.append(DemEl_com_op[h,k].X)
            array_DemEl_pub_op.append(DemEl_pub_op[h,k].X)
           
    array_P_imp=np.array(array_P_imp)
    array_P_exp=np.array(array_P_exp)
    array_Econd=np.array(array_Econd)
    # res
    array_SocEES_res=np.array(array_SocEES_res)
    array_PcEES_res=np.array(array_PcEES_res)
    array_PdEES_res=np.array(array_PdEES_res)
    array_deltaEES_res=np.array(array_deltaEES_res)
    array_QHP_res=np.array(array_QHP_res)
    array_deltaHP_res=np.array(array_deltaHP_res)
    array_SocTES_res=np.array(array_SocTES_res)
    array_QcTES_res=np.array(array_QcTES_res)
    array_QdTES_res=np.array(array_QdTES_res)
    array_deltaTES_res=np.array(array_deltaTES_res)
    array_P_imp_res=np.array(array_P_imp_res)
    array_P_exp_res=np.array(array_P_exp_res)
    # ter
    array_QBOIL_ter=np.array(array_QBOIL_ter)
    array_P_imp_ter=np.array(array_P_imp_ter)
    # com
    array_QBOIL_com=np.array(array_QBOIL_com)
    array_P_imp_com=np.array(array_P_imp_com)
    # pub
    array_SocEES_pub=np.array(array_SocEES_pub)
    array_PcEES_pub=np.array(array_PcEES_pub)
    array_PdEES_pub=np.array(array_PdEES_pub)
    array_deltaEES_pub=np.array(array_deltaEES_pub)
    array_QHP_pub=np.array(array_QHP_pub)
    array_deltaHP_pub=np.array(array_deltaHP_pub)
    array_SocTES_pub=np.array(array_SocTES_pub)
    array_QcTES_pub=np.array(array_QcTES_pub)
    array_QdTES_pub=np.array(array_QdTES_pub)
    array_deltaTES_pub=np.array(array_deltaTES_pub)
    array_P_imp_pub=np.array(array_P_imp_pub)
    array_P_exp_pub=np.array(array_P_exp_pub)
    # DR
    array_DemEl_res_op=np.array(array_DemEl_res_op)
    array_DemEl_ter_op=np.array(array_DemEl_ter_op)
    array_DemEl_com_op=np.array(array_DemEl_com_op)
    array_DemEl_pub_op=np.array(array_DemEl_pub_op)
    
    # reshape of arrays
    array_P_imp=np.reshape(array_P_imp,(H,K),order='C')
    array_P_exp=np.reshape(array_P_exp,(H,K),order='C')
    array_Econd=np.reshape(array_Econd,(H,K),order='C')
    # res
    array_SocEES_res=np.reshape(array_SocEES_res,(H+1,K),order='C')
    array_PcEES_res=np.reshape(array_PcEES_res,(H,K),order='C')
    array_PdEES_res=np.reshape(array_PdEES_res,(H,K),order='C')
    array_deltaEES_res=np.reshape(array_deltaEES_res,(H,K),order='C')
    array_QHP_res=np.reshape(array_QHP_res,(H,K),order='C')
    array_deltaHP_res=np.reshape(array_deltaHP_res,(H,K),order='C')
    array_SocTES_res=np.reshape(array_SocTES_res,(H+1,K),order='C')
    array_QcTES_res=np.reshape(array_QcTES_res,(H,K),order='C')
    array_QdTES_res=np.reshape(array_QdTES_res,(H,K),order='C')
    array_deltaTES_res=np.reshape(array_deltaTES_res,(H,K),order='C')
    array_P_imp_res=np.reshape(array_P_imp_res,(H,K),order='C')
    array_P_exp_res=np.reshape(array_P_exp_res,(H,K),order='C')
    # ter
    array_QBOIL_ter=np.reshape(array_QBOIL_ter,(H,K),order='C')
    array_P_imp_ter=np.reshape(array_P_imp_ter,(H,K),order='C')
    # com
    array_QBOIL_com=np.reshape(array_QBOIL_com,(H,K),order='C')
    array_P_imp_com=np.reshape(array_P_imp_com,(H,K),order='C')
    # pub
    array_SocEES_pub=np.reshape(array_SocEES_pub,(H+1,K),order='C')
    array_PcEES_pub=np.reshape(array_PcEES_pub,(H,K),order='C')
    array_PdEES_pub=np.reshape(array_PdEES_pub,(H,K),order='C')
    array_deltaEES_pub=np.reshape(array_deltaEES_pub,(H,K),order='C')
    array_QHP_pub=np.reshape(array_QHP_pub,(H,K),order='C')
    array_deltaHP_pub=np.reshape(array_deltaHP_pub,(H,K),order='C')
    array_SocTES_pub=np.reshape(array_SocTES_pub,(H+1,K),order='C')
    array_QcTES_pub=np.reshape(array_QcTES_pub,(H,K),order='C')
    array_QdTES_pub=np.reshape(array_QdTES_pub,(H,K),order='C')
    array_deltaTES_pub=np.reshape(array_deltaTES_pub,(H,K),order='C')
    array_P_imp_pub=np.reshape(array_P_imp_pub,(H,K),order='C')
    array_P_exp_pub=np.reshape(array_P_exp_pub,(H,K),order='C')
    # DR
    array_DemEl_res_op=np.reshape(array_DemEl_res_op,(H,K),order='C')
    array_DemEl_ter_op=np.reshape(array_DemEl_ter_op,(H,K),order='C')
    array_DemEl_com_op=np.reshape(array_DemEl_com_op,(H,K),order='C')
    array_DemEl_pub_op=np.reshape(array_DemEl_pub_op,(H,K),order='C')
    
    # Output
    return  C_TOT, design_cost, operation_cost, array_P_imp, array_P_exp, array_Econd, \
            cap_PV_res, cap_EES_res, cap_HP_res, cap_TES_res, array_SocTES_res, array_QcTES_res, array_QdTES_res, \
            cap_BOIL_ter, cap_BOIL_com, \
            cap_PV_pub, cap_EES_pub, cap_HP_pub, cap_TES_pub, opt_gap, csell_el, cbuy_el

#%% AVERAGE SEASONAL DAYS AND STOCHASTIC SCENARIOS

# Length of a day
H = 24

print('Choice of the city')
address=input()  # input location
PATH_clustering = "C:/Users/utente/Desktop/Università Gabriele/Assegno_biennale/Models/Model_ML2/Input_data/Timeseries/"+address
timeseries_w = np.load(os.path.join(PATH_clustering,'Clustering_winter_'+address+'.npz'))
timeseries_sp = np.load(os.path.join(PATH_clustering,'Clustering_spring_'+address+'.npz'))
timeseries_su = np.load(os.path.join(PATH_clustering,'Clustering_summer_'+address+'.npz'))
timeseries_a = np.load(os.path.join(PATH_clustering,'Clustering_autumn_'+address+'.npz'))

# Number of years of the whole dataset (2005-2020)
years = timeseries_w['years']

# Average seasonal days and their weights
DF_winter_mean = timeseries_w['DF_winter_mean']
DF_spring_mean = timeseries_sp['DF_spring_mean']
DF_summer_mean = timeseries_su['DF_summer_mean']
DF_autumn_mean = timeseries_a['DF_autumn_mean']
N_days_winter = timeseries_w['N_days_winter']
N_days_spring = timeseries_sp['N_days_spring']
N_days_summer = timeseries_su['N_days_summer']
N_days_autumn = timeseries_a['N_days_autumn']
# Average seasonal days as (H,4,years-6) aggregating the different seasons in the training dataset (2005-2014)
SolRad_av = np.zeros((H*4,years-6))       # solar irradiance
Tamb_av = np.zeros((H*4,years-6))         # ambient temperature
DemEl_res_av = np.zeros((H*4,years-6))    # electricity demand of Res
DemEl_ter_av = np.zeros((H*4,years-6))    # electricity demand of Ter
DemEl_com_av = np.zeros((H*4,years-6))    # electricity demand of Com
DemEl_pub_av = np.zeros((H*4,years-6))    # electricity demand of Pub
DemTh_res_av = np.zeros((H*4,years-6))    # heating demand of Res 
DemTh_ter_av = np.zeros((H*4,years-6))    # heating demand of Ter
DemTh_com_av = np.zeros((H*4,years-6))    # heating demand of Com
DemTh_pub_av = np.zeros((H*4,years-6))    # heating demand of Pub
csell_el_av = np.zeros((H*4,years-6))     # PUN sale prices
for year in range(2005,2015):    # different years of the training dataset
    SolRad_av[:,year-2005] = np.concatenate([DF_winter_mean[year-2005,:,3],DF_spring_mean[year-2005,:,3],DF_summer_mean[year-2005,:,3],DF_autumn_mean[year-2005,:,3]])/1000 # [kW/m^2]
    Tamb_av[:,year-2005] = np.concatenate([DF_winter_mean[year-2005,:,4],DF_spring_mean[year-2005,:,4],DF_summer_mean[year-2005,:,4],DF_autumn_mean[year-2005,:,4]])
    DemEl_res_av[:,year-2005] = np.concatenate([DF_winter_mean[year-2005,:,6],DF_spring_mean[year-2005,:,6],DF_summer_mean[year-2005,:,6],DF_autumn_mean[year-2005,:,6]])
    DemEl_ter_av[:,year-2005] = np.concatenate([DF_winter_mean[year-2005,:,7],DF_spring_mean[year-2005,:,7],DF_summer_mean[year-2005,:,7],DF_autumn_mean[year-2005,:,7]])
    DemEl_com_av[:,year-2005] = np.concatenate([DF_winter_mean[year-2005,:,8],DF_spring_mean[year-2005,:,8],DF_summer_mean[year-2005,:,8],DF_autumn_mean[year-2005,:,8]])
    DemEl_pub_av[:,year-2005] = np.concatenate([DF_winter_mean[year-2005,:,9],DF_spring_mean[year-2005,:,9],DF_summer_mean[year-2005,:,9],DF_autumn_mean[year-2005,:,9]])
    DemTh_res_av[:,year-2005] = np.concatenate([DF_winter_mean[year-2005,:,10],DF_spring_mean[year-2005,:,10],DF_summer_mean[year-2005,:,10],DF_autumn_mean[year-2005,:,10]])
    DemTh_ter_av[:,year-2005] = np.concatenate([DF_winter_mean[year-2005,:,11],DF_spring_mean[year-2005,:,11],DF_summer_mean[year-2005,:,11],DF_autumn_mean[year-2005,:,11]])
    DemTh_com_av[:,year-2005] = np.concatenate([DF_winter_mean[year-2005,:,12],DF_spring_mean[year-2005,:,12],DF_summer_mean[year-2005,:,12],DF_autumn_mean[year-2005,:,12]])
    DemTh_pub_av[:,year-2005] = np.concatenate([DF_winter_mean[year-2005,:,13],DF_spring_mean[year-2005,:,13],DF_summer_mean[year-2005,:,13],DF_autumn_mean[year-2005,:,13]])
    csell_el_av[:,year-2005] = np.concatenate([DF_winter_mean[year-2005,:,14],DF_spring_mean[year-2005,:,14],DF_summer_mean[year-2005,:,14],DF_autumn_mean[year-2005,:,14]])/1000*0.45 # [€/kWh]
SolRad_av=np.reshape(SolRad_av,(H,4,years-6),order='F')
Tamb_av=np.reshape(Tamb_av,(H,4,years-6),order='F')
DemEl_res_av=np.reshape(DemEl_res_av,(H,4,years-6),order='F')
DemEl_ter_av=np.reshape(DemEl_ter_av,(H,4,years-6),order='F')
DemEl_com_av=np.reshape(DemEl_com_av,(H,4,years-6),order='F')
DemEl_pub_av=np.reshape(DemEl_pub_av,(H,4,years-6),order='F')
DemTh_res_av=np.reshape(DemTh_res_av,(H,4,years-6),order='F')
DemTh_ter_av=np.reshape(DemTh_ter_av,(H,4,years-6),order='F')
DemTh_com_av=np.reshape(DemTh_com_av,(H,4,years-6),order='F')
DemTh_pub_av=np.reshape(DemTh_pub_av,(H,4,years-6),order='F')
csell_el_av=np.reshape(csell_el_av,(H,4,years-6),order='F')
cbuy_el_av = csell_el_av + 0.20 # purchase price, 0.15 €/kWh is the grid tariff

### Stochastic scenarios of the uncertain parameters are found by clustering each season of each year of the training dataset (2005-2014)
### The typical days found by clustering are stochastic scenarios
### The other timeseries (electricity and heating demands, PUN prices) are associated with the typical days of the uncertain parameters

# Different numbers of typical days assessed
n_cases = timeseries_w['n_cases']
# Number of columns for the complete dataset 
lenghth_features_c = timeseries_w['lenghth_features_c']
# Training dataset with all timeseries
train_test_c_w = timeseries_w['train_test_c_w']
train_test_c_sp = timeseries_sp['train_test_c_sp']
train_test_c_su = timeseries_su['train_test_c_su']
train_test_c_a = timeseries_a['train_test_c_a']
# Indexes of typical days
index_representative_r_w = timeseries_w['index_representative_r']
index_representative_r_sp = timeseries_sp['index_representative_r']
index_representative_r_su = timeseries_su['index_representative_r']
index_representative_r_a = timeseries_a['index_representative_r']
# Weights of the typical days
number_elements_each_cluster_w = timeseries_w['number_elements_each_cluster']
number_elements_each_cluster_sp = timeseries_sp['number_elements_each_cluster']
number_elements_each_cluster_su = timeseries_su['number_elements_each_cluster']
number_elements_each_cluster_a = timeseries_a['number_elements_each_cluster']

# Typical days with all timeseries, varying the number of typical days for each season of each year of the training dataset (10 years)
typical_days_w = np.zeros((n_cases+1,lenghth_features_c,n_cases,years-6))
typical_days_sp = np.zeros((n_cases+1,lenghth_features_c,n_cases,years-6))
typical_days_su = np.zeros((n_cases+1,lenghth_features_c,n_cases,years-6))
typical_days_a = np.zeros((n_cases+1,lenghth_features_c,n_cases,years-6))
# Typical days and weights aggregating the different seasons
typical_days = np.zeros(((n_cases+1)*4,lenghth_features_c,n_cases,years-6))
number_elements_each_cluster = np.zeros((n_cases,(n_cases+1)*4,years-6))
for year in range(2005,2015):    # different years of the training dataset
    for c in range(2,n_cases+2): # different numbers of clusters
        for yi in range(c):      # different clusters
            typical_days_w[yi,:,c-2,year-2005] = train_test_c_w[index_representative_r_w[c-2,yi,year-2005],:,year-2005]
            typical_days_sp[yi,:,c-2,year-2005] = train_test_c_sp[index_representative_r_sp[c-2,yi,year-2005],:,year-2005]
            typical_days_su[yi,:,c-2,year-2005] = train_test_c_su[index_representative_r_su[c-2,yi,year-2005],:,year-2005]
            typical_days_a[yi,:,c-2,year-2005] = train_test_c_a[index_representative_r_a[c-2,yi,year-2005],:,year-2005]
        typical_days[0:c*4,:,c-2,year-2005] = np.concatenate([typical_days_w[0:c,:,c-2,year-2005],
                                                              typical_days_sp[0:c,:,c-2,year-2005],
                                                              typical_days_su[0:c,:,c-2,year-2005],
                                                              typical_days_a[0:c,:,c-2,year-2005]])
        number_elements_each_cluster[c-2,0:c*4,year-2005] = np.concatenate([number_elements_each_cluster_w[c-2,0:c,year-2005],
                                                                            number_elements_each_cluster_sp[c-2,0:c,year-2005],
                                                                            number_elements_each_cluster_su[c-2,0:c,year-2005],
                                                                            number_elements_each_cluster_a[c-2,0:c,year-2005]])
# To check the aggregation in a specific year
A = typical_days_w[:,:,:,0]
B = typical_days_sp[:,:,:,0]
C = typical_days_su[:,:,:,0]
D = typical_days_a[:,:,:,0]
E = typical_days[:,:,:,0]

#%% SOLVING THE MODEL WITHOUT UNCERTAINTY
### Solving the model in each year without uncertainty: 4 average seasonal days

start_time = time.time()

# Results for 4 average seasonal days for each year of the training dataset
# Optimal costs
C_TOT_d = np.zeros(years-6)
design_cost_d = np.zeros(years-6)
operation_cost_d = np.zeros(years-6)
# Optimal operational variables
P_imp_d = np.zeros((H,4,years-6))
P_exp_d = np.zeros((H,4,years-6))
E_cond_d = np.zeros((H,4,years-6))
SocTES_res_d = np.zeros((H+1,4,years-6))
QcTES_res_d = np.zeros((H,4,years-6))
QdTES_res_d = np.zeros((H,4,years-6))
# Optimal capacities
capPV_res_d = np.zeros(years-6)
capEES_res_d = np.zeros(years-6)
capHP_res_d = np.zeros(years-6)
capTES_res_d = np.zeros(years-6)
capBOIL_ter_d = np.zeros(years-6)
capBOIL_com_d = np.zeros(years-6)
capPV_pub_d = np.zeros(years-6)
capEES_pub_d = np.zeros(years-6)
capHP_pub_d = np.zeros(years-6)
capTES_pub_d = np.zeros(years-6)
# Optimal gaps
opt_gap_d = np.zeros(years-6)
# Sale/purchase prices 
c_sell_save_d = np.zeros((H,4,years-6)) 
c_buy_save_d = np.zeros((H,4,years-6)) 

# The input timeseries are defined for each year of the training dataset (from 2005 to 2014)
# and each average seasonal day
SolRad = np.zeros((H,4))       # solar irradiance
Tamb = np.zeros((H,4))         # ambient temperature
DemEl_res = np.zeros((H,4))    # electricity demand of Res
DemEl_ter = np.zeros((H,4))    # electricity demand of Ter
DemEl_com = np.zeros((H,4))    # electricity demand of Com
DemEl_pub = np.zeros((H,4))    # electricity demand of Pub
DemTh_res = np.zeros((H,4))    # heating demand of Res 
DemTh_ter = np.zeros((H,4))    # heating demand of Ter
DemTh_com = np.zeros((H,4))    # heating demand of Com
DemTh_pub = np.zeros((H,4))    # heating demand of Pub
csell_el = np.zeros((H,4))     # PUN sale prices

for year in range(2005,2015):          # different years of the training dataset
    SolRad = SolRad_av[:,:,year-2005]
    Tamb = Tamb_av[:,:,year-2005]
    DemEl_res = DemEl_res_av[:,:,year-2005]
    DemEl_ter = DemEl_ter_av[:,:,year-2005]
    DemEl_com = DemEl_com_av[:,:,year-2005]
    DemEl_pub = DemEl_pub_av[:,:,year-2005]
    DemTh_res = DemTh_res_av[:,:,year-2005]
    DemTh_ter = DemTh_ter_av[:,:,year-2005]
    DemTh_com = DemTh_com_av[:,:,year-2005]
    DemTh_pub = DemTh_pub_av[:,:,year-2005]
    csell_el = csell_el_av[:,:,year-2005]
    cbuy_el = cbuy_el_av[:,:,year-2005]
    w_scenarios = np.array([N_days_winter,N_days_spring,N_days_summer,N_days_autumn])
    
    C_TOT_d[year-2005], design_cost_d[year-2005], operation_cost_d[year-2005], \
    P_imp_d[:,:,year-2005], P_exp_d[:,:,year-2005], E_cond_d[:,:,year-2005], \
    capPV_res_d[year-2005], capEES_res_d[year-2005], capHP_res_d[year-2005], capTES_res_d[year-2005], SocTES_res_d[:,:,year-2005], QcTES_res_d[:,:,year-2005], QdTES_res_d[:,:,year-2005], \
    capBOIL_ter_d[year-2005], capBOIL_com_d[year-2005], \
    capPV_pub_d[year-2005], capEES_pub_d[year-2005], capHP_pub_d[year-2005], capTES_pub_d[year-2005], opt_gap_d[year-2005], \
    c_sell_save_d[:,:,year-2005], c_buy_save_d[:,:,year-2005] = REC_MES_typical_day(4)

### TIME TO SOLVE THE MODEL ###
print("--- %s seconds ---" % (time.time() - start_time))
    
#%% SOLVING THE MODEL WITH UNCERTAINTY
### Solving the model in each year with uncertainty: from 2 to 7 stochastic scenarios for each season

start_time = time.time()

# Results for each case (different numbers of typical days for each season) and each year of the training dataset
# Optimal costs
C_TOT_sto_c = np.zeros((n_cases,years-6))
design_cost_c = np.zeros((n_cases,years-6))
operation_cost_c = np.zeros((n_cases,years-6))
# Optimal operational variables
P_imp_c = np.zeros((H,(n_cases+1)*4,n_cases,years-6))
P_exp_c = np.zeros((H,(n_cases+1)*4,n_cases,years-6))
E_cond_c = np.zeros((H,(n_cases+1)*4,n_cases,years-6))
SocTES_res_c = np.zeros((H+1,(n_cases+1)*4,n_cases,years-6))
QcTES_res_c = np.zeros((H,(n_cases+1)*4,n_cases,years-6))
QdTES_res_c = np.zeros((H,(n_cases+1)*4,n_cases,years-6))
# Optimal capacities
capPV_res_c = np.zeros((n_cases,years-6))
capEES_res_c = np.zeros((n_cases,years-6))
capHP_res_c = np.zeros((n_cases,years-6))
capTES_res_c = np.zeros((n_cases,years-6))
capBOIL_ter_c = np.zeros((n_cases,years-6))
capBOIL_com_c = np.zeros((n_cases,years-6))
capPV_pub_c = np.zeros((n_cases,years-6))
capEES_pub_c = np.zeros((n_cases,years-6))
capHP_pub_c = np.zeros((n_cases,years-6))
capTES_pub_c = np.zeros((n_cases,years-6))
# Optimal gaps
opt_gap_c = np.zeros((n_cases,years-6))
# Sale/purchase prices 
c_sell_save_c = np.zeros((H,(n_cases+1)*4,n_cases,years-6)) 
c_buy_save_c = np.zeros((H,(n_cases+1)*4,n_cases,years-6)) 

# The input timeseries are defined for each year of the training dataset (from 2005 to 2014)
# and each number of typical days "c" for each season (from 2 to 7)
for year in range(2005,2015):          # different years of the training dataset
    for c in range(2,n_cases+2):       # different numbers of clusters
        SolRad = np.zeros((H,c*4))       # solar irradiance
        Tamb = np.zeros((H,c*4))         # ambient temperature
        DemEl_res = np.zeros((H,c*4))    # electricity demand of Res
        DemEl_ter = np.zeros((H,c*4))    # electricity demand of Ter
        DemEl_com = np.zeros((H,c*4))    # electricity demand of Com
        DemEl_pub = np.zeros((H,c*4))    # electricity demand of Pub
        DemTh_res = np.zeros((H,c*4))    # heating demand of Res 
        DemTh_ter = np.zeros((H,c*4))    # heating demand of Ter
        DemTh_com = np.zeros((H,c*4))    # heating demand of Com
        DemTh_pub = np.zeros((H,c*4))    # heating demand of Pub
        csell_el = np.zeros((H,c*4))     # PUN sale prices
        w_scenarios = np.zeros(c*4)      # Weights of scenarios in one year
        for yi in range(c*4):            # different clusters
            SolRad[:,yi] = (typical_days[yi,:,c-2,year-2005][0:24])/1000 # [kW/m^2]
            Tamb[:,yi] = typical_days[yi,:,c-2,year-2005][24:48]
            DemEl_res[:,yi] = typical_days[yi,:,c-2,year-2005][72:96]
            DemEl_ter[:,yi] = typical_days[yi,:,c-2,year-2005][96:120]
            DemEl_com[:,yi] = typical_days[yi,:,c-2,year-2005][120:144]
            DemEl_pub[:,yi] = typical_days[yi,:,c-2,year-2005][144:168]
            DemTh_res[:,yi] = typical_days[yi,:,c-2,year-2005][168:192]
            DemTh_ter[:,yi] = typical_days[yi,:,c-2,year-2005][192:216]
            DemTh_com[:,yi] = typical_days[yi,:,c-2,year-2005][216:240]
            DemTh_pub[:,yi] = typical_days[yi,:,c-2,year-2005][240:264]
            csell_el[:,yi] = (typical_days[yi,:,c-2,year-2005][264:288])/1000*0.45   # [€/kWh]
            w_scenarios[yi] = number_elements_each_cluster[c-2,yi,year-2005]
        cbuy_el = csell_el + 0.20 # purchase price, 0.15 €/kWh is the grid tariff
    
        # Solve the stochastic model for different numbers of typical days "c" for each season
        C_TOT_sto_c[c-2,year-2005], design_cost_c[c-2,year-2005], operation_cost_c[c-2,year-2005], \
        P_imp_c[:,0:c*4,c-2,year-2005], P_exp_c[:,0:c*4,c-2,year-2005], E_cond_c[:,0:c*4,c-2,year-2005], \
        capPV_res_c[c-2,year-2005], capEES_res_c[c-2,year-2005], capHP_res_c[c-2,year-2005], capTES_res_c[c-2,year-2005], SocTES_res_c[:,0:c*4,c-2,year-2005], QcTES_res_c[:,0:c*4,c-2,year-2005], QdTES_res_c[:,0:c*4,c-2,year-2005], \
        capBOIL_ter_c[c-2,year-2005], capBOIL_com_c[c-2,year-2005], \
        capPV_pub_c[c-2,year-2005], capEES_pub_c[c-2,year-2005], capHP_pub_c[c-2,year-2005], capTES_pub_c[c-2,year-2005], opt_gap_c[c-2,year-2005], \
        c_sell_save_c[:,0:c*4,c-2,year-2005], c_buy_save_c[:,0:c*4,c-2,year-2005] = REC_MES_typical_day(c*4)

### TIME TO SOLVE THE MODEL ###
print("--- %s seconds ---" % (time.time() - start_time))

#%% SAVE AND EXPORT RESULTS

PATH = "C:/Users/utente/Desktop/Università Gabriele/Assegno_biennale/Models/Model_ML2/Results/"+address

### Without uncertainty ###

# Save results
np.savez(os.path.join(PATH,'Deterministic_solution_'+address), 
         SolRad_av = SolRad_av, Tamb_av = Tamb_av,
         DemEl_res_av = DemEl_res_av, DemEl_ter_av = DemEl_ter_av, DemEl_com_av = DemEl_com_av, DemEl_pub_av = DemEl_pub_av,
         DemTh_res_av = DemTh_res_av, DemTh_ter_av = DemTh_ter_av, DemTh_com_av = DemTh_com_av, DemTh_pub_av = DemTh_pub_av,
         C_TOT_d = C_TOT_d, design_cost_d = design_cost_d, operation_cost_d = operation_cost_d,
         capPV_res_d = capPV_res_d, capEES_res_d = capEES_res_d, capHP_res_d = capHP_res_d, capTES_res_d = capTES_res_d, 
         capBOIL_ter_d = capBOIL_ter_d, capBOIL_com_d = capBOIL_com_d, 
         capPV_pub_d = capPV_pub_d, capEES_pub_d = capEES_pub_d, capHP_pub_d = capHP_pub_d, capTES_pub_d = capTES_pub_d, 
         opt_gap_d = opt_gap_d, c_sell_save_d = c_sell_save_d, c_buy_save_d = c_buy_save_d,
         SocTES_res_d = SocTES_res_d, QcTES_res_d = QcTES_res_d, QdTES_res_d = QdTES_res_d)

# Export optimal capacities
PATH_export_det = "C:/Users/utente/Desktop/Università Gabriele/Assegno_biennale/Models/Model_ML2/Results/"+address+"/Deterministic_solution_"+address+".xlsx"

index_years = [2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014]

df1=pd.DataFrame(C_TOT_d/1000, index=index_years)
df2=pd.DataFrame(design_cost_d/1000, index=index_years)
df3=pd.DataFrame(operation_cost_d/1000, index=index_years)
df4=pd.DataFrame(opt_gap_d, index=index_years)
df1.index.name = 'Years'
df2.index.name = 'Years'
df3.index.name = 'Years'
df4.index.name = 'Years'

writer = pd.ExcelWriter(PATH_export_det, engine="xlsxwriter")   
df1.to_excel(writer, sheet_name='Total_costs', float_format="%.13f")
df2.to_excel(writer, sheet_name='Design_costs', float_format="%.13f")
df3.to_excel(writer, sheet_name='Operation_costs', float_format="%.13f")
df4.to_excel(writer, sheet_name='Optimal_gaps', float_format="%.13f")

writer.save()


### With uncertainty ###

# Save results
np.savez(os.path.join(PATH,'Stochastic_solutions_'+address), 
          typical_days = typical_days, number_elements_each_cluster = number_elements_each_cluster,
          C_TOT_sto_c = C_TOT_sto_c, design_cost_c = design_cost_c, operation_cost_c = operation_cost_c,
          capPV_res_c = capPV_res_c, capEES_res_c = capEES_res_c, capHP_res_c = capHP_res_c, capTES_res_c = capTES_res_c, 
          capBOIL_ter_c = capBOIL_ter_c, capBOIL_com_c = capBOIL_com_c, 
          capPV_pub_c = capPV_pub_c, capEES_pub_c = capEES_pub_c, capHP_pub_c = capHP_pub_c, capTES_pub_c = capTES_pub_c, 
          opt_gap_c = opt_gap_c, c_sell_save_c = c_sell_save_c, c_buy_save_c = c_buy_save_c,
          SocTES_res_c = SocTES_res_c, QcTES_res_c = QcTES_res_c, QdTES_res_c = QdTES_res_c)

# Export optimal capacities
PATH_export_sto = "C:/Users/utente/Desktop/Università Gabriele/Assegno_biennale/Models/Model_ML2/Results/"+address+"/Stochastic_solutions_"+address+".xlsx"

column_names = ['2005','2006','2007','2008','2009','2010','2011','2012','2013','2014']
index_cases = [8, 12, 16, 20, 24, 28] # different numbers of stochastic scenarios: 8 means 2 scenarios for each season

df1=pd.DataFrame(C_TOT_sto_c/1000, columns=column_names, index=index_cases)
df2=pd.DataFrame(design_cost_c/1000, columns=column_names, index=index_cases)
df3=pd.DataFrame(operation_cost_c/1000, columns=column_names, index=index_cases)
df4=pd.DataFrame(opt_gap_c, columns=column_names, index=index_cases)
df1.index.name = 'Number of stochastic scenarios'
df2.index.name = 'Number of stochastic scenarios'
df3.index.name = 'Number of stochastic scenarios'
df4.index.name = 'Number of stochastic scenarios'

writer = pd.ExcelWriter(PATH_export_sto, engine="xlsxwriter")   
df1.to_excel(writer, sheet_name='Total_costs', float_format="%.13f")
df2.to_excel(writer, sheet_name='Design_costs', float_format="%.13f")
df3.to_excel(writer, sheet_name='Operation_costs', float_format="%.13f")
df4.to_excel(writer, sheet_name='Optimal_gaps', float_format="%.13f")

writer.save()

#%% PLOTS

hours=np.arange(1,25)
hours_new=np.arange(1,26)

### With uncertainty
PATH_plots_sto = "C:/Users/utente/Desktop/Università Gabriele/Assegno_biennale/Models/Model_ML2/Results/"+address+"/Plots/Stochastic_scenarios"
# Choose a year and the number of stochastic scenarios for each season of each year 
print('Choose the year')
Y=int(input())  
print('Choose the number of stochastic scenarios for each season of each year: from 2 to 7')
C=int(input())  

# Solar irradiance
fig1 = plt.figure(constrained_layout=True, figsize=(10,8))
ax1 = fig1.add_subplot()
ax1.set_xlabel('Hour [h]', fontsize=24)
ax1.xaxis.set_tick_params(labelsize = 16)
ax1.set_ylabel('Global solar irradiance [kW/$m^2$]', fontsize=24)
ax1.yaxis.set_tick_params(labelsize = 20)
ax1.grid(True, which='both')
for yi in range(C*4): # number of scenarios in one year
    ax1.plot(hours,(typical_days[yi,:,C-2,Y-2005][0:24])/1000, color = "b")
plt.savefig(os.path.join(PATH_plots_sto,'Solar irradiance.png'))

fontsz = 16
labelsz = 14
# Electricity demands
fig2, ax2 = plt.subplots(figsize=(12, 10), nrows=2, ncols=2)
fig2.tight_layout(pad = 3)
ax2[0,0].grid(True, which='both') 
ax2[0,1].grid(True, which='both') 
ax2[1,0].grid(True, which='both') 
ax2[1,1].grid(True, which='both') 
ax2[0,0].xaxis.set_tick_params(labelsize = labelsz)
ax2[0,0].yaxis.set_tick_params(labelsize = labelsz)
ax2[0,1].xaxis.set_tick_params(labelsize = labelsz)
ax2[0,1].yaxis.set_tick_params(labelsize = labelsz)
ax2[1,0].xaxis.set_tick_params(labelsize = labelsz)
ax2[1,0].yaxis.set_tick_params(labelsize = labelsz)
ax2[1,1].xaxis.set_tick_params(labelsize = labelsz)
ax2[1,1].yaxis.set_tick_params(labelsize = labelsz)

for yi in range(C*4):
    res_el=ax2[0,0].plot(hours, typical_days[yi,:,C-2,Y-2005][72:96], "b", linewidth=2)
    ter_el=ax2[0,1].plot(hours, typical_days[yi,:,C-2,Y-2005][96:120], "b", linewidth=2)
    com_el=ax2[1,0].plot(hours, typical_days[yi,:,C-2,Y-2005][120:144], "b", linewidth=2)
    pub_el=ax2[1,1].plot(hours, typical_days[yi,:,C-2,Y-2005][144:168], "b", linewidth=2)

ax2[0,0].set_ylim(0,40)
ax2[0,0].set_title('Res', fontsize=20)
ax2[0,0].set_ylabel('Electricity demand [kWh]', fontsize=fontsz)
    
ax2[0,1].set_ylim(0,40)
ax2[0,1].set_title('Ter', fontsize=20)

ax2[1,0].set_ylim(0,60)
ax2[1,0].set_title('Com', fontsize=20)
ax2[1,0].set_xlabel('Hour [h]', fontsize=fontsz)
ax2[1,0].set_ylabel('Electricity demand [kWh]', fontsize=fontsz)

ax2[1,1].set_ylim(0,400)
ax2[1,1].set_xlabel('Hour [h]', fontsize=fontsz)
ax2[1,1].set_title('Pub', fontsize=20)

plt.savefig(os.path.join(PATH_plots_sto,'Electricity demand.png'))

# Heating demands
fig3, ax3 = plt.subplots(figsize=(12, 10), nrows=2, ncols=2)
fig3.tight_layout(pad = 3)
ax3[0,0].grid(True, which='both') 
ax3[0,1].grid(True, which='both') 
ax3[1,0].grid(True, which='both') 
ax3[1,1].grid(True, which='both') 
ax3[0,0].xaxis.set_tick_params(labelsize = labelsz)
ax3[0,0].yaxis.set_tick_params(labelsize = labelsz)
ax3[0,1].xaxis.set_tick_params(labelsize = labelsz)
ax3[0,1].yaxis.set_tick_params(labelsize = labelsz)
ax3[1,0].xaxis.set_tick_params(labelsize = labelsz)
ax3[1,0].yaxis.set_tick_params(labelsize = labelsz)
ax3[1,1].xaxis.set_tick_params(labelsize = labelsz)
ax3[1,1].yaxis.set_tick_params(labelsize = labelsz)

for yi in range(C*4):
    res_th=ax3[0,0].plot(hours, typical_days[yi,:,C-2,Y-2005][168:192], "b", linewidth=0.5)
    ter_th=ax3[0,1].plot(hours, typical_days[yi,:,C-2,Y-2005][192:216], "b", linewidth=0.5)
    com_th=ax3[1,0].plot(hours, typical_days[yi,:,C-2,Y-2005][216:240], "b", linewidth=0.5)
    pub_th=ax3[1,1].plot(hours, typical_days[yi,:,C-2,Y-2005][240:264], "b", linewidth=0.5)

ax3[0,0].set_ylim(0,60)
ax3[0,0].set_title('Res', fontsize=20)
ax3[0,0].set_ylabel('Heating demand [kWh]', fontsize=fontsz)
    
ax3[0,1].set_ylim(0,60)
ax3[0,1].set_title('Ter', fontsize=20)

ax3[1,0].set_ylim(0,100)
ax3[1,0].set_title('Com', fontsize=20)
ax3[1,0].set_xlabel('Hour [h]', fontsize=fontsz)
ax3[1,0].set_ylabel('Heating demand [kWh]', fontsize=fontsz)

ax3[1,1].set_ylim(0,1000)
ax3[1,1].set_xlabel('Hour [h]', fontsize=fontsz)
ax3[1,1].set_title('Pub', fontsize=20)

plt.savefig(os.path.join(PATH_plots_sto,'Heating demand.png'))

# Ambient temperature
fig4 = plt.figure(constrained_layout=True, figsize=(10,8))
ax4 = fig4.add_subplot()
ax4.set_xlabel('Hour [h]', fontsize=24)
ax4.xaxis.set_tick_params(labelsize = 16)
ax4.set_ylabel('Ambient temperature [°C]', fontsize=24)
ax4.yaxis.set_tick_params(labelsize = 20)
ax4.grid(True, which='both')
for yi in range(C*4):
    ax4.plot(hours, typical_days[yi,:,C-2,Y-2005][24:48], color = "b")
plt.savefig(os.path.join(PATH_plots_sto,'Ambient temperature.png'))

# Grid sale price
fig5 = plt.figure(constrained_layout=True, figsize=(10,8))
ax5 = fig5.add_subplot()
ax5.set_xlabel('Hour [h]', fontsize=24)
ax5.xaxis.set_tick_params(labelsize = 16)
ax5.set_ylabel('Grid sale price [€/kWh]', fontsize=24)
ax5.yaxis.set_tick_params(labelsize = 20)
ax5.grid(True, which='both')
for yi in range(C*4):
    ax5.plot(hours, (typical_days[yi,:,C-2,Y-2005][264:288])/1000*0.5, color = "b")
plt.savefig(os.path.join(PATH_plots_sto,'Grid sale price.png'))

### Without uncertainty: 4 average seasonal days
PATH_plots_det = "C:/Users/utente/Desktop/Università Gabriele/Assegno_biennale/Models/Model_ML2/Results/"+address+"/Plots/Average_seasonal_days"

# Solar irradiance
fig7 = plt.figure(constrained_layout=True, figsize=(10,8))
ax7 = fig7.add_subplot()
ax7.set_xlabel('Hour [h]', fontsize=24)
ax7.xaxis.set_tick_params(labelsize = 16)
ax7.set_ylabel('Global solar irradiance [kW/$m^2$]', fontsize=24)
ax7.yaxis.set_tick_params(labelsize = 20)
ax7.grid(True, which='both')
for yi in range(4):  # 4 average seasonal days
    ax7.plot(hours, SolRad_av[:,yi,Y-2005], color = "b")
plt.savefig(os.path.join(PATH_plots_det,'Solar irradiance.png'))

# Electricity demands
fig8, ax8 = plt.subplots(figsize=(12, 10), nrows=2, ncols=2)
fig8.tight_layout(pad = 3)
ax8[0,0].grid(True, which='both') 
ax8[0,1].grid(True, which='both') 
ax8[1,0].grid(True, which='both') 
ax8[1,1].grid(True, which='both') 
ax8[0,0].xaxis.set_tick_params(labelsize = labelsz)
ax8[0,0].yaxis.set_tick_params(labelsize = labelsz)
ax8[0,1].xaxis.set_tick_params(labelsize = labelsz)
ax8[0,1].yaxis.set_tick_params(labelsize = labelsz)
ax8[1,0].xaxis.set_tick_params(labelsize = labelsz)
ax8[1,0].yaxis.set_tick_params(labelsize = labelsz)
ax8[1,1].xaxis.set_tick_params(labelsize = labelsz)
ax8[1,1].yaxis.set_tick_params(labelsize = labelsz)

for yi in range(4):
    res_el=ax8[0,0].plot(hours, DemEl_res_av[:,yi,Y-2005], "b", linewidth=2)
    ter_el=ax8[0,1].plot(hours, DemEl_ter_av[:,yi,Y-2005], "b", linewidth=2)
    com_el=ax8[1,0].plot(hours, DemEl_com_av[:,yi,Y-2005], "b", linewidth=2)
    pub_el=ax8[1,1].plot(hours, DemEl_pub_av[:,yi,Y-2005], "b", linewidth=2)

ax8[0,0].set_ylim(0,40)
ax8[0,0].set_title('Res', fontsize=20)
ax8[0,0].set_ylabel('Electricity demand [kWh]', fontsize=fontsz)
    
ax8[0,1].set_ylim(0,40)
ax8[0,1].set_title('Ter', fontsize=20)

ax8[1,0].set_ylim(0,60)
ax8[1,0].set_title('Com', fontsize=20)
ax8[1,0].set_xlabel('Hour [h]', fontsize=fontsz)
ax8[1,0].set_ylabel('Electricity demand [kWh]', fontsize=fontsz)

ax8[1,1].set_ylim(0,400)
ax8[1,1].set_xlabel('Hour [h]', fontsize=fontsz)
ax8[1,1].set_title('Pub', fontsize=20)

plt.savefig(os.path.join(PATH_plots_det,'Electricity demand.png'))

# Heating demands
fig9, ax9 = plt.subplots(figsize=(12, 10), nrows=2, ncols=2)
fig9.tight_layout(pad = 3)
ax9[0,0].grid(True, which='both') 
ax9[0,1].grid(True, which='both') 
ax9[1,0].grid(True, which='both') 
ax9[1,1].grid(True, which='both') 
ax9[0,0].xaxis.set_tick_params(labelsize = labelsz)
ax9[0,0].yaxis.set_tick_params(labelsize = labelsz)
ax9[0,1].xaxis.set_tick_params(labelsize = labelsz)
ax9[0,1].yaxis.set_tick_params(labelsize = labelsz)
ax9[1,0].xaxis.set_tick_params(labelsize = labelsz)
ax9[1,0].yaxis.set_tick_params(labelsize = labelsz)
ax9[1,1].xaxis.set_tick_params(labelsize = labelsz)
ax9[1,1].yaxis.set_tick_params(labelsize = labelsz)

for yi in range(4):
    res_th=ax9[0,0].plot(hours, DemTh_res_av[:,yi,Y-2005], "b", linewidth=0.5)
    ter_th=ax9[0,1].plot(hours, DemTh_ter_av[:,yi,Y-2005], "b", linewidth=0.5)
    com_th=ax9[1,0].plot(hours, DemTh_com_av[:,yi,Y-2005], "b", linewidth=0.5)
    pub_th=ax9[1,1].plot(hours, DemTh_pub_av[:,yi,Y-2005], "b", linewidth=0.5)

ax9[0,0].set_ylim(0,60)
ax9[0,0].set_title('Res', fontsize=20)
ax9[0,0].set_ylabel('Heating demand [kWh]', fontsize=fontsz)
    
ax9[0,1].set_ylim(0,60)
ax9[0,1].set_title('Ter', fontsize=20)

ax9[1,0].set_ylim(0,100)
ax9[1,0].set_title('Com', fontsize=20)
ax9[1,0].set_xlabel('Hour [h]', fontsize=fontsz)
ax9[1,0].set_ylabel('Heating demand [kWh]', fontsize=fontsz)

ax9[1,1].set_ylim(0,1000)
ax9[1,1].set_xlabel('Hour [h]', fontsize=fontsz)
ax9[1,1].set_title('Pub', fontsize=20)

plt.savefig(os.path.join(PATH_plots_det,'Heating demand.png'))

# Ambient temperature
fig10 = plt.figure(constrained_layout=True, figsize=(10,8))
ax10 = fig10.add_subplot()
ax10.set_xlabel('Hour [h]', fontsize=24)
ax10.xaxis.set_tick_params(labelsize = 16)
ax10.set_ylabel('Ambient temperature [°C]', fontsize=24)
ax10.yaxis.set_tick_params(labelsize = 20)
ax10.grid(True, which='both')
for yi in range(4):
    ax10.plot(hours, Tamb_av[:,yi,Y-2005], color = "b")
plt.savefig(os.path.join(PATH_plots_det,'Ambient temperature.png'))

# Grid sale price
fig11 = plt.figure(constrained_layout=True, figsize=(10,8))
ax11 = fig11.add_subplot()
ax11.set_xlabel('Hour [h]', fontsize=24)
ax11.xaxis.set_tick_params(labelsize = 16)
ax11.set_ylabel('Grid sale price [€/kWh]', fontsize=24)
ax11.yaxis.set_tick_params(labelsize = 20)
ax11.grid(True, which='both')
for yi in range(4):
    ax11.plot(hours, csell_el_av[:,yi,Y-2005], color = "b")
plt.savefig(os.path.join(PATH_plots_det,'Grid sale price.png'))