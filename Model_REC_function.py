

# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 09:37:45 2024

@author: utente
"""

### REC MODEL

# Optimization of a Renewable Energy Community (REC)
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

### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
### Clustering-based design-operation optimization
### The model is solved in the training dataset (2005-2017)
### With uncertainty in PV capacity factor, ambient temperature and zonal price
#### - Annual (from 2-4 to 50-52 clusters) or seasonal (from 2-4 to 6-8 clusters for each season) clustering 
#### - Replacing or adding approaches for extreme days 

### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Full timeseries operation-only optimization
### The model is solved in the testing dataset (2018-2023)
### - fixed the optimal sizes according to the "best" clustering-based design solutions
### - investment costs refer to the entire time horizon of the testing dataset
###   to align with the operational costs


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# TWO FUNCTIONS
# "clustering_design_operation": clustering-based design-operation optimization
# "fullts_operation": full timeseries operation-only optimization

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# CLUSTERING-BASED DESIGN-OPERATION OPTIMIZATION

# Inputs:
# - loc: name of location for weather variables (solar irradiance, ambient temperature) 
#   and zonal electricity prices (type: string, "Padova", "Palermo", "Trieste")
# - seasons: parameter deciding whether to consider or not division of the
#   available dataset into seasons (type: int, 0: no, 1: yes, default: 0)
# - extreme_method: criterion to take into account extreme days
#   (type: string, "A": adding criterion, "R": replacing criterion, default: 'A')

def clustering_design_operation(loc="",
                                seasons=0,
                                extreme_method='A'):

    ### PYTHON PACKAGES
    
    import gurobipy as grb            # Gurobi, for optimization
    import numpy as np                # Numpy, for multi-dimensional array
    import pandas as pd               # Pandas, for data analysis
    import matplotlib.pyplot as plt   # Pyplot, for plots, graphs, diagrams
    import time                       # time required to solve the optimization
    import os                         # to select the folder for saving
    
    cfold = os.getcwd()  # working directory

    ######################## TECHNO-ECONOMIC PARAMETERS #######################
    
    ifold = os.path.join(cfold, 'Input_data') # folder with input data
    
    # TECHNO-ECONOMIC PARAMETERS OF EACH TECHNOLOGY
    TEparam = pd.read_excel(ifold+"/TechnoEconomicData.xlsx",
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
    ECparam = pd.read_excel(ifold+"/TechnoEconomicData.xlsx",
                            sheet_name='Carriers',
                            index_col=[0,1])
    # Natural gas
    # Purchasing price of natural gas in [€/Mwh]
    cbuy_ng = ECparam.loc["NatGas","CostHouseholds"]["Value"]
    cbuy_ng = cbuy_ng/1000 # from [€/MWh] to [€/kWh]
    
    # FINANCIAL PARAMETERS (INTEREST RATE)
    Fparam = pd.read_excel(ifold+"/TechnoEconomicData.xlsx",
                           sheet_name='Parameters',
                           index_col=[0])
    r = Fparam.loc["InterestRate"]["Value"]
    
    # MAXIMUM HOURLY VARIATION OF FLEXIBLE DEMANDS 
    D_var = 0.1

    ############ CLUSTERING-BASED DESIGN-OPERATION OPTIMIZATION MODEL #########
    
    array=[] # to save the optimal values of the decision variables
        
    def REC_clustering_design_operation(K): # K is the number of days considered in one year
    
        # Initialize the Gurobi model
        m = grb.Model()
        
        ### DECISION VARIABLES for the community
        # Total net electrical power exchanged with the grid
        # Imported power
        P_imp = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="P_imp")
        # Exported power
        P_exp = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS, name="P_exp")
        # Electricity shared (that is incentivised), minimum between "P_imp" and "P_exp"
        # Shared energy according to the italian legislation
        Econd = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS,name='Econd') 
        
        # TIP (Tariffa premio incentivante)
        TIP = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS,name='TIP')
        # Incentive for shared energy
        inc = m.addVars(H, K, vtype=grb.GRB.CONTINUOUS,name='inc')
        
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
        PPV_res = [[CapPV_res*SolRad[h,k] for k in range(K)] for h in range(H)] # SolRad in [kW/kWp]
        
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
        PPV_pub = [[CapPV_pub*SolRad[h,k] for k in range(K)] for h in range(H)] # SolRad in [kW/kWp]
            
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
        # m.addConstrs(Econd[h,k] == min_(P_imp[h,k],P_exp[h,k]) for h in range(H) for k in range(K))
        ### As alternative ###
        # Relaxation of the min constraint (it works if Econd is maximized in the obj function)
        for h in range(H):
            for k in range(K):
                m.addConstr(Econd[h,k] <= P_imp[h,k])
                m.addConstr(Econd[h,k] <= P_exp[h,k])

        ### Definition of the incentive: TIP + restituzione oneri di trasmissione (10.57) [€/MWh]
        # PZ: prezzo zonale [€/MWh]
        # inc: total incentive [€/kWh]
        # Total generation capacity
        CapPV_tot = CapPV_res + CapPV_pub
        
        # Define binary variables for the three cases
        z1 = m.addVar(vtype=grb.GRB.BINARY, name="z1")  # CapPV_tot <= 200
        z2 = m.addVar(vtype=grb.GRB.BINARY, name="z2")  # 200 < CapPV_tot <= 600
        z3 = m.addVar(vtype=grb.GRB.BINARY, name="z3")  # CapPV_tot > 600
    
        # Ensure only one of the binary variables is active
        m.addConstr(z1 + z2 + z3 == 1)
        
        # Define Big-M (large enough number to deactivate constraints when not active)
        M = 10000  
        
        # Activate the correct range
        m.addConstr(CapPV_tot <= 200 + M * (1 - z1))
        m.addConstr(CapPV_tot >= 201 - M * (1 - z2))
        m.addConstr(CapPV_tot <= 600 + M * (1 - z2))
        m.addConstr(CapPV_tot >= 601 - M * (1 - z3))
        
        # Define TIP_base and TIP_max using the binary variables
        TIP_base = m.addVar(name="TIP_base")
        TIP_max = m.addVar(name="TIP_max")
        m.addConstr(TIP_base == 80 * z1 + 70 * z2 + 60 * z3)
        m.addConstr(TIP_max == 120 * z1 + 110 * z2 + 100 * z3)
    
        # Constraints using Gurobi's built-in min/max functions
        for h in range(H):
            for k in range(K):
                TIP_var = m.addVar(name=f"TIP_{h}_{k}")  # Decision variable for TIP
                max_var = m.addVar(name=f"maxVar_{h}_{k}")  # Auxiliary variable for max constraint
        
                # max_var = max(0, 180 - PZ[h,k])
                m.addGenConstrMax(max_var, [0, 180 - PZ[h,k]])
        
                # TIP_var = min(TIP_max, TIP_base + max_var) 
                temp_var = m.addVar(name=f"tempVar_{h}_{k}")  # Auxiliary variable for TIP_base + max_var
                m.addConstr(temp_var == TIP_base + max_var)
                m.addGenConstrMin(TIP_var, [TIP_max, temp_var])
                ### Sum with location parameter 
                if loc=='Padova' or loc=='Trieste':
                   m.addConstr(TIP[h,k] == TIP_var + 10)
                elif loc=='Palermo':   
                     m.addConstr(TIP[h,k] == TIP_var + 0)
                     
                # Constraint for `inc`
                m.addConstr(inc[h,k] == (TIP[h,k] + 10.57) / 1000)
        
        m.update()
        
        ### DEFINE COSTS
        
        ### The solution is based on typical days 
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
                                         + P_imp[h,k]*cbuy_el[h,k] - P_exp[h,k]*csell_el[h,k] - Econd[h,k]*inc[h,k]
                                         for h in range(H)) for k in range(K))
        
        m.update()
            
        ### SET THE OBJECTIVE FUNCTION
        
        m.setObjective(CostINV+CostOP, grb.GRB.MINIMIZE)    
        
        opt_gap = []
        def data_cb(m, where):
            if where == grb.GRB.Callback.MIPSOL:
               cur_obj = m.cbGet(grb.GRB.Callback.MIPSOL_OBJBST)  # best incumbent solution (upper bound, primal solution)
               cur_bd = m.cbGet(grb.GRB.Callback.MIPSOL_OBJBND)   # best bound solution (lower bound, dual solution)
               opt_gap.append(((cur_obj-cur_bd)/cur_obj)*100)     # all gaps during iterations
        
        m.Params.NonConvex = 2           # non-convex MIQCP model
        # Set extra termination criteria
        m.Params.MIPGap = 2e-2           # tolerance of the optimal solution in []
        
        m.optimize(callback=data_cb)     # optimization with callback to retrieve optimal gaps
        
        if len(opt_gap) == 0:  # Handle empty case
           print("Warning: Callback did not run, retrieving gap from solver.")
           if m.status == grb.GRB.OPTIMAL:
              opt_gap = m.Params.MIPGap * 100  # Get gap directly from solver
           else:
              opt_gap = np.nan  # Assign NaN if no solution was found
        else:
              opt_gap = np.array(opt_gap)
              print(f"opt_gap contents: {opt_gap}")
              opt_gap = opt_gap[-1]  # Get final gap
    
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
        array_TIP=[]
        array_inc=[]
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
                array_TIP.append(TIP[h,k].X)
                array_inc.append(inc[h,k].X)
        
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
        array_TIP=np.array(array_TIP)
        array_inc=np.array(array_inc)
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
        array_TIP=np.reshape(array_TIP,(H,K),order='C')
        array_inc=np.reshape(array_inc,(H,K),order='C')
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
        
        return  (C_TOT, design_cost, operation_cost, 
        array_P_imp, array_P_exp, 
        array_Econd, PZ, array_TIP, array_inc, 
        cap_PV_res, cap_EES_res, cap_HP_res, cap_TES_res, 
        cap_BOIL_ter, cap_BOIL_com, 
        cap_PV_pub, cap_EES_pub, cap_HP_pub, cap_TES_pub, 
        opt_gap, csell_el, cbuy_el, 
        array_DemEl_res_op, array_DemEl_ter_op, array_DemEl_com_op, array_DemEl_pub_op)

    ################## AVERAGE DAYS AND STOCHASTIC SCENARIOS ##################
    
    clust_fold = os.path.join(cfold, 'Clustering')   # folder with clustering results
    PATH_clustering = clust_fold+"/"+loc
    
    # Length of a day
    H = 24
    
    if loc=='Padova' or loc=='Trieste':
       c_price = 0.5
    elif loc=='Palermo':
         c_price = 0.45
    
    if seasons==0: ### WITHOUT DIVISION BY SEASONS
       timeseries = np.load(os.path.join(PATH_clustering,'Clustering_train_'+extreme_method+'.npz'), allow_pickle=True)
       years_train = timeseries['years_train']
       
       # # Average annual days with their weights in the training and testing datasets
       # DF_no_s_mean = timeseries['DF_no_s_mean']
       # # Average annual days as (H,years) 
       # SolRad_av = np.zeros((H,years))       # PV capacity factor
       # Tamb_av = np.zeros((H,years))         # ambient temperature
       # csell_el_av = np.zeros((H,years))     # zonal price
       # DemEl_res_av = np.zeros((H,years))    # electricity demand of Res
       # DemEl_ter_av = np.zeros((H,years))    # electricity demand of Ter
       # DemEl_com_av = np.zeros((H,years))    # electricity demand of Com
       # DemEl_pub_av = np.zeros((H,years))    # electricity demand of Pub
       # DemTh_res_av = np.zeros((H,years))    # heating demand of Res 
       # DemTh_ter_av = np.zeros((H,years))    # heating demand of Ter
       # DemTh_com_av = np.zeros((H,years))    # heating demand of Com
       # DemTh_pub_av = np.zeros((H,years))    # heating demand of Pub
       # for year in range(year_start,year_end+1):    # different years of the training dataset
       #     SolRad_av[:,year-year_start] = DF_no_s_mean[year-year_start,:,1]/1000 # [kW/kWp]
       #     Tamb_av[:,year-year_start] = DF_no_s_mean[year-year_start,:,2]
       #     csell_el_av[:,year-year_start] = DF_no_s_mean[year-year_start,:,3]/1000*c_price # [€/kWh]
       #     DemEl_res_av[:,year-year_start] = DF_no_s_mean[year-year_start,:,4]
       #     DemEl_ter_av[:,year-year_start] = DF_no_s_mean[year-year_start,:,5]
       #     DemEl_com_av[:,year-year_start] = DF_no_s_mean[year-year_start,:,6]
       #     DemEl_pub_av[:,year-year_start] = DF_no_s_mean[year-year_start,:,7]
       #     DemTh_res_av[:,year-year_start] = DF_no_s_mean[year-year_start,:,8]
       #     DemTh_ter_av[:,year-year_start] = DF_no_s_mean[year-year_start,:,9]
       #     DemTh_com_av[:,year-year_start] = DF_no_s_mean[year-year_start,:,10]
       #     DemTh_pub_av[:,year-year_start] = DF_no_s_mean[year-year_start,:,11]
       # cbuy_el_av = csell_el_av + 0.20 # purchase price, 0.15 €/kWh is the grid tariff
       
       # Typical days
       if extreme_method == 'R':  # Replace
          # Weights of the typical days in one year (weight in years_train years divided by the number of years)
          number_elements_each_cluster = timeseries['number_elements_each_cluster']/years_train
       elif extreme_method == 'A': # Adding
          # Weights of the typical days in one year (weight in years_train years divided by the number of years)
          number_elements_each_cluster = timeseries['number_elements_each_cluster_new']/years_train
       
       # Number of columns in the training dataset
       lenghth_features_c = timeseries['lenghth_features_c']
       # Indexes of typical days (0-4745)
       index_representative_r_new = timeseries['index_representative_r_new']
       # Training dataset with all timeseries in series
       train_c = timeseries['train_c']
       # Different numbers of clusters (typical days) assessed
       n_cases = timeseries['n_cases']    
       
       # Results for each case (number of typical days "c") and each year of the training dataset
       if extreme_method == 'R':  # Replace
          typ_days = n_cases+1
       elif extreme_method == 'A': # Adding
            typ_days = n_cases+3
       
       # Typical days of PV capacity factor, ambient temperature and zonal prices
       # found over different years
       if extreme_method == 'R':  # Replace
           typical_days = np.zeros((typ_days,lenghth_features_c,n_cases))
           for c in range(2,n_cases+2): # different numbers of clusters
               for yi in range(c):      # different clusters
                   typical_days[yi,:,c-2] = train_c[int(index_representative_r_new[c-2,yi])]
       elif extreme_method == 'A': # Adding
             typical_days = np.zeros((typ_days,lenghth_features_c,n_cases))
             for c in range(2,n_cases+2): # different numbers of clusters
                 for yi in range(c+2):      # different clusters
                     typical_days[yi,:,c-2] = train_c[int(index_representative_r_new[c-2,yi])]
       typical_days = np.array(typical_days)
    
    elif seasons==1: ### WITH DIVISION BY SEASONS
         timeseries_w = np.load(os.path.join(PATH_clustering,'Clustering_train_'+extreme_method+'_'+'winter.npz'), allow_pickle=True)
         timeseries_sp = np.load(os.path.join(PATH_clustering,'Clustering_train_'+extreme_method+'_'+'spring.npz'), allow_pickle=True)
         timeseries_su = np.load(os.path.join(PATH_clustering,'Clustering_train_'+extreme_method+'_'+'summer.npz'), allow_pickle=True)
         timeseries_a = np.load(os.path.join(PATH_clustering,'Clustering_train_'+extreme_method+'_'+'autumn.npz'), allow_pickle=True)
         years_train = timeseries_w['years_train']
         n_cases = timeseries_w['n_cases']# Different numbers of clusters (typical days) assessed
         lenghth_features_c = timeseries_w['lenghth_features_c']
         
         # # Average seasonal days with their weights in the training and testing datasets
         # DF_winter_mean = timeseries_w['DF_winter_mean']
         # DF_spring_mean = timeseries_sp['DF_spring_mean']
         # DF_summer_mean = timeseries_su['DF_summer_mean']
         # DF_autumn_mean = timeseries_a['DF_autumn_mean']
         # N_days_winter = timeseries_w['N_days_winter']
         # N_days_spring = timeseries_sp['N_days_spring']
         # N_days_summer = timeseries_su['N_days_summer']
         # N_days_autumn = timeseries_a['N_days_autumn']
         # # Average seasonal days as (H,4,years) aggregating the different seasons in all years
         # SolRad_av = np.zeros((H*4,years))       # PV capacity factor
         # Tamb_av = np.zeros((H*4,years))         # ambient temperature
         # csell_el_av = np.zeros((H*4,years))     # PUN sale prices
         # DemEl_res_av = np.zeros((H*4,years))    # electricity demand of Res
         # DemEl_ter_av = np.zeros((H*4,years))    # electricity demand of Ter
         # DemEl_com_av = np.zeros((H*4,years))    # electricity demand of Com
         # DemEl_pub_av = np.zeros((H*4,years))    # electricity demand of Pub
         # DemTh_res_av = np.zeros((H*4,years))    # heating demand of Res 
         # DemTh_ter_av = np.zeros((H*4,years))    # heating demand of Ter
         # DemTh_com_av = np.zeros((H*4,years))    # heating demand of Com
         # DemTh_pub_av = np.zeros((H*4,years))    # heating demand of Pub
         # for year in range(year_start,year_end+1):    # different years of the training dataset
         #     SolRad_av[:,year-year_start] = np.concatenate([DF_winter_mean[year-year_start,:,1],DF_spring_mean[year-year_start,:,1],DF_summer_mean[year-year_start,:,1],DF_autumn_mean[year-year_start,:,1]])/1000 # [kW/kWp]
         #     Tamb_av[:,year-year_start] = np.concatenate([DF_winter_mean[year-year_start,:,2],DF_spring_mean[year-year_start,:,2],DF_summer_mean[year-year_start,:,2],DF_autumn_mean[year-year_start,:,2]])
         #     csell_el_av[:,year-year_start] = np.concatenate([DF_winter_mean[year-year_start,:,3],DF_spring_mean[year-year_start,:,3],DF_summer_mean[year-year_start,:,3],DF_autumn_mean[year-year_start,:,3]])/1000*c_price # [€/kWh]
         #     DemEl_res_av[:,year-year_start] = np.concatenate([DF_winter_mean[year-year_start,:,4],DF_spring_mean[year-year_start,:,4],DF_summer_mean[year-year_start,:,4],DF_autumn_mean[year-year_start,:,4]])
         #     DemEl_ter_av[:,year-year_start] = np.concatenate([DF_winter_mean[year-year_start,:,5],DF_spring_mean[year-year_start,:,5],DF_summer_mean[year-year_start,:,5],DF_autumn_mean[year-year_start,:,5]])
         #     DemEl_com_av[:,year-year_start] = np.concatenate([DF_winter_mean[year-year_start,:,6],DF_spring_mean[year-year_start,:,6],DF_summer_mean[year-year_start,:,6],DF_autumn_mean[year-year_start,:,6]])
         #     DemEl_pub_av[:,year-year_start] = np.concatenate([DF_winter_mean[year-year_start,:,7],DF_spring_mean[year-year_start,:,7],DF_summer_mean[year-year_start,:,7],DF_autumn_mean[year-year_start,:,7]])
         #     DemTh_res_av[:,year-year_start] = np.concatenate([DF_winter_mean[year-year_start,:,8],DF_spring_mean[year-year_start,:,8],DF_summer_mean[year-year_start,:,8],DF_autumn_mean[year-year_start,:,8]])
         #     DemTh_ter_av[:,year-year_start] = np.concatenate([DF_winter_mean[year-year_start,:,9],DF_spring_mean[year-year_start,:,9],DF_summer_mean[year-year_start,:,9],DF_autumn_mean[year-year_start,:,9]])
         #     DemTh_com_av[:,year-year_start] = np.concatenate([DF_winter_mean[year-year_start,:,10],DF_spring_mean[year-year_start,:,10],DF_summer_mean[year-year_start,:,10],DF_autumn_mean[year-year_start,:,10]])
         #     DemTh_pub_av[:,year-year_start] = np.concatenate([DF_winter_mean[year-year_start,:,11],DF_spring_mean[year-year_start,:,11],DF_summer_mean[year-year_start,:,11],DF_autumn_mean[year-year_start,:,11]])
         # SolRad_av=np.reshape(SolRad_av,(H,4,years),order='F')
         # Tamb_av=np.reshape(Tamb_av,(H,4,years),order='F')
         # csell_el_av=np.reshape(csell_el_av,(H,4,years),order='F')
         # DemEl_res_av=np.reshape(DemEl_res_av,(H,4,years),order='F')
         # DemEl_ter_av=np.reshape(DemEl_ter_av,(H,4,years),order='F')
         # DemEl_com_av=np.reshape(DemEl_com_av,(H,4,years),order='F')
         # DemEl_pub_av=np.reshape(DemEl_pub_av,(H,4,years),order='F')
         # DemTh_res_av=np.reshape(DemTh_res_av,(H,4,years),order='F')
         # DemTh_ter_av=np.reshape(DemTh_ter_av,(H,4,years),order='F')
         # DemTh_com_av=np.reshape(DemTh_com_av,(H,4,years),order='F')
         # DemTh_pub_av=np.reshape(DemTh_pub_av,(H,4,years),order='F')
         # cbuy_el_av = csell_el_av + 0.20 # purchase price, 0.15 €/kWh is the grid tariff
         
         # Typical days
         if extreme_method == 'R':  # Replace
            # Weights of the typical days in one year (weight in years_train years divided by the number of years)
            number_elements_each_cluster_w = timeseries_w['number_elements_each_cluster']/years_train
            number_elements_each_cluster_sp = timeseries_sp['number_elements_each_cluster']/years_train
            number_elements_each_cluster_su = timeseries_su['number_elements_each_cluster']/years_train
            number_elements_each_cluster_a = timeseries_a['number_elements_each_cluster']/years_train
         elif extreme_method == 'A': # Adding
              # Weights of the typical days in one year (weight in years_train years divided by the number of years)
              number_elements_each_cluster_w = timeseries_w['number_elements_each_cluster_new']/years_train
              number_elements_each_cluster_sp = timeseries_sp['number_elements_each_cluster_new']/years_train
              number_elements_each_cluster_su = timeseries_su['number_elements_each_cluster_new']/years_train
              number_elements_each_cluster_a = timeseries_a['number_elements_each_cluster_new']/years_train
              
         # Seasonal datasets with all timeseries over all years (2D)
         train_winter_c = timeseries_w['train_winter_c']
         train_spring_c = timeseries_sp['train_spring_c']
         train_summer_c = timeseries_su['train_summer_c']
         train_autumn_c = timeseries_a['train_autumn_c']
    
         # Indexes of typical days (considering all years)
         index_representative_r_new_w = timeseries_w['index_representative_r_new']
         index_representative_r_new_sp = timeseries_sp['index_representative_r_new']
         index_representative_r_new_su = timeseries_su['index_representative_r_new']
         index_representative_r_new_a = timeseries_a['index_representative_r_new']
         
         # Results for each case (number of typical days "c") and each year of the training dataset
         if extreme_method == 'R':  # Replace
            typ_days = (n_cases+1)*4
         elif extreme_method == 'A': # Adding
              typ_days = (n_cases+3)*4
        
         # Different sets of stochastic scenarios
         # Fix the scenarios of solar irradiance, ambient temperature and zonal prices found across the different years 
         if extreme_method == 'R':  # Replace
            sets_w = np.zeros((n_cases+1,len(train_winter_c[1]),n_cases))
            sets_sp = np.zeros((n_cases+1,len(train_spring_c[1]),n_cases))
            sets_su = np.zeros((n_cases+1,len(train_summer_c[1]),n_cases))
            sets_a = np.zeros((n_cases+1,len(train_autumn_c[1]),n_cases))
            # Stochastic scenarios in all seasons with their weights
            sets = np.zeros((typ_days,len(train_winter_c[1]),n_cases))
            number_elements_each_cluster = np.zeros((n_cases,typ_days))
            for c in range(2,n_cases+2):   # different numbers of clusters
                for yi in range(c):        # different clusters
                    sets_w[yi,:,c-2] = train_winter_c[int(index_representative_r_new_w[c-2,yi])]   # first row contains 2 scenarios
                    sets_sp[yi,:,c-2] = train_spring_c[int(index_representative_r_new_sp[c-2,yi])]
                    sets_su[yi,:,c-2] = train_summer_c[int(index_representative_r_new_su[c-2,yi])]
                    sets_a[yi,:,c-2] = train_autumn_c[int(index_representative_r_new_a[c-2,yi])]
                sets[0:c*4,:,c-2] = np.concatenate([sets_w[0:c,:,c-2],
                                                    sets_sp[0:c,:,c-2],
                                                    sets_su[0:c,:,c-2],
                                                    sets_a[0:c,:,c-2]])
                number_elements_each_cluster[c-2,0:c*4] = np.concatenate([number_elements_each_cluster_w[c-2,0:c],
                                                                          number_elements_each_cluster_sp[c-2,0:c],
                                                                          number_elements_each_cluster_su[c-2,0:c],
                                                                          number_elements_each_cluster_a[c-2,0:c]])
        
         elif extreme_method == 'A': # Adding
              sets_w = np.zeros((n_cases+3,len(train_winter_c[1]),n_cases))
              sets_sp = np.zeros((n_cases+3,len(train_spring_c[1]),n_cases))
              sets_su = np.zeros((n_cases+3,len(train_summer_c[1]),n_cases))
              sets_a = np.zeros((n_cases+3,len(train_autumn_c[1]),n_cases))
              # Stochastic scenarios in all seasons with their weights
              sets = np.zeros((typ_days,len(train_winter_c[1]),n_cases))
              number_elements_each_cluster = np.zeros((n_cases,typ_days))
              for c in range(2,n_cases+2): # different numbers of clusters
                  for yi in range(c+2):  # different clusters
                      sets_w[yi,:,c-2] = train_winter_c[int(index_representative_r_new_w[c-2,yi])]   # first row contains 2 scenarios
                      sets_sp[yi,:,c-2] = train_spring_c[int(index_representative_r_new_sp[c-2,yi])]
                      sets_su[yi,:,c-2] = train_summer_c[int(index_representative_r_new_su[c-2,yi])]
                      sets_a[yi,:,c-2] = train_autumn_c[int(index_representative_r_new_a[c-2,yi])]
                  sets[0:(c+2)*4,:,c-2] = np.concatenate([sets_w[0:c+2,:,c-2],
                                                          sets_sp[0:c+2,:,c-2],
                                                          sets_su[0:c+2,:,c-2],
                                                          sets_a[0:c+2,:,c-2]])
                  number_elements_each_cluster[c-2,0:(c+2)*4] = np.concatenate([number_elements_each_cluster_w[c-2,0:c+2],
                                                                                number_elements_each_cluster_sp[c-2,0:c+2],
                                                                                number_elements_each_cluster_su[c-2,0:c+2],
                                                                                number_elements_each_cluster_a[c-2,0:c+2]])
        
    ### SOLVING THE CLUSTERING-BASED DESIGN-OPERATION OPTIMIZATION MODEL ######
    ### IN THE TRAINING DATASET
    ### Annual clustering:
    ### - from 2 to 50 stochastic scenarios, replacing criterion
    ### - from 4 to 52 stochastic scenarios, adding criterion
    ### Seasonal clustering:
    ### - from 2 to 6 stochastic scenarios for each season, replacing criterion
    ### - from 4 to 8 stochastic scenarios for each season, adding criterion
    
    start_time = time.time()
    
    # Results for each case (different numbers of typical days)
    # Optimal costs
    C_TOT_sto_c = np.zeros(n_cases)
    design_cost_c = np.zeros(n_cases)
    operation_cost_c = np.zeros(n_cases)
    # Optimal operational variables
    P_imp_c = np.zeros((H,typ_days,n_cases))
    P_exp_c = np.zeros((H,typ_days,n_cases))
    E_cond_c = np.zeros((H,typ_days,n_cases))
    SocTES_res_c = np.zeros((H+1,typ_days,n_cases))
    QcTES_res_c = np.zeros((H,typ_days,n_cases))
    QdTES_res_c = np.zeros((H,typ_days,n_cases))
    PZ_c = np.zeros((H,typ_days,n_cases))
    TIP_c = np.zeros((H,typ_days,n_cases))
    inc_c = np.zeros((H,typ_days,n_cases))
    # Optimal capacities
    capPV_res_c = np.zeros(n_cases)
    capEES_res_c = np.zeros(n_cases)
    capHP_res_c = np.zeros(n_cases)
    capTES_res_c = np.zeros(n_cases)
    capBOIL_ter_c = np.zeros(n_cases)
    capBOIL_com_c = np.zeros(n_cases)
    capPV_pub_c = np.zeros(n_cases)
    capEES_pub_c = np.zeros(n_cases)
    capHP_pub_c = np.zeros(n_cases)
    capTES_pub_c = np.zeros(n_cases)
    # Optimal gaps
    opt_gap_c = np.zeros(n_cases)
    # Sale/purchase prices 
    c_sell_save_c = np.zeros((H,typ_days,n_cases)) 
    c_buy_save_c = np.zeros((H,typ_days,n_cases)) 
    # Optimal user demands
    DemEl_res_op_c = np.zeros((H,typ_days,n_cases))
    DemEl_ter_op_c = np.zeros((H,typ_days,n_cases))
    DemEl_com_op_c = np.zeros((H,typ_days,n_cases))
    DemEl_pub_op_c = np.zeros((H,typ_days,n_cases))
    
    # # Coefficients to correct electricity demands
    # corr_res = 1.5
    # corr_ter = 1.5 
    # corr_com = 1.5 
    # corr_pub = 0.4 
    
    # The input timeseries are defined for each number of typical days "c" for each season/year
    
    if seasons==0: ### WITHOUT DIVISION BY SEASONS
       
        for c in range(2,n_cases+2):       # different numbers of clusters
            if extreme_method == 'R':      # Replace
                SolRad = np.zeros((H,c))       # PV capacity factor
                Tamb = np.zeros((H,c))         # ambient temperature
                csell_el = np.zeros((H,c))     # zonal prices
                DemEl_res = np.zeros((H,c))    # electricity demand of Res
                DemEl_ter = np.zeros((H,c))    # electricity demand of Ter
                DemEl_com = np.zeros((H,c))    # electricity demand of Com
                DemEl_pub = np.zeros((H,c))    # electricity demand of Pub
                DemTh_res = np.zeros((H,c))    # heating demand of Res 
                DemTh_ter = np.zeros((H,c))    # heating demand of Ter
                DemTh_com = np.zeros((H,c))    # heating demand of Com
                DemTh_pub = np.zeros((H,c))    # heating demand of Pub
                w_scenarios = np.zeros(c)      # Weights of scenarios in one year
                PZ = np.zeros((H,c))           # zonal market prices
                for yi in range(c):            # different clusters
                    SolRad[:,yi] = (typical_days[yi,:,c-2][0:24])/1000 # [kW/kWp]
                    Tamb[:,yi] = typical_days[yi,:,c-2][24:48]
                    csell_el[:,yi] = (typical_days[yi,:,c-2][48:72])/1000*c_price
                    DemEl_res[:,yi] = typical_days[yi,:,c-2][72:96]
                    DemEl_ter[:,yi] = typical_days[yi,:,c-2][96:120]
                    DemEl_com[:,yi] = typical_days[yi,:,c-2][120:144]
                    DemEl_pub[:,yi] = typical_days[yi,:,c-2][144:168]
                    DemTh_res[:,yi] = typical_days[yi,:,c-2][168:192]
                    DemTh_ter[:,yi] = typical_days[yi,:,c-2][192:216]
                    DemTh_com[:,yi] = typical_days[yi,:,c-2][216:240]
                    DemTh_pub[:,yi] = typical_days[yi,:,c-2][240:264]
                    w_scenarios[yi] = number_elements_each_cluster[c-2,yi]
                    PZ[:,yi] = (typical_days[yi,:,c-2][48:72])  # €/MWh
                cbuy_el = csell_el + 0.20 # purchase price, 0.20 €/kWh is the grid tariff
                
                # Solve the stochastic model for different numbers of typical days "c" for each year
                
                C_TOT_sto_c[c-2], design_cost_c[c-2], operation_cost_c[c-2], \
                P_imp_c[:,0:c,c-2], P_exp_c[:,0:c,c-2], \
                E_cond_c[:,0:c,c-2], PZ_c[:,0:c,c-2], TIP_c[:,0:c,c-2], inc_c[:,0:c,c-2], \
                capPV_res_c[c-2], capEES_res_c[c-2], capHP_res_c[c-2], capTES_res_c[c-2], \
                capBOIL_ter_c[c-2], capBOIL_com_c[c-2], \
                capPV_pub_c[c-2], capEES_pub_c[c-2], capHP_pub_c[c-2], capTES_pub_c[c-2], \
                opt_gap_c[c-2], c_sell_save_c[:,0:c,c-2], c_buy_save_c[:,0:c,c-2], \
                DemEl_res_op_c[:,0:c,c-2], DemEl_ter_op_c[:,0:c,c-2], DemEl_com_op_c[:,0:c,c-2], DemEl_pub_op_c[:,0:c,c-2] = REC_clustering_design_operation(c)
    
            elif extreme_method == 'A': # Adding
                  SolRad = np.zeros((H,c+2))       # PV capacity factor
                  Tamb = np.zeros((H,c+2))         # ambient temperature
                  csell_el = np.zeros((H,c+2))     # zonal prices
                  DemEl_res = np.zeros((H,c+2))    # electricity demand of Res
                  DemEl_ter = np.zeros((H,c+2))    # electricity demand of Ter
                  DemEl_com = np.zeros((H,c+2))    # electricity demand of Com
                  DemEl_pub = np.zeros((H,c+2))    # electricity demand of Pub
                  DemTh_res = np.zeros((H,c+2))    # heating demand of Res 
                  DemTh_ter = np.zeros((H,c+2))    # heating demand of Ter
                  DemTh_com = np.zeros((H,c+2))    # heating demand of Com
                  DemTh_pub = np.zeros((H,c+2))    # heating demand of Pub
                  w_scenarios = np.zeros(c+2)      # Weights of scenarios in one year
                  PZ = np.zeros((H,c+2))           # zonal market prices
                  for yi in range(c+2):            # different clusters
                      SolRad[:,yi] = (typical_days[yi,:,c-2][0:24])/1000 # [kW/kWp]
                      Tamb[:,yi] = typical_days[yi,:,c-2][24:48]
                      csell_el[:,yi] = (typical_days[yi,:,c-2][48:72])/1000*c_price
                      DemEl_res[:,yi] = typical_days[yi,:,c-2][72:96]
                      DemEl_ter[:,yi] = typical_days[yi,:,c-2][96:120]
                      DemEl_com[:,yi] = typical_days[yi,:,c-2][120:144]
                      DemEl_pub[:,yi] = typical_days[yi,:,c-2][144:168]
                      DemTh_res[:,yi] = typical_days[yi,:,c-2][168:192]
                      DemTh_ter[:,yi] = typical_days[yi,:,c-2][192:216]
                      DemTh_com[:,yi] = typical_days[yi,:,c-2][216:240]
                      DemTh_pub[:,yi] = typical_days[yi,:,c-2][240:264]
                      w_scenarios[yi] = number_elements_each_cluster[c-2,yi]
                      PZ[:,yi] = (typical_days[yi,:,c-2][48:72])  # €/MWh
                  cbuy_el = csell_el + 0.20 # purchase price, 0.20 €/kWh is the grid tariff
                  
                  # Solve the stochastic model for different numbers of typical days "c" for each year
                  
                  C_TOT_sto_c[c-2], design_cost_c[c-2], operation_cost_c[c-2], \
                  P_imp_c[:,0:c+2,c-2], P_exp_c[:,0:c+2,c-2], \
                  E_cond_c[:,0:c+2,c-2], PZ_c[:,0:c+2,c-2], TIP_c[:,0:c+2,c-2], inc_c[:,0:c+2,c-2], \
                  capPV_res_c[c-2], capEES_res_c[c-2], capHP_res_c[c-2], capTES_res_c[c-2], \
                  capBOIL_ter_c[c-2], capBOIL_com_c[c-2], \
                  capPV_pub_c[c-2], capEES_pub_c[c-2], capHP_pub_c[c-2], capTES_pub_c[c-2], \
                  opt_gap_c[c-2], c_sell_save_c[:,0:c+2,c-2], c_buy_save_c[:,0:c+2,c-2], \
                  DemEl_res_op_c[:,0:c+2,c-2], DemEl_ter_op_c[:,0:c+2,c-2], DemEl_com_op_c[:,0:c+2,c-2], DemEl_pub_op_c[:,0:c+2,c-2] = REC_clustering_design_operation(c+2)
    
    elif seasons==1: ### WITH DIVISION BY SEASONS
    
        for c in range(2,n_cases+2):       # different numbers of clusters
            
            if extreme_method == 'R':  # Replace
               SolRad = np.zeros((H,c*4))       # PV capacity factor
               Tamb = np.zeros((H,c*4))         # ambient temperature
               csell_el = np.zeros((H,c*4))     # zonal prices
               DemEl_res = np.zeros((H,c*4))    # electricity demand of Res
               DemEl_ter = np.zeros((H,c*4))    # electricity demand of Ter
               DemEl_com = np.zeros((H,c*4))    # electricity demand of Com
               DemEl_pub = np.zeros((H,c*4))    # electricity demand of Pub
               DemTh_res = np.zeros((H,c*4))    # heating demand of Res 
               DemTh_ter = np.zeros((H,c*4))    # heating demand of Ter
               DemTh_com = np.zeros((H,c*4))    # heating demand of Com
               DemTh_pub = np.zeros((H,c*4))    # heating demand of Pub
               w_scenarios = np.zeros(c*4)      # Weights of scenarios in one year
               PZ = np.zeros((H,c*4))           # zonal market prices
               for yi in range(c*4):            # different clusters
                   SolRad[:,yi] = (typical_days[yi,:,c-2][0:24])/1000 # [kW/kWp]
                   Tamb[:,yi] = typical_days[yi,:,c-2][24:48]
                   csell_el[:,yi] = (typical_days[yi,:,c-2][48:72])/1000*c_price   # [€/kWh]
                   DemEl_res[:,yi] = typical_days[yi,:,c-2][72:96]
                   DemEl_ter[:,yi] = typical_days[yi,:,c-2][96:120]
                   DemEl_com[:,yi] = typical_days[yi,:,c-2][120:144]
                   DemEl_pub[:,yi] = typical_days[yi,:,c-2][144:168]
                   DemTh_res[:,yi] = typical_days[yi,:,c-2][168:192]
                   DemTh_ter[:,yi] = typical_days[yi,:,c-2][192:216]
                   DemTh_com[:,yi] = typical_days[yi,:,c-2][216:240]
                   DemTh_pub[:,yi] = typical_days[yi,:,c-2][240:264]
                   w_scenarios[yi] = number_elements_each_cluster[c-2,yi]
                   PZ[:,yi] = (typical_days[yi,:,c-2][48:72])  # €/MWh
               cbuy_el = csell_el + 0.20 # purchase price, 0.15 €/kWh is the grid tariff
               
               # Solve the stochastic model for different numbers of typical days "c" for each season
    
               C_TOT_sto_c[c-2], design_cost_c[c-2], operation_cost_c[c-2], \
               P_imp_c[:,0:c*4,c-2], P_exp_c[:,0:c*4,c-2], \
               E_cond_c[:,0:c*4,c-2], PZ_c[:,0:c*4,c-2], TIP_c[:,0:c*4,c-2], inc_c[:,0:c*4,c-2], \
               capPV_res_c[c-2], capEES_res_c[c-2], capHP_res_c[c-2], capTES_res_c[c-2], \
               capBOIL_ter_c[c-2], capBOIL_com_c[c-2], \
               capPV_pub_c[c-2], capEES_pub_c[c-2], capHP_pub_c[c-2], capTES_pub_c[c-2],  
               opt_gap_c[c-2], c_sell_save_c[:,0:c*4,c-2], c_buy_save_c[:,0:c*4,c-2], \
               DemEl_res_op_c[:,0:c*4,c-2], DemEl_ter_op_c[:,0:c*4,c-2], DemEl_com_op_c[:,0:c*4,c-2], DemEl_pub_op_c[:,0:c*4,c-2] = REC_clustering_design_operation(c*4)
    
            elif extreme_method == 'A': # Additional
                 SolRad = np.zeros((H,(c+2)*4))       # PV capacity factor
                 Tamb = np.zeros((H,(c+2)*4))         # ambient temperature
                 csell_el = np.zeros((H,(c+2)*4))     # zonal prices
                 DemEl_res = np.zeros((H,(c+2)*4))    # electricity demand of Res
                 DemEl_ter = np.zeros((H,(c+2)*4))    # electricity demand of Ter
                 DemEl_com = np.zeros((H,(c+2)*4))    # electricity demand of Com
                 DemEl_pub = np.zeros((H,(c+2)*4))    # electricity demand of Pub
                 DemTh_res = np.zeros((H,(c+2)*4))    # heating demand of Res 
                 DemTh_ter = np.zeros((H,(c+2)*4))    # heating demand of Ter
                 DemTh_com = np.zeros((H,(c+2)*4))    # heating demand of Com
                 DemTh_pub = np.zeros((H,(c+2)*4))    # heating demand of Pub
                 w_scenarios = np.zeros((c+2)*4)      # Weights of scenarios in one year
                 PZ = np.zeros((H,(c+2)*4))           # zonal market prices
                 for yi in range((c+2)*4):            # different clusters
                     SolRad[:,yi] = (typical_days[yi,:,c-2][0:24])/1000 # [kW/kWp]
                     Tamb[:,yi] = typical_days[yi,:,c-2][24:48]
                     csell_el[:,yi] = (typical_days[yi,:,c-2][48:72])/1000*c_price   # [€/kWh]
                     DemEl_res[:,yi] = typical_days[yi,:,c-2][72:96]
                     DemEl_ter[:,yi] = typical_days[yi,:,c-2][96:120]
                     DemEl_com[:,yi] = typical_days[yi,:,c-2][120:144]
                     DemEl_pub[:,yi] = typical_days[yi,:,c-2][144:168]
                     DemTh_res[:,yi] = typical_days[yi,:,c-2][168:192]
                     DemTh_ter[:,yi] = typical_days[yi,:,c-2][192:216]
                     DemTh_com[:,yi] = typical_days[yi,:,c-2][216:240]
                     DemTh_pub[:,yi] = typical_days[yi,:,c-2][240:264]
                     w_scenarios[yi] = number_elements_each_cluster[c-2,yi]
                     PZ[:,yi] = (typical_days[yi,:,c-2][48:72])  # €/MWh
                 cbuy_el = csell_el + 0.20 # purchase price, 0.15 €/kWh is the grid tariff
                 
                 # Solve the stochastic model for different numbers of typical days "c" for each season
    
                 C_TOT_sto_c[c-2], design_cost_c[c-2], operation_cost_c[c-2], \
                 P_imp_c[:,0:(c+2)*4,c-2], P_exp_c[:,0:(c+2)*4,c-2], \
                 E_cond_c[:,0:(c+2)*4,c-2], PZ_c[:,0:(c+2)*4,c-2], TIP_c[:,0:(c+2)*4,c-2], inc_c[:,0:(c+2)*4,c-2], \
                 capPV_res_c[c-2], capEES_res_c[c-2], capHP_res_c[c-2], capTES_res_c[c-2], \
                 capBOIL_ter_c[c-2], capBOIL_com_c[c-2], \
                 capPV_pub_c[c-2], capEES_pub_c[c-2], capHP_pub_c[c-2], capTES_pub_c[c-2], \
                 opt_gap_c[c-2], c_sell_save_c[:,0:(c+2)*4,c-2], c_buy_save_c[:,0:(c+2)*4,c-2], \
                 DemEl_res_op_c[:,0:(c+2)*4,c-2], DemEl_ter_op_c[:,0:(c+2)*4,c-2], DemEl_com_op_c[:,0:(c+2)*4,c-2], DemEl_pub_op_c[:,0:(c+2)*4,c-2] = REC_clustering_design_operation((c+2)*4)
    
    ### TIME TO SOLVE THE MODEL ###
    print("--- %s seconds ---" % (time.time() - start_time))

    ################### SAVE AND EXPORT RESULTS OF TRAINING ###################
    
    cfold = os.getcwd()  # working directory
    results_fold = os.path.join(cfold, 'Results')
    PATH_results = results_fold+"/"+loc
    
    if seasons==0: ### WITHOUT DIVISION BY SEASONS
    
        # Save results
        np.savez(os.path.join(PATH_results,'Stochastic_solutions_annual_'+extreme_method+'_'+loc), 
                  typical_days = typical_days, number_elements_each_cluster = number_elements_each_cluster,
                  C_TOT_sto_c = C_TOT_sto_c, design_cost_c = design_cost_c, operation_cost_c = operation_cost_c,
                  E_cond_c = E_cond_c, PZ_c = PZ_c, TIP_c = TIP_c, inc_c = inc_c,
                  capPV_res_c = capPV_res_c, capEES_res_c = capEES_res_c, capHP_res_c = capHP_res_c, capTES_res_c = capTES_res_c, 
                  capBOIL_ter_c = capBOIL_ter_c, capBOIL_com_c = capBOIL_com_c, 
                  capPV_pub_c = capPV_pub_c, capEES_pub_c = capEES_pub_c, capHP_pub_c = capHP_pub_c, capTES_pub_c = capTES_pub_c, 
                  opt_gap_c = opt_gap_c, c_sell_save_c = c_sell_save_c, c_buy_save_c = c_buy_save_c,
                  DemEl_res_op_c = DemEl_res_op_c, DemEl_ter_op_c = DemEl_ter_op_c,
                  DemEl_com_op_c = DemEl_com_op_c, DemEl_pub_op_c = DemEl_pub_op_c)
        
        # Export optimal capacities
        PATH_export_sto = PATH_results+"/Stochastic_solutions_annual_"+extreme_method+'_'+loc+".xlsx"
        
        # Different numbers of stochastic scenarios
        if extreme_method == 'R':  # Replace
           index_cases = [a for a in range(2,n_cases+2)] 
        elif extreme_method == 'A': # Additional
             index_cases = [a+2 for a in range(2,n_cases+2)] 
             
        df1=pd.DataFrame(C_TOT_sto_c/1000, index=index_cases)
        df2=pd.DataFrame(design_cost_c/1000, index=index_cases)
        df3=pd.DataFrame(operation_cost_c/1000, index=index_cases)
        df4=pd.DataFrame(data=np.array([capPV_res_c, capEES_res_c, capHP_res_c, capTES_res_c,
                                        capBOIL_ter_c, capBOIL_com_c,
                                        capPV_pub_c, capEES_pub_c, capHP_pub_c, capTES_pub_c]).T,
                                        index=index_cases,
                                        columns=['capPV_res [kWp]', 'capEES_res [kWh]', 'capHP_res [kW]', 'capTES_res [kWh]',
                                                 'capBOIL_ter [kW]', 'capBOIL_com [kW]',
                                                 'capPV_pub [kWp]', 'capEES_pub [kWh]', 'capHP_pub [kW]', 'capTES_pub [kWh]'])
        
        df1.index.name = 'Number of stochastic scenarios'
        df2.index.name = 'Number of stochastic scenarios'
        df3.index.name = 'Number of stochastic scenarios'
        df4.index.name = 'Number of stochastic scenarios'
        
        writer = pd.ExcelWriter(PATH_export_sto, engine="xlsxwriter")   
        df1.to_excel(writer, sheet_name='Total_costs', float_format="%.13f")
        df2.to_excel(writer, sheet_name='Design_costs', float_format="%.13f")
        df3.to_excel(writer, sheet_name='Operation_costs', float_format="%.13f")
        df4.to_excel(writer, sheet_name='Optimal sizes', float_format="%.13f")
        
        writer.close()
    
    elif seasons==1: ### WITH DIVISION BY SEASONS
    
        # Save results
        np.savez(os.path.join(PATH_results,'Stochastic_solutions_seasons_'+extreme_method+'_'+loc), 
                  typical_days = typical_days, number_elements_each_cluster = number_elements_each_cluster,
                  C_TOT_sto_c = C_TOT_sto_c, design_cost_c = design_cost_c, operation_cost_c = operation_cost_c,
                  E_cond_c = E_cond_c, PZ_c = PZ_c, TIP_c = TIP_c, inc_c = inc_c,
                  capPV_res_c = capPV_res_c, capEES_res_c = capEES_res_c, capHP_res_c = capHP_res_c, capTES_res_c = capTES_res_c, 
                  capBOIL_ter_c = capBOIL_ter_c, capBOIL_com_c = capBOIL_com_c, 
                  capPV_pub_c = capPV_pub_c, capEES_pub_c = capEES_pub_c, capHP_pub_c = capHP_pub_c, capTES_pub_c = capTES_pub_c, 
                  opt_gap_c = opt_gap_c, c_sell_save_c = c_sell_save_c, c_buy_save_c = c_buy_save_c,
                  DemEl_res_op_c = DemEl_res_op_c, DemEl_ter_op_c = DemEl_ter_op_c,
                  DemEl_com_op_c = DemEl_com_op_c, DemEl_pub_op_c = DemEl_pub_op_c)
        
        # Export optimal capacities
        PATH_export_sto = PATH_results+"/Stochastic_solutions_seasons_"+extreme_method+'_'+loc+".xlsx"
        
        # Different numbers of stochastic scenarios: 
        # R: 8 means 2 scenarios for each season
        # A: 16 means 4 scenarios for each season
        if extreme_method == 'R':  # Replace
           index_cases = [a*4+4 for a in range(1,n_cases+1)] 
        elif extreme_method == 'A': # Additional
             index_cases = [a*4 for a in range(4,n_cases+4)] 
             
        df1=pd.DataFrame(C_TOT_sto_c/1000, index=index_cases)
        df2=pd.DataFrame(design_cost_c/1000, index=index_cases)
        df3=pd.DataFrame(operation_cost_c/1000, index=index_cases)
        df4=pd.DataFrame(data=np.array([capPV_res_c, capEES_res_c, capHP_res_c, capTES_res_c,
                                        capBOIL_ter_c, capBOIL_com_c,
                                        capPV_pub_c, capEES_pub_c, capHP_pub_c, capTES_pub_c]).T,
                                        index=index_cases,
                                        columns=['capPV_res [kWp]', 'capEES_res [kWh]', 'capHP_res [kW]', 'capTES_res [kWh]',
                                                 'capBOIL_ter [kW]', 'capBOIL_com [kW]',
                                                 'capPV_pub [kWp]', 'capEES_pub [kWh]', 'capHP_pub [kW]', 'capTES_pub [kWh]'])
    
        df1.index.name = 'Number of stochastic scenarios'
        df2.index.name = 'Number of stochastic scenarios'
        df3.index.name = 'Number of stochastic scenarios'
        df4.index.name = 'Number of stochastic scenarios'
        
        writer = pd.ExcelWriter(PATH_export_sto, engine="xlsxwriter")   
        df1.to_excel(writer, sheet_name='Total_costs', float_format="%.13f")
        df2.to_excel(writer, sheet_name='Design_costs', float_format="%.13f")
        df3.to_excel(writer, sheet_name='Operation_costs', float_format="%.13f")
        df4.to_excel(writer, sheet_name='Optimal sizes', float_format="%.13f")
        
        writer.close()
    
    ################### PLOTS OF OPTIMAL SIZES AND TOTAL COSTS ################
    
    cfold = os.getcwd()  # working directory
    
    ### Import the optimal sizes and total costs of the clustering-based design solutions
    results_fold = os.path.join(cfold, 'Results')
    PATH_results = results_fold+"/"+loc
    if seasons==0:
       results_training = np.load(os.path.join(PATH_results,'Stochastic_solutions_annual_'+extreme_method+'_'+loc+'.npz'), allow_pickle=True)
    else:
       results_training = np.load(os.path.join(PATH_results,'Stochastic_solutions_seasons_'+extreme_method+'_'+loc+'.npz'), allow_pickle=True)
    
    C_TOT_sto_c = results_training['C_TOT_sto_c']/1000 # k€
    capPV_res_c = results_training['capPV_res_c']
    capEES_res_c = results_training['capEES_res_c']
    capHP_res_c = results_training['capHP_res_c']
    capTES_res_c = results_training['capTES_res_c']
    capPV_pub_c = results_training['capPV_pub_c']
    capEES_pub_c = results_training['capEES_pub_c']
    capHP_pub_c = results_training['capHP_pub_c']
    capTES_pub_c = results_training['capTES_pub_c']
    n_cases = len(capPV_res_c)  # number of different sets of representative days

    fontsz = 24
    labelsz = 20
    
    # Total costs
    plt.figure(figsize=(10, 6))
    plt.plot(np.arange(4,n_cases+4), C_TOT_sto_c, "b", marker='o', linestyle='--', linewidth=2, alpha=0.7)
    plt.xlabel("Number of representative days", fontsize=fontsz)
    plt.ylabel("Optimal total costs [k€]", fontsize=fontsz)
    plt.title(f"Location: {loc}", fontsize=fontsz, fontweight='bold')
    plt.xticks(fontsize=labelsz)
    plt.yticks(fontsize=labelsz)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(PATH_results, "Clustering_based_costs.jpg"), bbox_inches='tight')
    
    # Optimal sizes
    fig1, ax1 = plt.subplots(figsize=(20, 10), nrows=2, ncols=2)
    fig1.tight_layout(pad = 3)
    fig1.subplots_adjust(left=0.06, wspace=0.13, hspace=0.25, bottom=0.08, top=0.90) 
    plt.suptitle(f"Location: {loc}", fontsize=fontsz, fontweight='bold') 
    for ax in ax1.flatten():
        ax.grid(True)
        ax.xaxis.set_tick_params(labelsize=labelsz)
        ax.yaxis.set_tick_params(labelsize=labelsz)

    ax1[0,0].plot(np.arange(4,n_cases+4), capPV_res_c, "b", marker='o', linestyle='--', linewidth=2, alpha=0.7)
    ax1[0,1].plot(np.arange(4,n_cases+4), capHP_res_c, "b", marker='o', linestyle='--', linewidth=2, alpha=0.7)
    ax1[1,0].plot(np.arange(4,n_cases+4), capPV_pub_c, "b", marker='o', linestyle='--', linewidth=2, alpha=0.7)
    ax1[1,1].plot(np.arange(4,n_cases+4), capHP_pub_c, "b", marker='o', linestyle='--', linewidth=2, alpha=0.7)
    
    titles = [
    "Residential PV", "Residential HP",
    "Public PV", "Public HP"]
    axes = [ax1[0,0], ax1[0,1], ax1[1,0], ax1[1,1]]

    # Loop through subplots to set titles and labels
    for i, ax in enumerate(axes):
        ax.set_title(titles[i], fontsize=fontsz)
        if i >= 2:  
            ax.set_xlabel("Number of representative days", fontsize=fontsz)
    
    ax1[0,0].set_ylabel('Size [kW]', fontsize=fontsz)
    ax1[1,0].set_ylabel('Size [kW]', fontsize=fontsz)
    
    fig1.savefig(os.path.join(PATH_results, "Clustering_based_sizes_PV-HP.jpg"), bbox_inches='tight')
    
    # Optimal sizes of TES
    fig2, ax2 = plt.subplots(figsize=(20, 10), nrows=1, ncols=2)
    fig2.tight_layout(pad = 3)
    fig2.subplots_adjust(left=0.06, wspace=0.13, hspace=0.25, bottom=0.08, top=0.90) 
    plt.suptitle(f"Location: {loc}", fontsize=fontsz, fontweight='bold') 
    for ax in ax2.flatten():
        ax.grid(True)
        ax.xaxis.set_tick_params(labelsize=labelsz)
        ax.yaxis.set_tick_params(labelsize=labelsz)

    ax2[0].plot(np.arange(4,n_cases+4), capTES_res_c, "b", marker='o', linestyle='--', linewidth=2, alpha=0.7)
    ax2[1].plot(np.arange(4,n_cases+4), capTES_pub_c, "b", marker='o', linestyle='--', linewidth=2, alpha=0.7)

    titles = [
    "Residential TES", "Public TES"]
    axes = [ax2[0], ax2[1]]

    # Loop through subplots to set titles and labels
    for i, ax in enumerate(axes):
        ax.set_title(titles[i], fontsize=fontsz) 
        ax.set_xlabel("Number of representative days", fontsize=fontsz)

    ax2[0].set_ylabel('Size [kWh]', fontsize=fontsz)

    fig2.savefig(os.path.join(PATH_results, "Clustering_based_sizes_TES.jpg"), bbox_inches='tight')

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### FULL TIMESERIES OPERATION ONLY OPTIMIZATION
### TESTING THE "BEST" CLUSTERING-BASED DESIGN SOLUTIONS

# Inputs:
# - loc: name of location for weather variables (solar irradiance, ambient temperature) 
#   and zonal electricity prices (type: string, "Padova", "Palermo", "Trieste")
# - seasons: parameter deciding whether to consider or not division of the
#   available dataset into seasons (type: int, 0: no, 1: yes, default: 0)
# - extreme_method: criterion to take into account extreme days
#   (type: string, "A": adding criterion, "R": replacing criterion, default: 'A')
# - "opt_K": number of representative days associated with the "best" 
#    clustering-based design solution to be tested (type: int)

def fullts_operation(loc="",
                     seasons=0,
                     extreme_method='A',
                     opt_K=int()):

    ### PYTHON PACKAGES
    
    import gurobipy as grb            # Gurobi, for optimization
    import numpy as np                # Numpy, for multi-dimensional array
    import pandas as pd               # Pandas, for data analysis
    import matplotlib.pyplot as plt   # Pyplot, for plots, graphs, diagrams
    import time                       # time required to solve the optimization
    import os                         # to select the folder for saving
    
    cfold = os.getcwd()  # working directory

    ######################## TECHNO-ECONOMIC PARAMETERS #######################
    
    ifold = os.path.join(cfold, 'Input_data') # folder with input data
    
    # TECHNO-ECONOMIC PARAMETERS OF EACH TECHNOLOGY
    TEparam = pd.read_excel(ifold+"/TechnoEconomicData.xlsx",
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
    ECparam = pd.read_excel(ifold+"/TechnoEconomicData.xlsx",
                            sheet_name='Carriers',
                            index_col=[0,1])
    # Natural gas
    # Purchasing price of natural gas in [€/Mwh]
    cbuy_ng = ECparam.loc["NatGas","CostHouseholds"]["Value"]
    cbuy_ng = cbuy_ng/1000 # from [€/MWh] to [€/kWh]
    
    # FINANCIAL PARAMETERS (INTEREST RATE)
    Fparam = pd.read_excel(ifold+"/TechnoEconomicData.xlsx",
                           sheet_name='Parameters',
                           index_col=[0])
    r = Fparam.loc["InterestRate"]["Value"]
    
    # MAXIMUM HOURLY VARIATION OF FLEXIBLE DEMANDS 
    D_var = 0.1

    ########### FULL TIMESERIES OPERATION ONLY OPTIMIZATION MODEL #############
    
    array=[] # to save the optimal values of the decision variables
        
    def REC_fullts_operation(K, PV_res, EES_res, HP_res, TES_res,
                             BOIL_ter, BOIL_com,
                             PV_pub, EES_pub, HP_pub, TES_pub, years_test): 
    
        # Initialize the Gurobi model
        m = grb.Model()
        
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
        PPV_res = [[CapPV_res*SolRad[h,k] for k in range(K)] for h in range(H)] # SolRad in [kW/kWp]
        
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
        PPV_pub = [[CapPV_pub*SolRad[h,k] for k in range(K)] for h in range(H)] # SolRad in [kW/kWp]
            
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
        
        ### Warm start: Gurobi uses the clustering-based solutions 
        ### as an incumbent to guide the MIP search
        # Example: setting start values for technology capacities
        CapPV_res.Start   = PV_res
        CapEES_res.Start  = EES_res
        CapHP_res.Start   = HP_res
        CapTES_res.Start  = TES_res
        CapBOIL_ter.Start = BOIL_ter
        CapBOIL_com.Start = BOIL_com
        CapPV_pub.Start   = PV_pub
        CapEES_pub.Start  = EES_pub
        CapHP_pub.Start   = HP_pub
        CapTES_pub.Start  = TES_pub
        
        # Fix the lower bounds of capacities according to the "best" clustering-based design solution 
        # res
        m.addConstr(CapPV_res >= PV_res, name="LB_CapPV_res")
        m.addConstr(CapEES_res >= EES_res, name="LB_CapEES_res")
        m.addConstr(CapHP_res >= HP_res, name="LB_CapHP_res")
        m.addConstr(CapTES_res >= TES_res, name="LB_CapTES_res")
        # ter
        m.addConstr(CapBOIL_ter >= BOIL_ter, name="LB_CapBOIL_ter")
        # com
        m.addConstr(CapBOIL_com >= BOIL_com, name="LB_CapBOIL_com")
        # pub
        m.addConstr(CapPV_pub >= PV_pub, name="LB_CapPV_pub")
        m.addConstr(CapEES_pub >= EES_pub, name="LB_CapEES_pub")
        m.addConstr(CapHP_pub >= HP_pub, name="LB_CapHP_pub")
        m.addConstr(CapTES_pub >= TES_pub, name="LB_CapTES_pub")
        
        # Total net electrical power imported from (P_imp) and exported to (P_exp) the grid 
        m.addConstrs(P_imp[h,k] == P_imp_res[h,k]+P_imp_ter[h,k]+P_imp_com[h,k]+P_imp_pub[h,k] for h in range(H) for k in range(K))
        m.addConstrs(P_exp[h,k] == P_exp_res[h,k]+P_exp_pub[h,k] for h in range(H) for k in range(K))
        
        ### Calculation of electricity shared
        # m.addConstrs(Econd[h,k] == min_(P_imp[h,k],P_exp[h,k]) for h in range(H) for k in range(K))
        ### As alternative ###
        # Relaxation of the min constraint (it works if Econd is maximized in the obj function)
        for h in range(H):
            for k in range(K):
                m.addConstr(Econd[h,k] <= P_imp[h,k])
                m.addConstr(Econd[h,k] <= P_exp[h,k])
        
        # Fixed optimal capacities
        CapPV_tot = PV_res + PV_pub
        # Maximum possible incentive (already considering zonal prices)
        # equal for each hour
        if CapPV_tot <= 200:
            if loc=='Padova' or loc=='Trieste':
               inc = (130 + 10.57)/1000
            elif loc=='Palermo':
                 inc = (120 + 10.57)/1000
        elif CapPV_tot > 200 and CapPV_tot <= 600:
            if loc=='Padova' or loc=='Trieste':
               inc = (120 + 10.57)/1000
            elif loc=='Palermo':
                 inc = (110 + 10.57)/1000
        elif CapPV_tot > 600:
            if loc=='Padova' or loc=='Trieste':
               inc = (110 + 10.57)/1000
            elif loc=='Palermo':
                 inc = (100 + 10.57)/1000
        
        
        m.update()
        
        ### DEFINE COSTS
        
        ### The solution is based on full timeseries
        ### Time horizon of 6 years (whole testing dataset)    
        ### The investment cost is actualized on 1 year and reported to 6 years
        ### The operational cost is calculated in each day of the testing dataset
        
        # Investment cost [€]
        CostINV = years_test*((r*np.power(1+r,nPV)/(np.power(1+r,nPV)-1)+OnM_fix_PV)*invPV*(CapPV_res+CapPV_pub) +\
                  (r*np.power(1+r,nBOIL)/(np.power(1+r,nBOIL)-1)+OnM_fix_BOIL)*invBOIL*(CapBOIL_ter+CapBOIL_com) +\
                  (r*np.power(1+r,nHP)/(np.power(1+r,nHP)-1)+OnM_fix_HP)*invHP*(CapHP_res+CapHP_pub) +\
                  (r*np.power(1+r,nTES)/(np.power(1+r,nTES)-1)+OnM_fix_TES)*invTES*(CapTES_res+CapTES_pub) +\
                  (r*np.power(1+r,nEES)/(np.power(1+r,nEES)-1)+OnM_fix_EES)*invEES*(CapEES_res+CapEES_pub))
        
        # Operational costs [€], including the REC incentive
        ### With maximum possible incentive
        CostOP = sum(sum((FBOIL_ter[h][k]+FBOIL_com[h][k])*cbuy_ng 
                         + P_imp[h,k]*cbuy_el[h,k] - P_exp[h,k]*csell_el[h,k] - Econd[h,k]*inc
                         for h in range(H)) for k in range(K))
        
        m.update()
            
        ### SET THE OBJECTIVE FUNCTION
        
        m.setObjective(CostINV+CostOP, grb.GRB.MINIMIZE)    
        
        opt_gap = []
        def data_cb(m, where):
            if where == grb.GRB.Callback.MIPSOL:
               cur_obj = m.cbGet(grb.GRB.Callback.MIPSOL_OBJBST)  # best incumbent solution (upper bound, primal solution)
               cur_bd = m.cbGet(grb.GRB.Callback.MIPSOL_OBJBND)   # best bound solution (lower bound, dual solution)
               opt_gap.append(((cur_obj-cur_bd)/cur_obj)*100)     # all gaps during iterations
        
        ### Solver parameters
        m.Params.NonConvex = 2         # non-convex MIQCP model
        
        ### Padova
        
        ### Trials for opt_K = 10
        
        ### gap of 43.9%
        # m.Params.MIPGap = 2e-2         # tolerance of the optimal solution in []
        # m.Params.MIPFocus = 1          # prioritize feasibility over optimality (2), default is a trade-off (0)
        # m.Params.IntegralityFocus = 1  # emphasize integrality
        # m.Params.Heuristics = 0.2      # spend 20% of time on heuristics
        
        ### gap of 59.7%
        # m.Params.MIPGap = 10e-2        # tolerance of the optimal solution in []
        # m.Params.MIPFocus = 3          # focus on improving bound for hard models
        # m.Params.Heuristics = 0.05     # back to default for speed
        # m.Params.Cuts = 2              # aggressive cuts
        
        ### So far, the best gap of 11.57%
        # m.Params.MIPGap = 10e-2          # tolerance of the optimal solution in []
        # m.Params.MIPFocus = 1            # focus on improving bound
        # m.Params.Heuristics = 0.05       # back to default for speed
        # m.Params.Cuts = 2                # aggressive cuts
        
        ### So far, the best: gap of 11.57%
        # m.Params.MIPGap = 5e-2           # tolerance of the optimal solution in []
        # m.Params.MIPFocus = 1            # focus on improving bound
        # m.Params.Heuristics = 0.05       # back to default for speed
        # m.Params.Cuts = 2                # aggressive cuts
        
        ### So far, the best: gap of 11.57%
        # m.Params.MIPGap = 5e-2           # tolerance of the optimal solution in []
        # m.Params.MIPFocus = 1            # focus on improving bound
        # m.Params.Heuristics = 0.05       # back to default for speed
        # m.Params.Cuts = 2                # aggressive cuts
        # m.Params.IntegralityFocus = 1    # emphasize integrality
        
        ### Trials for opt_K = 15
        
        ### gap of 45.62%
        # m.Params.MIPGap = 5e-2           # tolerance of the optimal solution in []
        # m.Params.MIPFocus = 1            # focus on improving bound
        # m.Params.Heuristics = 0.05       # back to default for speed
        # m.Params.Cuts = 2                # aggressive cuts
        
        ### Trials for opt_K = 20
        
        ### gap of 81.59%
        # m.Params.MIPGap = 5e-2           # tolerance of the optimal solution in []
        # m.Params.MIPFocus = 1            # focus on improving bound
        # m.Params.Heuristics = 0.05       # back to default for speed
        # m.Params.Cuts = 2                # aggressive cuts
        
        ### Trials for opt_K = 30: no solution
        
        # m.Params.MIPGap = 5e-2           # tolerance of the optimal solution in []
        # m.Params.MIPFocus = 1            # focus on improving bound
        # m.Params.Heuristics = 0.05       # back to default for speed
        # m.Params.Cuts = 2                # aggressive cuts
        
        
        ### Palermo
        
        ### Trials for opt_K = 10: no solution
        
        ### Trials for opt_K = 15
        
        ### Gap of 11.7%
        # m.Params.MIPGap = 5e-2           # tolerance of the optimal solution in []
        # m.Params.MIPFocus = 1            # focus on improving bound
        # m.Params.Heuristics = 0.05       # back to default for speed
        # m.Params.Cuts = 2                # aggressive cuts
        
        ### Trials for opt_K = 20
        
        ### Gap of 81.93%
        m.Params.MIPGap = 5e-2           # tolerance of the optimal solution in []
        m.Params.MIPFocus = 1            # focus on improving bound
        m.Params.Heuristics = 0.05       # back to default for speed
        m.Params.Cuts = 2                # aggressive cuts
        
        m.optimize(callback=data_cb)   # optimization with callback to retrieve optimal gaps
        
        # Check model status
        if m.status == grb.GRB.INFEASIBLE:
            print("Model is infeasible. Computing IIS...")
            m.computeIIS()
            m.write("model.ilp")
            m.write("model.ilp.json")
            print("IIS written to model.ilp and model.ilp.json")
            return None  # exit gracefully, don't try to access objVal
        
        elif m.status == grb.GRB.UNBOUNDED:
            print("Model is unbounded. You may need to add bounds or constraints.")
            return None
        
        # elif m.status == grb.GRB.OPTIMAL:
        #     # Assess the objectives
        #     C_TOT = m.objVal  # optimal value of the objective function [€]
        #     design_cost = CostINV.getValue()
        #     operation_cost = CostOP.getValue()
            
        #     for v in m.getVars():
        #         array.append(v.x)
        
        if len(opt_gap) == 0:  # Handle empty case
           print("Warning: Callback did not run, retrieving gap from solver.")
           if m.status == grb.GRB.OPTIMAL:
              opt_gap = m.Params.MIPGap * 100  # Get gap directly from solver
           else:
              opt_gap = np.nan  # Assign NaN if no solution was found
        else:
              opt_gap = np.array(opt_gap)
              print(f"opt_gap contents: {opt_gap}")
              opt_gap = opt_gap[-1]  # Get final gap
    
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
        # array_TIP=[]
        # array_inc=[]
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
                # array_TIP.append(TIP[h,k].X)
                # array_inc.append(inc[h,k].X)
        
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
        # array_TIP=np.array(array_TIP)
        # array_inc=np.array(array_inc)
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
        # array_TIP=np.reshape(array_TIP,(H,K),order='C')
        # array_inc=np.reshape(array_inc,(H,K),order='C')
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
        
        # return  (C_TOT, design_cost, operation_cost, 
        # array_P_imp, array_P_exp, 
        # array_Econd, PZ, array_TIP, array_inc, 
        # cap_PV_res, cap_EES_res, cap_HP_res, cap_TES_res, 
        # cap_BOIL_ter, cap_BOIL_com, 
        # cap_PV_pub, cap_EES_pub, cap_HP_pub, cap_TES_pub, 
        # opt_gap, csell_el, cbuy_el, 
        # array_DemEl_res_op, array_DemEl_ter_op, array_DemEl_com_op, array_DemEl_pub_op)
        
        return  (C_TOT, design_cost, operation_cost, 
        array_P_imp, array_P_exp, 
        array_Econd, PZ, 
        cap_PV_res, cap_EES_res, cap_HP_res, cap_TES_res, 
        cap_BOIL_ter, cap_BOIL_com, 
        cap_PV_pub, cap_EES_pub, cap_HP_pub, cap_TES_pub, 
        opt_gap, csell_el, cbuy_el, 
        array_DemEl_res_op, array_DemEl_ter_op, array_DemEl_com_op, array_DemEl_pub_op)

    ##### SOLVING THE FULL TIMESERIES OPERATION-ONLY OPTIMIZATION MODEL #######
    ##### IN THE TESTING DATASET
    
    cfold = os.getcwd()  # working directory
    
    ### Import the full testing dataset (2018-2023)
    clust_fold = os.path.join(cfold, 'Clustering')   # folder with clustering results
    PATH_clustering = clust_fold+"/"+loc
    timeseries = np.load(os.path.join(PATH_clustering,'Clustering_train_'+extreme_method+'.npz'), allow_pickle=True)
    
    years_test = timeseries['years_test']
    test_c = timeseries['test_c'] # testing dataset with all timeseries over years (2D)
    K = len(test_c)   # number of days in the testing dataset
    
    ### Import the optimal sizes of the "best" clustering-based design solution
    # Clustering-based design solutions
    results_fold = os.path.join(cfold, 'Results')
    PATH_results = results_fold+"/"+loc
    if seasons==0:
       opt_sol = np.load(os.path.join(PATH_results,'Stochastic_solutions_annual_'+extreme_method+'_'+loc+'.npz'), allow_pickle=True)
    else:
       opt_sol = np.load(os.path.join(PATH_results,'Stochastic_solutions_seasons_'+extreme_method+'_'+loc+'.npz'), allow_pickle=True)
    
    # Optimal sizes
    if extreme_method=='A':
        PV_res = opt_sol['capPV_res_c'][opt_K-4]     # "opt_K" is the "best" set of typical days
        EES_res = opt_sol['capEES_res_c'][opt_K-4]
        HP_res = opt_sol['capHP_res_c'][opt_K-4]
        TES_res = opt_sol['capTES_res_c'][opt_K-4]
        BOIL_ter = opt_sol['capBOIL_ter_c'][opt_K-4]
        BOIL_com = opt_sol['capBOIL_com_c'][opt_K-4]
        PV_pub = opt_sol['capPV_pub_c'][opt_K-4]
        EES_pub = opt_sol['capEES_pub_c'][opt_K-4]
        HP_pub = opt_sol['capHP_pub_c'][opt_K-4]
        TES_pub = opt_sol['capTES_pub_c'][opt_K-4]
    elif extreme_method=='R':
         PV_res = opt_sol['capPV_res_c'][opt_K-2]
         EES_res = opt_sol['capEES_res_c'][opt_K-2]
         HP_res = opt_sol['capHP_res_c'][opt_K-2]
         TES_res = opt_sol['capTES_res_c'][opt_K-2]
         BOIL_ter = opt_sol['capBOIL_ter_c'][opt_K-2]
         BOIL_com = opt_sol['capBOIL_com_c'][opt_K-2]
         PV_pub = opt_sol['capPV_pub_c'][opt_K-2]
         EES_pub = opt_sol['capEES_pub_c'][opt_K-2]
         HP_pub = opt_sol['capHP_pub_c'][opt_K-2]
         TES_pub = opt_sol['capTES_pub_c'][opt_K-2]
    
    start_time = time.time()
    
    # Length of a day
    H = 24
    
    if loc=='Padova' or loc=='Trieste':
       c_price = 0.5
    elif loc=='Palermo':
         c_price = 0.45
    
    SolRad = np.zeros((H,K))       # PV capacity factor
    Tamb = np.zeros((H,K))         # ambient temperature
    csell_el = np.zeros((H,K))     # zonal prices
    DemEl_res = np.zeros((H,K))    # electricity demand of Res
    DemEl_ter = np.zeros((H,K))    # electricity demand of Ter
    DemEl_com = np.zeros((H,K))    # electricity demand of Com
    DemEl_pub = np.zeros((H,K))    # electricity demand of Pub
    DemTh_res = np.zeros((H,K))    # heating demand of Res 
    DemTh_ter = np.zeros((H,K))    # heating demand of Ter
    DemTh_com = np.zeros((H,K))    # heating demand of Com
    DemTh_pub = np.zeros((H,K))    # heating demand of Pub
    PZ = np.zeros((H,K))           # zonal market prices
    for yi in range(K):            # different days
        SolRad[:,yi] = (test_c[yi,:][0:24])/1000 # [kW/kWp]
        Tamb[:,yi] = test_c[yi,:][24:48]
        csell_el[:,yi] = (test_c[yi,:][48:72])/1000*c_price
        DemEl_res[:,yi] = test_c[yi,:][72:96]
        DemEl_ter[:,yi] = test_c[yi,:][96:120]
        DemEl_com[:,yi] = test_c[yi,:][120:144]
        DemEl_pub[:,yi] = test_c[yi,:][144:168]
        DemTh_res[:,yi] = test_c[yi,:][168:192]
        DemTh_ter[:,yi] = test_c[yi,:][192:216]
        DemTh_com[:,yi] = test_c[yi,:][216:240]
        DemTh_pub[:,yi] = test_c[yi,:][240:264]
        PZ[:,yi] = (test_c[yi,:][48:72])  # €/MWh
    cbuy_el = csell_el + 0.20 # purchase price, 0.20 €/kWh is the grid tariff
    
    # Solve the full timeseries operation-only optimization
    
    # C_TOT_sto_c, design_cost_c, operation_cost_c, \
    # P_imp_c, P_exp_c, \
    # E_cond_c, PZ_c, TIP_c, inc_c, \
    # capPV_res_c, capEES_res_c, capHP_res_c, capTES_res_c, \
    # capBOIL_ter_c, capBOIL_com_c, \
    # capPV_pub_c, capEES_pub_c, capHP_pub_c, capTES_pub_c, \
    # opt_gap_c, c_sell_save_c, c_buy_save_c, \
    # DemEl_res_op_c, DemEl_ter_op_c, DemEl_com_op_c, DemEl_pub_op_c = REC_fullts_operation(K, PV_res, EES_res, HP_res, TES_res,
    #                                                                                       BOIL_ter, BOIL_com,
    #                                                                                       PV_pub, EES_pub, HP_pub, TES_pub, years_test)
    
    C_TOT_sto_c, design_cost_c, operation_cost_c, \
    P_imp_c, P_exp_c, \
    E_cond_c, PZ_c, \
    capPV_res_c, capEES_res_c, capHP_res_c, capTES_res_c, \
    capBOIL_ter_c, capBOIL_com_c, \
    capPV_pub_c, capEES_pub_c, capHP_pub_c, capTES_pub_c, \
    opt_gap_c, c_sell_save_c, c_buy_save_c, \
    DemEl_res_op_c, DemEl_ter_op_c, DemEl_com_op_c, DemEl_pub_op_c = REC_fullts_operation(K, PV_res, EES_res, HP_res, TES_res,
                                                                                          BOIL_ter, BOIL_com,
                                                                                          PV_pub, EES_pub, HP_pub, TES_pub, years_test)
    
    ### TIME TO SOLVE THE MODEL ###
    print("--- %s seconds ---" % (time.time() - start_time))

    ##################### SAVE AND EXPORT RESULTS OF TESTING ##################
    
    cfold = os.getcwd()  # working directory
    results_fold = os.path.join(cfold, 'Results')
    PATH_results = results_fold+"/"+loc
    
    # # Save results
    # np.savez(os.path.join(PATH_results,f"Stochastic_solution_{opt_K}_tested_"+extreme_method+'_'+loc), 
    #           C_TOT_sto_c = C_TOT_sto_c, design_cost_c = design_cost_c, operation_cost_c = operation_cost_c,
    #           E_cond_c = E_cond_c, PZ_c = PZ_c, TIP_c = TIP_c, inc_c = inc_c,
    #           capPV_res_c = capPV_res_c, capEES_res_c = capEES_res_c, capHP_res_c = capHP_res_c, capTES_res_c = capTES_res_c, 
    #           capBOIL_ter_c = capBOIL_ter_c, capBOIL_com_c = capBOIL_com_c, 
    #           capPV_pub_c = capPV_pub_c, capEES_pub_c = capEES_pub_c, capHP_pub_c = capHP_pub_c, capTES_pub_c = capTES_pub_c, 
    #           opt_gap_c = opt_gap_c, c_sell_save_c = c_sell_save_c, c_buy_save_c = c_buy_save_c,
    #           DemEl_res_op_c = DemEl_res_op_c, DemEl_ter_op_c = DemEl_ter_op_c,
    #           DemEl_com_op_c = DemEl_com_op_c, DemEl_pub_op_c = DemEl_pub_op_c)
    
    # Save results
    np.savez(os.path.join(PATH_results,f"Stochastic_solution_{opt_K}_tested_"+extreme_method+'_'+loc), 
              C_TOT_sto_c = C_TOT_sto_c, design_cost_c = design_cost_c, operation_cost_c = operation_cost_c,
              E_cond_c = E_cond_c, PZ_c = PZ_c,
              capPV_res_c = capPV_res_c, capEES_res_c = capEES_res_c, capHP_res_c = capHP_res_c, capTES_res_c = capTES_res_c, 
              capBOIL_ter_c = capBOIL_ter_c, capBOIL_com_c = capBOIL_com_c, 
              capPV_pub_c = capPV_pub_c, capEES_pub_c = capEES_pub_c, capHP_pub_c = capHP_pub_c, capTES_pub_c = capTES_pub_c, 
              opt_gap_c = opt_gap_c, c_sell_save_c = c_sell_save_c, c_buy_save_c = c_buy_save_c,
              DemEl_res_op_c = DemEl_res_op_c, DemEl_ter_op_c = DemEl_ter_op_c,
              DemEl_com_op_c = DemEl_com_op_c, DemEl_pub_op_c = DemEl_pub_op_c)

    # Export optimal capacities
    PATH_export_sto = PATH_results+"/"+f"Stochastic_solution_{opt_K}_tested_"+extreme_method+'_'+loc+".xlsx"
         
    df1=pd.DataFrame(C_TOT_sto_c/1000)
    df2=pd.DataFrame(design_cost_c/1000)
    df3=pd.DataFrame(operation_cost_c/1000)
    df4=pd.DataFrame(data=np.array([capPV_res_c, capEES_res_c, capHP_res_c, capTES_res_c,
                                    capBOIL_ter_c, capBOIL_com_c,
                                    capPV_pub_c, capEES_pub_c, capHP_pub_c, capTES_pub_c]).T,
                                    columns=['capPV_res [kWp]', 'capEES_res [kWh]', 'capHP_res [kW]', 'capTES_res [kWh]',
                                             'capBOIL_ter [kW]', 'capBOIL_com [kW]',
                                             'capPV_pub [kWp]', 'capEES_pub [kWh]', 'capHP_pub [kW]', 'capTES_pub [kWh]'])
    
    writer = pd.ExcelWriter(PATH_export_sto, engine="xlsxwriter")   
    df1.to_excel(writer, sheet_name='Total_costs', float_format="%.13f")
    df2.to_excel(writer, sheet_name='Design_costs', float_format="%.13f")
    df3.to_excel(writer, sheet_name='Operation_costs', float_format="%.13f")
    df4.to_excel(writer, sheet_name='Optimal sizes', float_format="%.13f")
    
    writer.close()