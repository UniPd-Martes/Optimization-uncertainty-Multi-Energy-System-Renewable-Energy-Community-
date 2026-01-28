# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 15:25:34 2024

@author: utente
"""

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# TWO FUNCTIONS
# "scenario_generation": to obtain different sets of typical days representative 
# of the full timeseries in the training dataset (2005-2017)
# "check_clustering": post-processing of clustering results

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# SCENARIO GENERATION

# Inputs:
# - loc: name of location for weather variables (solar irradiance, ambient temperature) 
#   and zonal electricity prices (type: string, "Padova", "Palermo", "Trieste")
# - year_start: first year of the available dataset of timeseries (type: int, default: 2005)
# - year_end: last year of the available dataset of timeseries (type: int, default: 2023)
# - seasons: parameter deciding whether to consider or not division of the
#   available dataset into seasons (type: int, 0: no, 1: yes, default: 0)
# - extreme_method: criterion to take into account extreme days
#   (type: string, "A": adding criterion, "R": replacing criterion, default: "A")
# - plot_ts: plot full timeseries of solar irradiance, ambient temperature and zonal electricity prices
#   (type: string, "yes" or "no")

def scenario_generation(loc="",
                        year_start=2005,
                        year_end=2023,
                        seasons=0,
                        extreme_method='A', 
                        plot_ts=""):


    ### Python packages
    import pandas as pd                              # Pandas, for data analysis 
    import numpy as np                               # Numpy, for multi-dimensional array
    import os                                        # to select the folder for saving
    # Clustering
    from sklearn.cluster import KMeans               # clustering algorithm 
    from sklearn import metrics                      # clustering metrics
    from sklearn.preprocessing import MinMaxScaler   # data normalization before clustering
    from scipy.spatial.distance import cdist         # distance and similarity metrics
    from scipy.signal import find_peaks              # analyze, filter, and transform timeseries       
    # Data visualization and plots
    import matplotlib.pyplot as plt                  # Pyplot, for plots, graphs, diagrams
    import matplotlib.dates as mdates               
    import time                                      # time required for a simulation
    
    ### The stochastic model is based on typical days (clustering-based design model)
    ### Multi-year clutering
    ### - Annual or seasonal clustering 
    ### - Replacing or adding approaches for extreme days
    ### Clustering attributes are: PV capacity factor, ambient temperature, zonal electricity price
    
    ### Input automatically defined: 
    ### - Uncertain parameters, number and profiles of stochastic scenarios of the uncertain parameters
    
    cfold = os.getcwd()  # working directory
    
    ### Training dataset: 2005 ---> 2017
    ### Testing dataset: 2018 ---> 2023
    years_train=13   # number of years in the training dataset
    years_test=6     # number of years in the testing dataset
    years = years_train + years_test
    
    ############################# DATA PROCESSING #############################
    
    ifold = os.path.join(cfold, 'Input_data')   # folder with input data
    
    ### Solar irradiance from PVGIS 5.3: radiation components and PV capacity factor
    DF_rad = pd.read_csv(ifold+"/"+loc+"_PV.csv", skiprows=10, nrows=166536)
    DF_rad.drop(columns=['time']) 
    
    ### Data of PUN (Prezzo Unico Nazionale) and zonal prices
    # The data are imported by a for cycle (data of each year in a separate file)
    DF_PUN=[]   # all years from 2005 to 2023
    for Year in range (year_start,year_end+1): 
        PATH_ts_PUN = ifold+"/Prices/Anno "+str(Year)+".xlsx"
        DF_PUN_Year = pd.read_excel(PATH_ts_PUN, sheet_name='Prezzi-Prices') #specific year
        DF_PUN.append(DF_PUN_Year)
    DF_PUN=pd.concat(DF_PUN, ignore_index=True)
    
    ### Hourly values of electricity and heating demands referring to one year
    PATH_ts_Dem = ifold+"/Demands/EC_HourlyYearValues.xlsx"
    # One year dataset of electricity and heating demands
    DemRes_EL = np.array(pd.read_excel(PATH_ts_Dem,sheet_name='NY_ResidentialLow',usecols='A')).flatten() \
                + np.array(pd.read_excel(PATH_ts_Dem,sheet_name='NY_ResidentialBase',usecols='A')).flatten() \
                + np.array(pd.read_excel(PATH_ts_Dem,sheet_name='NY_ResidentialHigh',usecols='A')).flatten()
    DemTer1_EL = np.array(pd.read_excel(PATH_ts_Dem,sheet_name='NY_SmallOffice',usecols='A')).flatten() 
    DemCom1_EL = np.array(pd.read_excel(PATH_ts_Dem,sheet_name='NY_QuickServiceRestaurant',usecols='A')).flatten() 
    DemPub1_EL = np.array(pd.read_excel(PATH_ts_Dem,sheet_name='NY_PrimarySchool',usecols='A')).flatten() 
    DemRes_TH = np.array(pd.read_excel(PATH_ts_Dem,sheet_name='NY_ResidentialLow',usecols='B')).flatten() \
                + np.array(pd.read_excel(PATH_ts_Dem,sheet_name='NY_ResidentialBase',usecols='B')).flatten() \
                + np.array(pd.read_excel(PATH_ts_Dem,sheet_name='NY_ResidentialHigh',usecols='B')).flatten()
    DemTer1_TH = np.array(pd.read_excel(PATH_ts_Dem,sheet_name='NY_SmallOffice',usecols='B')).flatten() 
    DemCom1_TH = np.array(pd.read_excel(PATH_ts_Dem,sheet_name='NY_QuickServiceRestaurant',usecols='B')).flatten() 
    DemPub1_TH = np.array(pd.read_excel(PATH_ts_Dem,sheet_name='NY_PrimarySchool',usecols='B')).flatten() 
    # The hourly values of the electricity and heating demands refer to one year
    # Values of the single year are repeated for all years of the dataset 2005-2023
    DemRes_el = np.concatenate([DemRes_EL]*19)
    DemTer_el = np.concatenate([DemTer1_EL]*19)
    DemCom_el = np.concatenate([DemCom1_EL]*19)
    DemPub_el = np.concatenate([DemPub1_EL]*19)
    DemRes_th = np.concatenate([DemRes_TH]*19)
    DemTer_th = np.concatenate([DemTer1_TH]*19)
    DemCom_th = np.concatenate([DemCom1_TH]*19)
    DemPub_th = np.concatenate([DemPub1_TH]*19)
    
    # Delete column "time" of "solar" dataframe and divide it into year, month, day, hour and date 
    TIME = DF_rad['time'] 
    TIME = TIME.tolist() 
    Date = []   
    Year = []   
    Month = []  
    Day = []    
    Hour = []    
    for i in range(len(TIME)):
        Year.append(TIME[i][:4])
        Month.append(TIME[i][4:6])
        Day.append(TIME[i][6:8])
        Hour.append(TIME[i][9:11])
        ymdh = TIME[i][:8]+TIME[i][9:11]
        ymdh = ymdh[:4] + "-" + ymdh[4:]
        ymdh = ymdh[:7] + "-" + ymdh[7:]
        ymdh = ymdh[:10] + " " + ymdh[10:]
        Date.append(ymdh)
    # Delete columns
    DF_rad.drop('time', inplace=True, axis=1) # delete column "time"
    # New columns
    Date = pd.date_range(Date[0], periods=len(TIME), freq="H")
    DF_rad.insert(0, "Date", Date, True)
    DF_rad.insert(1, "Year", Year, True)
    DF_rad.insert(2, "Month", Month, True)
    DF_rad.insert(3, "Day", Day, True)
    DF_rad.insert(4, "Hour", Hour, True)
    DF_PUN.insert(0, "Date", Date, True)
    DF_PUN.insert(1, "Year", Year, True)
    DF_PUN.insert(2, "Month", Month, True)
    DF_PUN.insert(3, "Day", Day, True)
    DF_PUN.insert(4, "Hour", Hour, True)
    
    # Remove 29/02 days
    to_rem = DF_rad[((DF_rad.Month == '02') & (DF_rad.Day == '29'))].index   # DF_rad
    DF_rad.drop(to_rem, inplace=True)
    to_remm = DF_PUN[((DF_PUN.Month == '02') & (DF_PUN.Day == '29'))].index   # DF_PUN
    DF_PUN.drop(to_remm, inplace=True)
    
    ### Adding zonal prices, electricity and heating demands to the whole dataset
    DF = DF_rad
    if loc=='Padova' or loc=='Trieste':
       DF['PZ'] = DF_PUN["NORD"].values  # zonal prices
    elif loc=='Palermo':
         DF['PZ'] = DF_PUN["SICI"].values
    ### Add here another "if" statement to define the zonal price timeseries in other locations
    
    DF['Res_el'] = DemRes_el
    DF['Ter_el'] = DemTer_el
    DF['Com_el'] = DemCom_el
    DF['Pub_el'] = DemPub_el
    DF['Res_th'] = DemRes_th
    DF['Ter_th'] = DemTer_th
    DF['Com_th'] = DemCom_th
    DF['Pub_th'] = DemPub_th
    
    DF.set_index("Year", inplace = True)
    
    ####################### SEASONAL DATASETS FOR EACH YEAR ###################
    
    # Year as index
    new_cols=["Date","Year","Month","Day","Hour","P","T2m","PZ","Res_el","Ter_el","Com_el","Pub_el","Res_th","Ter_th","Com_th","Pub_th"] #new columns order
    
    # DIVISION BY SEASONS
    DF_winter=[]
    DF_spring=[]
    DF_summer=[]
    DF_autumn=[]
    
    for year in range(year_start,year_end+1):     
            
        DF_train_test = DF.loc[str(year)]
        DF_train_test.insert(0, "Year", [year]*len(DF_train_test), True)    
        DF_train_test.set_index(['Month'], inplace = True)
        pd.options.mode.chained_assignment = None # to avoid warnings about view or copy of a dataframe
        DF_winter_Year=DF_train_test.loc['01':'03']
        DF_winter_Year.reset_index(inplace=True) 
        DF_spring_Year=DF_train_test.loc['04':'06']
        DF_spring_Year.reset_index(inplace=True) 
        DF_summer_year=DF_train_test.loc['07':'09']
        DF_summer_year.reset_index(inplace=True) 
        DF_autumn_Year=DF_train_test.loc['10':'12']
        DF_autumn_Year.reset_index(inplace=True) 
        
        ### Winter ###
        #remove the last days of March for Winter
        rem1 = DF_winter_Year[((DF_winter_Year.Month == '03') & (DF_winter_Year.Day >='21'))].index
        DF_winter_Year.drop(rem1, inplace=True)
        #add the last days of December for Winter
        DF_train_test.reset_index(inplace=True) 
        end_dec_Year=DF_train_test[((DF_train_test.Month == '12') & (DF_train_test.Day >='21'))]
        DF_winter_df=pd.concat([end_dec_Year,DF_winter_Year],ignore_index=True)
        DF_winter_df=DF_winter_df.reindex(columns=new_cols)
        DF_winter.append(DF_winter_df)
        
        ### Spring ###
        #remove the last days of June for Spring
        rem11 = DF_spring_Year[((DF_spring_Year.Month == '06') & (DF_spring_Year.Day >='21'))].index
        DF_spring_Year.drop(rem11, inplace=True)
        #add the last days of March for Spring
        end_mar_Year=DF_train_test[((DF_train_test.Month == '03') & (DF_train_test.Day >='21'))]
        DF_spring_df=pd.concat([end_mar_Year,DF_spring_Year],ignore_index=True)
        DF_spring_df=DF_spring_df.reindex(columns=new_cols)
        DF_spring.append(DF_spring_df)
        
        ### Summer ###
        #remove the last days of September for Summer
        rem111 = DF_summer_year[((DF_summer_year.Month == '09') & (DF_summer_year.Day >='23'))].index
        DF_summer_year.drop(rem111, inplace=True)
        #add the last days of June for Summer
        end_jun_Year=DF_train_test[((DF_train_test.Month == '06') & (DF_train_test.Day >='21'))]
        DF_summer_df=pd.concat([end_jun_Year,DF_summer_year],ignore_index=True)
        DF_summer_df=DF_summer_df.reindex(columns=new_cols)
        DF_summer.append(DF_summer_df)
        
        ### Autumn ###
        #remove the last days of December for Autumn
        rem1111 = DF_autumn_Year[((DF_autumn_Year.Month == '12') & (DF_autumn_Year.Day >='21'))].index
        DF_autumn_Year.drop(rem1111, inplace=True)
        #add the last days of September for Autumn
        end_sep_Year=DF_train_test[((DF_train_test.Month == '09') & (DF_train_test.Day >='23'))]
        DF_autumn_df=pd.concat([end_sep_Year,DF_autumn_Year],ignore_index=True)
        DF_autumn_df=DF_autumn_df.reindex(columns=new_cols)
        DF_autumn.append(DF_autumn_df)
    
    # NO DIVISION BY SEASONS
    DF_no_s=[]
    for year in range(year_start,year_end+1):
        DF_train_test = DF.loc[str(year)]
        DF_train_test.insert(0, "Year", [year]*len(DF_train_test), True)
        DF_train_test=DF_train_test.reindex(columns=new_cols)
        DF_no_s.append(DF_train_test)
        
    # Mean daily profiles for each season and year
    DF_winter_mean=[]
    DF_spring_mean=[]
    DF_summer_mean=[]
    DF_autumn_mean=[]
    DF_no_s_mean=[]
    for year in range(year_start,year_end+1):
        DF_winter_mean.append(DF_winter[year-year_start].groupby(["Hour"]).mean(numeric_only=True))  # average only in columns with numeric values
        DF_spring_mean.append(DF_spring[year-year_start].groupby(["Hour"]).mean(numeric_only=True))
        DF_summer_mean.append(DF_summer[year-year_start].groupby(["Hour"]).mean(numeric_only=True))
        DF_autumn_mean.append(DF_autumn[year-year_start].groupby(["Hour"]).mean(numeric_only=True))
        DF_no_s_mean.append(DF_no_s[year-year_start].groupby(["Hour"]).mean(numeric_only=True))
    
    # total number of days for season
    N_days_winter=int(len(DF_winter[0])/24) 
    N_days_spring=int(len(DF_spring[0])/24) 
    N_days_summer=int(len(DF_summer[0])/24) 
    N_days_autumn=int(len(DF_autumn[0])/24) 
    
    ################ TRAINING (train) and TESTING (test) datasets #############
    
    ### Training dataset: 2005 ---> 2017
    ### Testing dataset: 2018 ---> 2023
    
    # The possible uncertain parameters are:
    # "P","T2m","PZ","Res_el","Ter_el","Com_el","Pub_el","Res_th","Ter_th","Com_th","Pub_th"
    # print('What uncertain parameters do you want to consider?')
    # U=input("Enter names of uncertain parameters separated by spaces: ") # string with spaces
    
    U="P T2m PZ"                # PV capacity factor, temperature and zonal prices are the uncertain parameters
    U=U.split()                 # each word of the string becomes an item of a list
    U={u: None for u in U}      # the names of the uncertain parameters become keys of a dictionary
    keys=list(U.keys())         # keys of the dictionary as list
    
    # print("Do you want to divide the datasets by seasons? 0: no, 1: yes")
    # seasons = int(input())
    
    if seasons==0: ### WITHOUT DIVISION BY SEASONS
    
        DF_train = DF.loc['2005':'2017'] 
        DF_test = DF.loc['2018':'2023']
        DF_train_test = DF.loc['2005':'2023']
    
        train = []
        train_c= []
        test = []
        test_c= []
        train_test = []
        train_test_c = []
        
        x=0
        while x in range(len(DF_train)):
                    
              train_df = np.concatenate([DF_train[key].values[x:x+24] for key in keys])
              train.append(train_df)
                    
              train_df_c = np.concatenate([DF_train['P'].values[x:x+24],DF_train['T2m'].values[x:x+24],DF_train['PZ'].values[x:x+24],
                                           DF_train['Res_el'].values[x:x+24],DF_train['Ter_el'].values[x:x+24],
                                           DF_train['Com_el'].values[x:x+24],DF_train['Pub_el'].values[x:x+24],
                                           DF_train['Res_th'].values[x:x+24],DF_train['Ter_th'].values[x:x+24],
                                           DF_train['Com_th'].values[x:x+24],DF_train['Pub_th'].values[x:x+24]])  
              train_c.append(train_df_c)
            
              x=x+24
        
        x=0
        while x in range(len(DF_test)):
                    
              test_df = np.concatenate([DF_test[key].values[x:x+24] for key in keys])
              test.append(test_df)
                    
              test_df_c = np.concatenate([DF_test['P'].values[x:x+24],DF_test['T2m'].values[x:x+24],DF_test['PZ'].values[x:x+24],
                                           DF_test['Res_el'].values[x:x+24],DF_test['Ter_el'].values[x:x+24],
                                           DF_test['Com_el'].values[x:x+24],DF_test['Pub_el'].values[x:x+24],
                                           DF_test['Res_th'].values[x:x+24],DF_test['Ter_th'].values[x:x+24],
                                           DF_test['Com_th'].values[x:x+24],DF_test['Pub_th'].values[x:x+24]])  
              test_c.append(test_df_c)
            
              x=x+24     
        
        x=0
        while x in range(len(DF_train_test)):
                    
              train_test_df = np.concatenate([DF_train_test[key].values[x:x+24] for key in keys])
              train_test.append(train_test_df)
                    
              train_test_df_c = np.concatenate([DF_train_test['P'].values[x:x+24],DF_train_test['T2m'].values[x:x+24],DF_train_test['PZ'].values[x:x+24],
                                                DF_train_test['Res_el'].values[x:x+24],DF_train_test['Ter_el'].values[x:x+24],
                                                DF_train_test['Com_el'].values[x:x+24],DF_train_test['Pub_el'].values[x:x+24],
                                                DF_train_test['Res_th'].values[x:x+24],DF_train_test['Ter_th'].values[x:x+24],
                                                DF_train_test['Com_th'].values[x:x+24],DF_train_test['Pub_th'].values[x:x+24]])  
              train_test_c.append(train_test_df_c)
            
              x=x+24
              
        train = np.array(train)
        train_c = np.array(train_c)
        test = np.array(test)
        test_c = np.array(test_c)
        train_test = np.array(train_test)
        train_test_c = np.array(train_test_c)
        
        # Number of columns for the complete annual dataset
        lenghth_features_c = len(train_c[1])
        
        ### Training and testing datasets divided by years
        train_c_years = np.zeros((int(len(train_c)/years_train),lenghth_features_c,years_train))
        test_c_years = np.zeros((int(len(test_c)/years_test),lenghth_features_c,years_test))
        year=0
        while year in range(0,years_train):
              train_c_years[:,:,year] = train_c[year*365:year*365+365,:]
              year = year+1
        year=0
        while year in range(0,years_test):
              test_c_years[:,:,year] = test_c[year*365:year*365+365,:]
              year = year+1
              
    elif seasons==1: ### WITH DIVISION BY SEASONS
    
        # Output
        # seasonal training/testing datasets including all years (3D)
        train_test_w = np.zeros((N_days_winter,len(keys)*24,years))  # with only uncertain parameters (3 timeseries) for 19 years (training and testing)
        train_test_c_w = np.zeros((N_days_winter,11*24,years))       # with all 11 timeseries for 19 years (training and testing)
        train_test_sp = np.zeros((N_days_spring,len(keys)*24,years))  
        train_test_c_sp = np.zeros((N_days_spring,11*24,years))       
        train_test_su = np.zeros((N_days_summer,len(keys)*24,years))  
        train_test_c_su = np.zeros((N_days_summer,11*24,years))       
        train_test_a = np.zeros((N_days_autumn,len(keys)*24,years))  
        train_test_c_a = np.zeros((N_days_autumn,11*24,years))   
        
        # seasonal training/testing datasets with all years concatenated (2D)
        train_winter = np.zeros((N_days_winter*years_train,len(keys)*24))
        train_winter_c = np.zeros((N_days_winter*years_train,11*24))
        train_spring = np.zeros((N_days_spring*years_train,len(keys)*24))
        train_spring_c = np.zeros((N_days_spring*years_train,11*24))
        train_summer = np.zeros((N_days_summer*years_train,len(keys)*24))
        train_summer_c = np.zeros((N_days_summer*years_train,11*24))
        train_autumn = np.zeros((N_days_autumn*years_train,len(keys)*24))
        train_autumn_c = np.zeros((N_days_autumn*years_train,11*24))
        
        test_winter = np.zeros((N_days_winter*years_test,len(keys)*24))
        test_winter_c = np.zeros((N_days_winter*years_test,11*24))
        test_spring = np.zeros((N_days_spring*years_test,len(keys)*24))
        test_spring_c = np.zeros((N_days_spring*years_test,11*24))
        test_summer = np.zeros((N_days_summer*years_test,len(keys)*24))
        test_summer_c = np.zeros((N_days_summer*years_test,11*24))
        test_autumn = np.zeros((N_days_autumn*years_test,len(keys)*24))
        test_autumn_c = np.zeros((N_days_autumn*years_test,11*24))
        
        # Input
        DF_w=[]
        DF_sp=[]
        DF_su=[]
        DF_a=[]
        for year in range(year_start,year_end+1):
        
            DF_w=DF_winter[year-year_start]
            DF_sp=DF_spring[year-year_start]
            DF_su=DF_summer[year-year_start]
            DF_a=DF_autumn[year-year_start]
            
            # Seasonal datasets with only uncertain parameters and with all parameters
            train_test_year_w = []
            train_test_c_year_w= []
            train_test_year_sp = []
            train_test_c_year_sp= []
            train_test_year_su = []
            train_test_c_year_su= []
            train_test_year_a = []
            train_test_c_year_a = []
            
            x=0
            while x in range(len(DF_w)):
                    
                    train_test_df = np.concatenate([DF_w[key].values[x:x+24] for key in keys])
                    train_test_year_w.append(train_test_df)
                    
                    train_test_df_c = np.concatenate([DF_w['P'].values[x:x+24],DF_w['T2m'].values[x:x+24],DF_w['PZ'].values[x:x+24],
                                                      DF_w['Res_el'].values[x:x+24],DF_w['Ter_el'].values[x:x+24],
                                                      DF_w['Com_el'].values[x:x+24],DF_w['Pub_el'].values[x:x+24],
                                                      DF_w['Res_th'].values[x:x+24],DF_w['Ter_th'].values[x:x+24],
                                                      DF_w['Com_th'].values[x:x+24],DF_w['Pub_th'].values[x:x+24]])  
                    train_test_c_year_w.append(train_test_df_c)
                    
                    x=x+24
            
            train_test_w[:,:,year-year_start]=train_test_year_w
            train_test_c_w[:,:,year-year_start]=train_test_c_year_w
            
            x=0
            while x in range(len(DF_sp)):
                    
                    train_test_df = np.concatenate([DF_sp[key].values[x:x+24] for key in keys])
                    train_test_year_sp.append(train_test_df)
                    
                    train_test_df_c = np.concatenate([DF_sp['P'].values[x:x+24],DF_sp['T2m'].values[x:x+24],DF_sp['PZ'].values[x:x+24],
                                                      DF_sp['Res_el'].values[x:x+24],DF_sp['Ter_el'].values[x:x+24],
                                                      DF_sp['Com_el'].values[x:x+24],DF_sp['Pub_el'].values[x:x+24],
                                                      DF_sp['Res_th'].values[x:x+24],DF_sp['Ter_th'].values[x:x+24],
                                                      DF_sp['Com_th'].values[x:x+24],DF_sp['Pub_th'].values[x:x+24]])  
                    train_test_c_year_sp.append(train_test_df_c)
                    
                    x=x+24
            
            train_test_sp[:,:,year-year_start]=train_test_year_sp
            train_test_c_sp[:,:,year-year_start]=train_test_c_year_sp
            
            x=0
            while x in range(len(DF_su)):
                    
                    train_test_df = np.concatenate([DF_su[key].values[x:x+24] for key in keys])
                    train_test_year_su.append(train_test_df)
                    
                    train_test_df_c = np.concatenate([DF_su['P'].values[x:x+24],DF_su['T2m'].values[x:x+24],DF_su['PZ'].values[x:x+24],
                                                      DF_su['Res_el'].values[x:x+24],DF_su['Ter_el'].values[x:x+24],
                                                      DF_su['Com_el'].values[x:x+24],DF_su['Pub_el'].values[x:x+24],
                                                      DF_su['Res_th'].values[x:x+24],DF_su['Ter_th'].values[x:x+24],
                                                      DF_su['Com_th'].values[x:x+24],DF_su['Pub_th'].values[x:x+24]])  
                    train_test_c_year_su.append(train_test_df_c)
                    
                    x=x+24
            
            train_test_su[:,:,year-year_start]=train_test_year_su
            train_test_c_su[:,:,year-year_start]=train_test_c_year_su
            
            x=0
            while x in range(len(DF_a)):
                    
                    train_test_df = np.concatenate([DF_a[key].values[x:x+24] for key in keys])
                    train_test_year_a.append(train_test_df)
                    
                    train_test_df_c = np.concatenate([DF_a['P'].values[x:x+24],DF_a['T2m'].values[x:x+24],DF_a['PZ'].values[x:x+24],
                                                      DF_a['Res_el'].values[x:x+24],DF_a['Ter_el'].values[x:x+24],
                                                      DF_a['Com_el'].values[x:x+24],DF_a['Pub_el'].values[x:x+24],
                                                      DF_a['Res_th'].values[x:x+24],DF_a['Ter_th'].values[x:x+24],
                                                      DF_a['Com_th'].values[x:x+24],DF_a['Pub_th'].values[x:x+24]])  
                    train_test_c_year_a.append(train_test_df_c)
                    
                    x=x+24
            
            train_test_a[:,:,year-year_start]=train_test_year_a
            train_test_c_a[:,:,year-year_start]=train_test_c_year_a
                
        train_test_w = np.array(train_test_w)
        train_test_c_w = np.array(train_test_c_w)
        train_test_sp = np.array(train_test_sp)
        train_test_c_sp = np.array(train_test_c_sp)
        train_test_su = np.array(train_test_su)
        train_test_c_su = np.array(train_test_c_su)
        train_test_a = np.array(train_test_a)
        train_test_c_a = np.array(train_test_c_a)
        
        # Concatenate the arrays over the different years for the multi-year clustering 
        for year in range(years_train):
            train_winter[year*N_days_winter:(year+1)*N_days_winter] = np.concatenate([train_test_w[:,:,year]])
            train_winter_c[year*N_days_winter:(year+1)*N_days_winter] = np.concatenate([train_test_c_w[:,:,year]])
            train_spring[year*N_days_spring:(year+1)*N_days_spring] = np.concatenate([train_test_sp[:,:,year]])
            train_spring_c[year*N_days_spring:(year+1)*N_days_spring] = np.concatenate([train_test_c_sp[:,:,year]])
            train_summer[year*N_days_summer:(year+1)*N_days_summer] = np.concatenate([train_test_su[:,:,year]])
            train_summer_c[year*N_days_summer:(year+1)*N_days_summer] = np.concatenate([train_test_c_su[:,:,year]])
            train_autumn[year*N_days_autumn:(year+1)*N_days_autumn] = np.concatenate([train_test_a[:,:,year]])
            train_autumn_c[year*N_days_autumn:(year+1)*N_days_autumn] = np.concatenate([train_test_c_a[:,:,year]])
        
        for year in range(years_test):
            test_winter[year*N_days_winter:(year+1)*N_days_winter] = np.concatenate([train_test_w[:,:,years_train+year]])
            test_winter_c[year*N_days_winter:(year+1)*N_days_winter] = np.concatenate([train_test_c_w[:,:,years_train+year]])
            test_spring[year*N_days_spring:(year+1)*N_days_spring] = np.concatenate([train_test_sp[:,:,years_train+year]])
            test_spring_c[year*N_days_spring:(year+1)*N_days_spring] = np.concatenate([train_test_c_sp[:,:,years_train+year]])
            test_summer[year*N_days_summer:(year+1)*N_days_summer] = np.concatenate([train_test_su[:,:,years_train+year]])
            test_summer_c[year*N_days_summer:(year+1)*N_days_summer] = np.concatenate([train_test_c_su[:,:,years_train+year]])
            test_autumn[year*N_days_autumn:(year+1)*N_days_autumn] = np.concatenate([train_test_a[:,:,years_train+year]])
            test_autumn_c[year*N_days_autumn:(year+1)*N_days_autumn] = np.concatenate([train_test_c_a[:,:,years_train+year]])
        
        # Number of columns for the complete annual dataset
        lenghth_features_c = len(train_test_c_w[1])
        
    ######################## Extreme days of zonal price ######################
    
    ### Extreme days are identified in the training dataset (2005-2017)
    
    if seasons == 1:
        print("Choose the season for clustering")
        season = str(input())
    
    if seasons==0:   ### WITHOUT DIVISION BY SEASONS
        train_analysed = train
    elif seasons==1: ### WITH DIVISION BY SEASONS: extreme days for each season
          if season == "winter":
            train_analysed = train_winter
          elif season == "spring":
              train_analysed = train_spring
          elif season == "summer":
              train_analysed = train_summer
          elif season == "autumn":
              train_analysed = train_autumn     
    
    ### PRICES
    ### Identify the two extreme days with the lowest and highest hourly prices
    
    train_analyzed_PUN = train_analysed[:,48:72]
    
    ### Hourly min/max prices
    ### In case the extreme max price is higher than 1000 €/MWh
    ### Take the second highest price (as it occurs for Palermo)
    
    # Min price
    index_extreme_min_dh = np.argmin(train_analyzed_PUN)
    index_extreme_min_day, index_extreme_min_hour = divmod(index_extreme_min_dh, train_analyzed_PUN.shape[1])
    
    # Max price
    # Extract all hourly values for the zonal price time window
    all_hourly_prices = train_analyzed_PUN.flatten()
    # Sort values to get top 2
    sorted_indices = np.argsort(all_hourly_prices)  # ascending order
    # Maximum: if value is unusually high (e.g. > 1000 €/MWh), take second highest
    price_threshold = 1000
    if all_hourly_prices[sorted_indices[-1]] > price_threshold:
        # print(f"Extreme max price {all_hourly_prices[sorted_indices[-1]]:.2f} €/MWh exceeds threshold. Using second highest.")
        index_extreme_max_dh = sorted_indices[-2]
    else:
        index_extreme_max_dh = sorted_indices[-1]
    index_extreme_max_day, index_extreme_max_hour = divmod(index_extreme_max_dh, train_analyzed_PUN.shape[1])
    
    ### We take the days with minimum/maximum hourly values
    index_extreme_min = index_extreme_min_day
    index_extreme_max = index_extreme_max_day
    
    ######################### Plot timeseries in the whole dataset ############
    ############################# Without seasons #############################
        
    DF_train_test.set_index('Date', inplace=True)
    
    # Plot of zonal market prices to identify peaks in the whole dataset
    plots_fold = os.path.join(cfold, 'Plots')   # folder to save plots
    PATH_plots = plots_fold+"/Timeseries/"+loc
    
    # Prices of extreme days in the training dataset (excluding unusual prices higher than 1000 €/MWh)
    # Get timestamps and corresponding price values
    time_min = DF_train_test.index[index_extreme_min_dh]
    price_min = DF_train_test['PZ'].iloc[index_extreme_min_dh]
    time_max = DF_train_test.index[index_extreme_max_dh]
    price_max = DF_train_test['PZ'].iloc[index_extreme_max_dh]
    
    # Detect 1% top price peaks in the whole dataset (2005-2023)
    peaks, _ = find_peaks(DF_train_test['PZ'], height=DF_train_test['PZ'].quantile(0.99))  # Top 1% prices as peaks
    peak_prices = DF_train_test['PZ'].iloc[peaks]
    # Highest peak price
    index_top_peak = peaks[np.argmax(peak_prices)]
    time_top_peak = DF_train_test.index[index_top_peak]
    price_top_peak = DF_train_test['PZ'].iloc[index_top_peak]    
    
    # Initialize second peak price variables (set to None by default)
    time_second_peak = None
    price_second_peak = None
    second_peak_annotated = False
    
    # Check if top price exceeds threshold
    if price_top_peak > price_threshold:  
        sorted_peak_values = peak_prices.sort_values(ascending=False)
        time_second_peak = sorted_peak_values.index[1]  # Timestamp of 2nd highest peak
        price_second_peak = sorted_peak_values.iloc[1]  # Its value
        second_peak_annotated = True
    
    if plot_ts == "yes":
        
        ######################### Plot prices in the whole dataset ################
    
        # Plot full time series, 1% top price peaks and prices of extreme days
        plt.figure(figsize=(24,12))
        a1=plt.plot(DF_train_test.index, DF_train_test['PZ'], label='Hourly price', alpha=0.6)
        a2=plt.scatter(DF_train_test.index[peaks], peak_prices, color='red', label='Top 1% price peaks')
        a3=plt.scatter(time_min, price_min, color='blue', s=200, marker='v', label="Min hourly price training dataset, "+f"{year_start}-"+str(year_start+years_train-1))
        a4=plt.scatter(time_max, price_max, color='orange', s=200, marker='^', label="Max hourly price training dataset, "+f"{year_start}-"+str(year_start+years_train-1))
        a5=plt.scatter(time_top_peak, price_top_peak, color='green', s=200, marker='*', label='Highest peak (Top 1%)')
        
        if second_peak_annotated:
            a6 = plt.scatter(time_second_peak, price_second_peak, color='purple', s=200, marker='*',
                             label='Second highest peak (Top 1%)')
        
            plt.annotate(f"{time_second_peak.strftime('%Y-%m-%d %H:%M')}\n€/MWh {price_second_peak:.2f}",
                          (time_second_peak, price_second_peak),
                          textcoords="offset points", xytext=(70, -10), ha='center', fontsize=16, color='purple')
        
        
        # Add annotations
        plt.annotate(f"{time_min.strftime('%Y-%m-%d %H:%M')}\n€/MWh {price_min:.2f}",
                      (time_min, price_min),
                      textcoords="offset points", xytext=(80, -28), ha='center', fontsize=16, color='blue')
        
        plt.annotate(f"{time_max.strftime('%Y-%m-%d %H:%M')}\n€/MWh {price_max:.2f}",
                      (time_max, price_max),
                      textcoords="offset points", xytext=(0, 10), ha='center', fontsize=16, color='orange')
        
        plt.annotate(f"{time_top_peak.strftime('%Y-%m-%d %H:%M')}\n€/MWh {price_top_peak:.2f}",
                      (time_top_peak, price_top_peak),
                      textcoords="offset points", xytext=(75, -10), ha='center', fontsize=16, color='green')
        
        # Formatting x-axis to show all years
        plt.gca().xaxis.set_major_locator(mdates.YearLocator(1))  # One tick per year
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y'))  # Format as year
        plt.gca().xaxis.set_tick_params(labelsize = 16)
        plt.gca().yaxis.set_tick_params(labelsize = 16)
        
        if loc=='Padova' or loc=='Trieste':
           plt.title(f"Zonal market prices, {year_start}-{year_end}, {loc} in northern Italy", fontsize=24)
        elif loc=='Palermo':
             plt.title(f"Zonal market prices, {year_start}-{year_end}, {loc} in southern Italy", fontsize=24)
        ### Add here another "if" statement for other locations

        plt.xlabel("Year", fontsize=24)
        plt.ylabel("Zonal market prices [€/MWh]", fontsize=24)
        legend_handles = [a1[0], a2, a3, a4, a5]
        if second_peak_annotated:
            legend_handles.append(a6)
        plt.legend(handles=legend_handles, loc='upper left', prop=dict(size=20), bbox_to_anchor=(0.005,0.99))  
        plt.savefig(os.path.join(PATH_plots, "zonal_prices_"+f"{year_start}-{year_end}.jpg"), bbox_inches='tight')
        
        ############################# Plot ambient temperature in the whole dataset #############################
        plt.figure(figsize=(24,12))
        a1=plt.plot(DF_train_test.index, DF_train_test['T2m'], label='Ambient temperature', alpha=0.6)
        
        # Formatting x-axis to show all years
        plt.gca().xaxis.set_major_locator(mdates.YearLocator(1))  # One tick per year
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y'))  # Format as year
        plt.gca().xaxis.set_tick_params(labelsize = 16)
        plt.gca().yaxis.set_tick_params(labelsize = 16)
        
        if loc=='Padova' or loc=='Trieste':
           plt.title(f"Ambient temperature, {year_start}-{year_end}, {loc} in northern Italy", fontsize=24)
        elif loc=='Palermo':
             plt.title(f"Ambient temperature, {year_start}-{year_end}, {loc} in southern Italy", fontsize=24)
        plt.xlabel("Year", fontsize=24)
        plt.ylabel("Ambient temperature [°C]", fontsize=24)
        plt.savefig(os.path.join(PATH_plots, "ambient_temperature_"+f"{year_start}-{year_end}.jpg"), bbox_inches='tight')
        
        ############################# Plot PV capacity factor in the whole dataset #############################
        plt.figure(figsize=(24,12))
        a1=plt.plot(DF_train_test.index, DF_train_test['P']/1000, label='PV capacity', alpha=0.6)
        
        # Formatting x-axis to show all years
        plt.gca().xaxis.set_major_locator(mdates.YearLocator(1))  # One tick per year
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y'))  # Format as year
        plt.gca().xaxis.set_tick_params(labelsize = 16)
        plt.gca().yaxis.set_tick_params(labelsize = 16)
        
        if loc=='Padova' or loc=='Trieste':
           plt.title(f"Solar availability, {year_start}-{year_end}, {loc} in northern Italy", fontsize=24)
        elif loc=='Palermo':
             plt.title(f"Solar availability, {year_start}-{year_end}, {loc} in southern Italy", fontsize=24)
        plt.xlabel("Year", fontsize=24)
        plt.ylabel("PV capacity factor [-]", fontsize=24)
        plt.savefig(os.path.join(PATH_plots, "solar_"+f"{year_start}-{year_end}.jpg"), bbox_inches='tight')
    
    ####################### CLUSTERING THE TRAINING DATASET ###################
    
    ### MULTI-YEAR CLUSTERING ###
    ### Clustering is applied to the whole training dataset (years 2005--->2017)
    ### Without or with division by seasons
    
    # Number of rows: number of daily timeseries in the training dataset 
    # Number of columns: 24*3 (timesteps in a day)
    TS = train_analysed 
    nts = TS.shape[0] # number of rows of TS (number of observations)
    sz = TS.shape[1]  # number of columns of TS (daily segmentation for each attribute)
    
    # Normalization 
    # X_std = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0)) ---> minimum and maximum calculated for each column
    # X_scaled = X_std * (max - min) + min
    # min-max is the feature range ---> default is 0-1
    TS_I = MinMaxScaler().fit_transform(TS)   # normalized array
    
    clust_fold = os.path.join(cfold, 'Clustering')   # folder to save clustering results
    PATH_clustering = clust_fold+"/"+loc
    
    if seasons==0:    ### WITHOUT DIVISION BY SEASONS
       n_cases = 49   # 49 cases with different numbers of clusters (from 2 to 50 with replacing, from 4 to 52 with adding)
    elif seasons==1:  ### WITH DIVISION BY SEASONS
        n_cases = 5   # from 2 to 6 clusters for each season with replacing (from 8 to 24 annually), from 4 to 8 clusters for each season with adding (from 16 to 32 annually)
    
    index=np.zeros((nts,n_cases)) # cluster indexes assigned to all ts for each case (number of clusters)
    
    ### Metrics of clustering quality for each case
    ASV = np.zeros(n_cases)  # Average Silhouette Value (ASV) (measuring inertia, the sum of the distances between the ts and their closest centers)
    SSE = np.zeros(n_cases)  # Sum of Squared Errors (SSE) for elbow method
    DBI = np.zeros(n_cases)  # Davies-Bouldin Index (DBI), measuring separation and compactness, lower is better 
    Dunn = np.zeros(n_cases) # Dunn index, higher is better
    def dunn_index(X, labels): 
        clusters = np.unique(labels)
        intra_dists = []
        inter_dists = []
        
        for i in clusters:
            cluster_i = X[labels == i]
            intra_dists.append(np.max(cdist(cluster_i, cluster_i)))
            for j in clusters:
                if i < j:
                    cluster_j = X[labels == j]
                    inter_dists.append(np.min(cdist(cluster_i, cluster_j)))
        
        return min(inter_dists) / max(intra_dists) if intra_dists else 0
    
    centers=np.empty((n_cases+1,sz*n_cases)) # centroid for each cluster in each case
    centers[:]=np.nan
    index_representative_r=np.zeros((n_cases,n_cases+1)) # indexes (rows) of representative ts (rts) in each cluster for each number of clusters
    
    number_elements_each_cluster=np.zeros((n_cases,n_cases+1)) # number of ts of each cluster (column) for each number of clusters (row)
    d_ts=np.empty((nts,sz*(n_cases+1),n_cases))                # denormalized ts in each cluster (matrix 2-D) for each number of clusters (matrix 3-D)
    d_ts[:]=np.nan
    d_rts=np.zeros((n_cases+1,sz,n_cases))                     # denormalized rts in each cluster (matrix 2-D) for each number of clusters (matrix 3-D)
    p=np.zeros((n_cases,n_cases+1))                            # cluster frequency of each cluster (column) for each number of clusters (row)
    
    start_time = time.time()
    
    for k in range(2,n_cases+2):  # different numbers of clusters
        
        # clustering 
        clustering = KMeans(n_clusters=k, n_init=50).fit(TS_I) #clustering options
        index[:,k-2]=clustering.labels_
        ASV[k-2] = metrics.silhouette_score(TS_I, index[:,k-2], metric='euclidean') #ASV does not depend on the algorithm but only on the clusters found
        SSE[k-2] = clustering.inertia_
        DBI[k-2] = metrics.davies_bouldin_score(TS_I, index[:,k-2])
        Dunn[k-2] = dunn_index(TS_I, index[:,k-2])
        
        # centroids
        for yi in range(k):
            centers[yi,sz*(k-2):sz*(k-2)+sz]=clustering.cluster_centers_[yi].ravel() #barycenters calculated by arithmetic mean with "metric"=euclidean or by DBA with "metric"=dtw
        
        # representative timeseries
        distance_to_centers=clustering.transform(TS_I) #array with the distances of ts respect to the cluster centers
        min_distance_to_centers=distance_to_centers.min(axis=0) #minimum distances from the cluster centers
        for yi in range(k):
            rows = np.where(distance_to_centers[:,yi] == min_distance_to_centers[yi]) #index of the rts in each cluster, the one with minimum distance respect to its center
            index_representative_r[k-2,yi] = rows[0][0]  
            index_representative_r=index_representative_r.astype(int)
        
        # Denormalized ts and rts
        for yi in range(k):         # different clusters
            xx=np.where(index[:,k-2] == yi) #indexes of ts in each cluster yi for a number of clusters equal to "k"
            xx=np.array(xx)
            xx=xx.transpose()
            number_elements_each_cluster[k-2,yi]=len(xx)
            for xxx in range(len(xx)):
                d_ts[xxx,sz*yi:sz*yi+sz,k-2]=TS[xx[xxx]]       # denormalization of ts
            d_rts[yi,:,k-2]=TS[index_representative_r[k-2,yi]] # denormalization of rts
            
            # CLuster frequencies (weights of the typical days)
            p[k-2,yi]=number_elements_each_cluster[k-2,yi]/nts 
    
    ### TIME TO OBTAIN THE CLUSTERS VARYING THE NUMBER OF TYPICAL DAYS ###
    print("--- %s seconds ---" % (time.time() - start_time))
    
    # Save results for analysis
    np.save(os.path.join(PATH_clustering, 'ASV.npy'), ASV)
    np.save(os.path.join(PATH_clustering, 'SSE.npy'), SSE)
    np.save(os.path.join(PATH_clustering, 'DBI.npy'), DBI)
    np.save(os.path.join(PATH_clustering, 'Dunn.npy'), Dunn)
    
    ####################### Plot of ASV, SSE, DBI and Dunn ####################
    
    number_clusters = np.arange(2,n_cases+2,1)
    
    # Load the saved metrics
    ASV = np.load(os.path.join(PATH_clustering, 'ASV.npy'))
    SSE = np.load(os.path.join(PATH_clustering, 'SSE.npy'))
    DBI = np.load(os.path.join(PATH_clustering, 'DBI.npy'))
    Dunn = np.load(os.path.join(PATH_clustering, 'Dunn.npy'))
    
    fontsz = 20
    labelsz = 18
    
    # Create subplots for comparison
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))
    axs = axs.ravel()
    fig.subplots_adjust(left=0.06, right=0.98, top=0.90, bottom=0.06)
    
    plt.suptitle(f"Location: {loc}", fontsize=24, fontweight='bold')
    
    # ASV Plot
    axs[0].plot(number_clusters, ASV, marker='o', linewidth=2, color='blue')
    axs[0].set_xlabel("Number of clusters", fontsize=fontsz)
    axs[0].xaxis.set_tick_params(labelsize = labelsz)
    axs[0].set_ylabel("ASV", fontsize=fontsz)
    axs[0].yaxis.set_tick_params(labelsize = labelsz)
    axs[0].set_title("Average Silhouette Value (ASV)", fontsize=fontsz)
    axs[0].grid(True)
    
    # SSE Plot
    axs[1].plot(number_clusters, SSE, marker='o', linewidth=2, color='red')
    axs[1].set_xlabel("Number of clusters", fontsize=fontsz)
    axs[1].xaxis.set_tick_params(labelsize = labelsz)
    axs[1].set_ylabel("SSE", fontsize=fontsz)
    axs[1].yaxis.set_tick_params(labelsize = labelsz)
    axs[1].set_title("Sum of Squared Errors (SSE)", fontsize=fontsz)
    axs[1].grid(True)
    
    # DBI Plot
    axs[2].plot(number_clusters, DBI, marker='o', linewidth=2, color='green')
    axs[2].set_xlabel("Number of clusters", fontsize=fontsz)
    axs[2].xaxis.set_tick_params(labelsize = labelsz)
    axs[2].set_ylabel("DBI", fontsize=fontsz)
    axs[2].yaxis.set_tick_params(labelsize = labelsz)
    axs[2].set_title("Davies-Bouldin Index (DBI)", fontsize=fontsz)
    axs[2].grid(True)
    
    # Dunn Index Plot
    axs[3].plot(number_clusters, Dunn, marker='o', linewidth=2, color='purple')
    axs[3].set_xlabel("Number of clusters", fontsize=fontsz)
    axs[3].xaxis.set_tick_params(labelsize = labelsz)
    axs[3].set_ylabel("DI", fontsize=fontsz)
    axs[3].yaxis.set_tick_params(labelsize = labelsz)
    axs[3].set_title("Dunn Index (DI)", fontsize=fontsz)
    axs[3].grid(True)
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(os.path.join(PATH_clustering, 'Clustering_Metrics_Comparison.jpg'), dpi=1000)
    
    ################## EXTREME DAYS OF ZONAL ELECTRICITY PRICES ###############
    
    # print('Choose the method of considering the extreme days')
    # extreme_method=input()
    
    if extreme_method == 'R':  # Replace
    
        ### Method: "Replace representative period", Kotzur (2018)
        ### It could happen that the two extreme scenarios are
        ### in the same cluster: only the max extreme scenario becomes
        ### the representative of the cluster
        
        # Find the clusters of the extreme days for each different case
        cluster_extreme_min = np.zeros(n_cases)
        cluster_extreme_max = np.zeros(n_cases)
        for k in range(2,n_cases+2):
            cluster_extreme_min[k-2] = index[index_extreme_min,k-2]
            cluster_extreme_max[k-2] = index[index_extreme_max,k-2]
        
        # Replace the representative days of the clusters including the extreme days
        # with the extreme days
        index_representative_r_new = np.zeros((n_cases,n_cases+1))
        for k in range(2,n_cases+2):
            index_representative_r_new[k-2,:] = index_representative_r[k-2,:]
            index_representative_r_new[k-2,int(cluster_extreme_min[k-2])] = index_extreme_min
            index_representative_r_new[k-2,int(cluster_extreme_max[k-2])] = index_extreme_max
    
    
    elif extreme_method == 'A': # Adding
        
         ### Method: "Adding cluster center", Kotzur (2018)
         
         # Different cases  
         index_representative_r_new = np.zeros((n_cases,n_cases+3)) # new indexes of typical days
         centers_new = np.empty((n_cases+3,sz*n_cases)) # new centroids
         centers_new[:] = np.nan
         dist_new = np.zeros((nts,n_cases+3,n_cases)) # new distances of all elements to all centroids
         min_dist_new = np.empty((nts,n_cases))       # minimum distance of each element to a specific centroid
         min_dist_new[:] = np.nan
         index_new = np.zeros((nts,n_cases))          # new cluster assignment
         number_elements_each_cluster_new = np.zeros((n_cases,n_cases+3)) # new cluster weights
         
         for k in range(2,n_cases+2):
             
             # Update indexes of typical days
             index_representative_r_new[k-2,:] = np.append(index_representative_r[k-2,:],[int(index_extreme_min),int(index_extreme_max)], axis=0)
             # Get the last two elements of the row
             last_two_index = index_representative_r_new[k-2, -2:]
             # Insert the last two elements at positions k and k+1
             index_representative_r_new[k-2, k:k+2] = last_two_index
             
             # Append new centers (the two extreme scenarios become centroids)
             centers_new[0:k,sz*(k-2):sz*(k-2)+sz] = centers[0:k,sz*(k-2):sz*(k-2)+sz]
             centers_new[k, sz*(k-2):sz*(k-2)+sz] = TS_I[int(index_extreme_min)]
             centers_new[k+1, sz*(k-2):sz*(k-2)+sz] = TS_I[int(index_extreme_max)]
        
         # Remove duplicates
         index_representative_r_new[0:-1, -1] = 0 
         index_representative_r_new[0:2, -2] = 0
         
         # New cluster assignment
         
         for k in range(2,n_cases+2):
             for yi in range(k+2):
                 for i in range(nts):
                     
                     # Distance of each element "i" to each centroid "yi" for each case "k-2"
                     dist_new[i,yi,k-2] = np.sum(abs(TS_I[i]-centers_new[yi,sz*(k-2):sz*(k-2)+sz]))
                     
                     # Minimum distance of each element to a specific centroid/cluster for each case
                     min_dist_new[i,k-2] = np.min(dist_new[i,0:k+2,k-2])
                                      
                     # Find the index of the minimum distance for each element --> cluster assignment
                     index_new[i,k-2] = np.argmin(dist_new[i,0:k+2,k-2])
                   
         # New weights of clusters
         # The last two weights for each case are the weights of the extreme minimum and maximum scenarios, respectively
         
         for k in range(2,n_cases+2):
             for yi in range(k+2):
                 number_elements_each_cluster_new[k-2,yi] = np.sum(index_new[:,k-2] == yi)
                 
    ######################## Indexes of typical days ##########################
    ### "index_representative_r_new": indexes from 0 to 4745 (13 years)
    ### "index_representative_r_tot": indexes from 0 to 365 (1 year)
    
    if extreme_method == 'R':  # Replace
       index_representative_r_tot = np.zeros((n_cases,n_cases+1))
    elif extreme_method == 'A': # Adding
         index_representative_r_tot = np.zeros((n_cases,n_cases+3))
    
    for k in range(2,n_cases+2):
        
        if extreme_method == 'R':  # Replace
            for yi in range(k):
                
                if seasons==0: ### WITHOUT DIVISION BY SEASONS
                    
                   index_representative_r_tot[k - 2, yi] = index_representative_r_new[k - 2, yi] % 365
                
                elif seasons==1: ### WITH DIVISION BY SEASONS: extreme days for each season
                     
                     season_length = len(TS) / years_train  # Compute season length once
                     index_representative_r_tot[k - 2, yi] = index_representative_r_new[k - 2, yi] % season_length
        
        elif extreme_method == 'A':  # Adding
                for yi in range(k+2):
                    
                    if seasons==0: ### WITHOUT DIVISION BY SEASONS
                    
                       index_representative_r_tot[k - 2, yi] = index_representative_r_new[k - 2, yi] % 365
                    
                    elif seasons==1: ### WITH DIVISION BY SEASONS: extreme days for each season
                    
                         season_length = len(TS) / years_train  # Compute season length once
                         index_representative_r_tot[k - 2, yi] = index_representative_r_new[k - 2, yi] % season_length
                         
    ################################ Save results ##############################
    
    if extreme_method == 'R':  # Replace
    
       if seasons==0: ### WITHOUT DIVISION BY SEASONS
        
           np.savez(os.path.join(PATH_clustering,'Clustering_train_'+extreme_method),
                    lenghth_features_c = lenghth_features_c, TS_I = TS_I, nts = nts, sz = sz, n_cases = n_cases, index = index,
                    train = train, train_c = train_c, test = test, test_c = test_c,
                    d_ts = d_ts, d_rts = d_rts, index_representative_r = index_representative_r, p = p, number_elements_each_cluster = number_elements_each_cluster,
                    index_extreme_min = index_extreme_min, index_extreme_max = index_extreme_max, 
                    cluster_extreme_min = cluster_extreme_min, cluster_extreme_max = cluster_extreme_max, 
                    index_representative_r_new = index_representative_r_new, index_representative_r_tot = index_representative_r_tot,
                    train_c_years = train_c_years, test_c_years = test_c_years,
                    TS = TS, years_train = years_train, years_test = years_test, years = years, DF_no_s_mean = DF_no_s_mean,
                    extreme_method = extreme_method, seasons = seasons,
                    time_min = time_min, time_max = time_max, time_top_peak = time_top_peak, time_second_peak = time_second_peak,
                    price_min = price_min, price_max = price_max, price_top_peak = price_top_peak, price_second_peak = price_second_peak,
                    number_clusters = number_clusters, ASV = ASV, SSE = SSE, DBI = DBI, Dunn = Dunn)
        
       elif seasons==1: ### WITH DIVISION BY SEASONS: extreme days for each season
            
            np.savez(os.path.join(PATH_clustering,'Clustering_train_'+extreme_method+'_'+season),
                     lenghth_features_c = lenghth_features_c, TS_I = TS_I, nts = nts, sz = sz, n_cases = n_cases, index = index,
                     TS = TS, years_train = years_train, years_test = years_test, years = years, DF_no_s_mean = DF_no_s_mean,
                     DF_winter = DF_winter, DF_spring = DF_spring, DF_summer = DF_summer, DF_autumn = DF_autumn,
                     DF_winter_mean = DF_winter_mean, DF_spring_mean = DF_spring_mean, DF_summer_mean = DF_summer_mean, DF_autumn_mean = DF_autumn_mean,
                     N_days_winter = N_days_winter, N_days_spring = N_days_spring, N_days_summer = N_days_summer, N_days_autumn = N_days_autumn,
                     train_test_w = train_test_w, train_test_c_w = train_test_c_w,
                     train_test_sp = train_test_sp, train_test_c_sp = train_test_c_sp,
                     train_test_su = train_test_su, train_test_c_su = train_test_c_su,
                     train_test_a = train_test_a, train_test_c_a = train_test_c_a,
                     train_winter = train_winter, train_winter_c = train_winter_c,
                     train_spring = train_spring, train_spring_c = train_spring_c,
                     train_summer = train_summer, train_summer_c = train_summer_c,
                     train_autumn = train_autumn, train_autumn_c = train_autumn_c,
                     test_winter = test_winter, test_winter_c = test_winter_c,
                     test_spring = test_spring, test_spring_c = test_spring_c,
                     test_summer = test_summer, test_summer_c = test_summer_c,
                     test_autumn = test_autumn, test_autumn_c = test_autumn_c,
                     d_ts = d_ts, d_rts = d_rts, index_representative_r = index_representative_r, p = p, number_elements_each_cluster = number_elements_each_cluster,
                     index_extreme_min = index_extreme_min, index_extreme_max = index_extreme_max, 
                     cluster_extreme_min = cluster_extreme_min, cluster_extreme_max = cluster_extreme_max, 
                     index_representative_r_new = index_representative_r_new, index_representative_r_tot = index_representative_r_tot,
                     extreme_method = extreme_method, seasons = seasons,
                     time_min = time_min, time_max = time_max, time_top_peak = time_top_peak, time_second_peak = time_second_peak,
                     price_min = price_min, price_max = price_max, price_top_peak = price_top_peak, price_second_peak = price_second_peak,
                     number_clusters = number_clusters, ASV = ASV, SSE = SSE, DBI = DBI, Dunn = Dunn)
            
    elif extreme_method == 'A': # Additional
    
         if seasons==0: ### WITHOUT DIVISION BY SEASONS
        
             np.savez(os.path.join(PATH_clustering,'Clustering_train_'+extreme_method),
                      lenghth_features_c = lenghth_features_c, TS_I = TS_I, nts = nts, sz = sz, n_cases = n_cases, index_new = index_new,
                      train = train, train_c = train_c, test = test, test_c = test_c,
                      number_elements_each_cluster_new = number_elements_each_cluster_new, 
                      index_representative_r_new = index_representative_r_new, index_representative_r_tot = index_representative_r_tot,
                      train_c_years = train_c_years, test_c_years = test_c_years,
                      TS = TS, years_train = years_train, years_test = years_test, years = years, DF_no_s_mean = DF_no_s_mean,
                      extreme_method = extreme_method, seasons = seasons,
                      index_extreme_min = index_extreme_min, index_extreme_max = index_extreme_max,
                      time_min = time_min, time_max = time_max, time_top_peak = time_top_peak, time_second_peak = time_second_peak,
                      price_min = price_min, price_max = price_max, price_top_peak = price_top_peak, price_second_peak = price_second_peak,
                      number_clusters = number_clusters, ASV = ASV, SSE = SSE, DBI = DBI, Dunn = Dunn)
             
         elif seasons==1: ### WITH DIVISION BY SEASONS: extreme days for each season
         
              np.savez(os.path.join(PATH_clustering,'Clustering_train_'+extreme_method+'_'+season),
                       lenghth_features_c = lenghth_features_c, TS_I = TS_I, nts = nts, sz = sz, n_cases = n_cases, index_new = index_new,
                       TS = TS, years_train = years_train, years_test = years_test, years = years, DF_no_s_mean = DF_no_s_mean,
                       DF_winter = DF_winter, DF_spring = DF_spring, DF_summer = DF_summer, DF_autumn = DF_autumn,
                       DF_winter_mean = DF_winter_mean, DF_spring_mean = DF_spring_mean, DF_summer_mean = DF_summer_mean, DF_autumn_mean = DF_autumn_mean,
                       N_days_winter = N_days_winter, N_days_spring = N_days_spring, N_days_summer = N_days_summer, N_days_autumn = N_days_autumn,
                       train_test_w = train_test_w, train_test_c_w = train_test_c_w,
                       train_test_sp = train_test_sp, train_test_c_sp = train_test_c_sp,
                       train_test_su = train_test_su, train_test_c_su = train_test_c_su,
                       train_test_a = train_test_a, train_test_c_a = train_test_c_a,
                       train_winter = train_winter, train_winter_c = train_winter_c,
                       train_spring = train_spring, train_spring_c = train_spring_c,
                       train_summer = train_summer, train_summer_c = train_summer_c,
                       train_autumn = train_autumn, train_autumn_c = train_autumn_c,
                       test_winter = test_winter, test_winter_c = test_winter_c,
                       test_spring = test_spring, test_spring_c = test_spring_c,
                       test_summer = test_summer, test_summer_c = test_summer_c,
                       test_autumn = test_autumn, test_autumn_c = test_autumn_c,
                       number_elements_each_cluster_new = number_elements_each_cluster_new, 
                       index_representative_r_new = index_representative_r_new, index_representative_r_tot = index_representative_r_tot,
                       extreme_method = extreme_method, seasons = seasons,
                       index_extreme_min = index_extreme_min, index_extreme_max = index_extreme_max,
                       time_min = time_min, time_max = time_max, time_top_peak = time_top_peak, time_second_peak = time_second_peak,
                       price_min = price_min, price_max = price_max, price_top_peak = price_top_peak, price_second_peak = price_second_peak,
                       number_clusters = number_clusters, ASV = ASV, SSE = SSE, DBI = DBI, Dunn = Dunn)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# CHECK CLUSTERING

# Inputs:
# - some input are the same as in the previous function
# - opt_K: number of representative days for plots (type: int)

def post_clustering(loc="",
                    seasons=0,
                    extreme_method='A', 
                    opt_K=int()):
    
    ### Python packages
    import numpy as np                               # Numpy, for multi-dimensional array
    import math                                      # mathematical operations
    import os                                        # to select the folder for saving
    # Data visualization and plots
    import matplotlib.pyplot as plt                  # Pyplot, for plots, graphs, diagrams
    from matplotlib.lines import Line2D
    
    ### Import clustering results
    
    cfold = os.getcwd()  # working directory
    clust_fold = os.path.join(cfold, 'Clustering')   # folder to save clustering results
    PATH_clustering = clust_fold+"/"+loc
    
    if seasons==0: ### WITHOUT DIVISION BY SEASONS
       timeseries = np.load(os.path.join(PATH_clustering,'Clustering_train_'+extreme_method+'.npz'), allow_pickle=True)

       index_new = timeseries['index_new']        # cluster indeces
       index_representative_r_new = timeseries['index_representative_r_new'] # indeces of typical days
       number_elements_each_cluster_new = timeseries['number_elements_each_cluster_new'] # weights of typical days
       TS = timeseries['TS']     # training dataset with only uncertain parameters (clustering attributes) 
       train_c = timeseries['train_c']   # training dataset with all timeseries
       nts = timeseries['nts']   # number of timeseries in the training dataset
       sz = timeseries['sz']     # number of clustering attributs (3) X 24 hourly values
              
       # Timeseries in opt_K clusters and representative days
       d_ts=np.zeros((nts,sz*opt_K))  # all timeseries
       d_rts=np.zeros((opt_K,sz))     # representative days
       for yi in range(opt_K):         # different clusters 
           if extreme_method=='A':
              d_rts[yi,:]=TS[int(index_representative_r_new[int(opt_K-4),yi])] 
              xx=np.where(index_new[:,opt_K-4] == yi) 
           elif extreme_method=='R':
                d_rts[yi,:]=TS[int(index_representative_r_new[int(opt_K-2),yi])] 
                xx=np.where(index_new[:,opt_K-2] == yi)     
           xx=np.array(xx)
           xx=xx.transpose()
           for xxx in range(len(xx)):
               d_ts[xxx,sz*yi:sz*yi+sz]=TS[xx[xxx]]
       
       ### Plots
       
       hours = np.arange(1,25)
       
       ### Timeseries and representative days of prices divided by clusters
       
       # Layout logic
       if opt_K % 2 == 0 and opt_K <= 10:
          nrows = 2
          ncols = opt_K // 2
       elif opt_K % 2 == 0 and 12 <= opt_K <= 20:
            nrows = 4
            ncols = math.ceil(opt_K / 4)
       else:
            ncols = math.ceil(math.sqrt(opt_K))
            nrows = math.ceil(opt_K / ncols)
       
       fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 3*nrows), sharey=True)
       axes = axes.flatten()  # Flatten the 2D array of axes to index it with cluster number

       for yi in range(opt_K):
           ax = axes[yi] 

           # All timeseries for cluster yi — only zonal price (48:72)
           if extreme_method=='A':
              cluster_indices = np.where(index_new[:, opt_K-4] == yi)[0]
           elif extreme_method=='R':
                cluster_indices = np.where(index_new[:, opt_K-2] == yi)[0]
           for idx in cluster_indices:
               ax.plot(hours, TS[idx, 48:72], color='gray', alpha=0.3)

           # Representative day of zonal prices
           ax.plot(hours, d_rts[yi, 48:72], color='red', linewidth=2, label='Representative day')

           ax.set_title(f'Cluster {yi + 1}', fontsize=14)
           ax.grid(True)
           if yi % ncols == 0:
              ax.set_ylabel("Price [€/MWh]", fontsize=14)
           if yi >= ncols * (nrows - 1):
              ax.set_xlabel("Hour", fontsize=14)
           
           if loc=='Padova' or loc=='Trieste':
              ax.set_ylim(0,500) # for better visualization
           elif loc=='Palermo':  
                ax.set_ylim(0,1000) # for better visualization
           ax.xaxis.set_tick_params(labelsize = 14)
           ax.yaxis.set_tick_params(labelsize = 14)

       # Remove unused axes
       for j in range(opt_K, len(axes)):
           fig.delaxes(axes[j])

       plt.suptitle(f"{loc} - Timeseries of zonal market prices by cluster", fontsize=16)
       # Custom inline "legend" for the title
       fig.legend(
           handles=[Line2D([0], [0], color='red', lw=2, label='Representative days')],
           loc="upper center",
           bbox_to_anchor=(0.8, 0.98),  # slightly above title
           fontsize=14,
           frameon=False,
           ncol=1
       )
       plt.tight_layout(rect=[0, 0, 1, 0.97])
       plt.savefig(os.path.join(PATH_clustering, f"clusters_typical_prices_{opt_K}.jpg"), bbox_inches='tight')

       ### Representative days
       
       fontsz = 16
       labelsz = 14
       
       # PV capacity factor
       plt.figure(figsize=(10, 6))
       for i in range(opt_K):
           plt.plot(hours, d_rts[i, 0:24]/1000, label=f'Typical day {i+1}')
       plt.xlabel("Hour", fontsize=fontsz)
       plt.ylabel("PV capacity factor", fontsize=fontsz)
       plt.title(f"{loc} - Representative scenarios of PV capacity factor", fontsize=fontsz, fontweight='bold')
       # plt.legend(loc='upper right', fontsize=10)
       plt.xticks(fontsize=14)
       plt.yticks(fontsize=14)
       plt.grid(True)
       plt.tight_layout()
       plt.savefig(os.path.join(PATH_clustering, f"solar_typical_days_{opt_K}.jpg"), bbox_inches='tight')
        
       # Ambient temperature
       plt.figure(figsize=(10, 6))
       for i in range(opt_K):
           plt.plot(hours, d_rts[i, 24:48], label=f'Typical day {i+1}')
       plt.xlabel("Hour", fontsize=fontsz)
       plt.ylabel("Ambient temperature [°C]", fontsize=fontsz)
       plt.title(f"{loc} - Representative scenarios of ambient temperature", fontsize=fontsz, fontweight='bold')
       # plt.legend(loc='upper right', fontsize=10)
       plt.xticks(fontsize=14)
       plt.yticks(fontsize=14)
       plt.grid(True)
       plt.tight_layout()
       plt.savefig(os.path.join(PATH_clustering, f"temperature_typical_days_{opt_K}.jpg"), bbox_inches='tight')
      
       # Zonal market price 
       plt.figure(figsize=(10, 6))
       # representative days except extreme days
       for i in range(opt_K-2):
           plt.plot(hours, d_rts[i, 48:72], color='gray')
       # extreme days
       plt.plot(hours, d_rts[opt_K - 2, 48:72], color='red', linewidth=2.5, label='Extreme days')
       plt.plot(hours, d_rts[opt_K - 1, 48:72], color='red', linewidth=2.5)   
       plt.xlabel("Hour", fontsize=fontsz)
       plt.ylabel("Zonal market price [€/MWh]", fontsize=fontsz)
       plt.title(f"{loc} - Representative scenarios of zonal market price", fontsize=fontsz, fontweight='bold')
       plt.legend(loc='upper left', fontsize=fontsz)
       plt.xticks(fontsize=14)
       plt.yticks(fontsize=14)
       plt.grid(True)
       plt.tight_layout()
       plt.savefig(os.path.join(PATH_clustering, f"price_typical_days_{opt_K}.jpg"), bbox_inches='tight')
       
       ### Demands not considered as clustering attributes
       # Electricity demands
       fig2, ax2 = plt.subplots(figsize=(12, 10), nrows=2, ncols=2)
       fig2.tight_layout(pad = 3)
       fig2.subplots_adjust(left=0.06, right=0.98, top=0.90, bottom=0.06)
       plt.suptitle("Representative scenarios of electricity demand", fontsize=24, fontweight='bold')
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

       for yi in range(opt_K):
           ax2[0,0].plot(hours, train_c[yi,:][72:96], "b", linewidth=1)
           ax2[0,1].plot(hours, train_c[yi,:][96:120], "b", linewidth=1)
           ax2[1,0].plot(hours, train_c[yi,:][120:144], "b", linewidth=1)
           ax2[1,1].plot(hours, train_c[yi,:][144:168], "b", linewidth=1)

       ax2[0,0].set_ylim(0,40)
       ax2[0,0].set_title('Res', fontsize=20)
       ax2[0,0].set_ylabel('Electricity demand [kWh]', fontsize=fontsz)
           
       ax2[0,1].set_ylim(0,40)
       ax2[0,1].set_title('Ter', fontsize=20)

       ax2[1,0].set_ylim(0,40)
       ax2[1,0].set_title('Com', fontsize=20)
       ax2[1,0].set_xlabel('Hour [h]', fontsize=fontsz)
       ax2[1,0].set_ylabel('Electricity demand [kWh]', fontsize=fontsz)

       ax2[1,1].set_ylim(0,200)
       ax2[1,1].set_xlabel('Hour [h]', fontsize=fontsz)
       ax2[1,1].set_title('Pub', fontsize=20)

       plt.savefig(os.path.join(PATH_clustering,f"Electricity demand_{opt_K}.jpg"))

       # Heating demands
       fig3, ax3 = plt.subplots(figsize=(12, 10), nrows=2, ncols=2)
       fig3.tight_layout(pad = 3)
       fig3.subplots_adjust(left=0.07, right=0.98, top=0.90, bottom=0.06)
       plt.suptitle("Representative scenarios of heating demand", fontsize=24, fontweight='bold')
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

       for yi in range(opt_K):
           ax3[0,0].plot(hours, train_c[yi,:][168:192], "b", linewidth=1)
           ax3[0,1].plot(hours, train_c[yi,:][192:216], "b", linewidth=1)
           ax3[1,0].plot(hours, train_c[yi,:][216:240], "b", linewidth=1)
           ax3[1,1].plot(hours, train_c[yi,:][240:264], "b", linewidth=1)

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

       plt.savefig(os.path.join(PATH_clustering,f"Heating demand_{opt_K}.jpg"))
     
       return index_new, index_representative_r_new, number_elements_each_cluster_new, \
       TS, d_rts