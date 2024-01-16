# Optimization with and without uncertainty of a Multi-Energy System (MES) meeting the energy demands of a Renewable Energy Community (REC)


## Description of the system 
### REC model:
  - Members with electricity and heating demands:
    - residential, house (res)
    - tertiary, small office (ter)
    - commercial, small restaurant (com)
    - public, primary school (pub)
  - Individual self-consumption for each member
  - Virtual energy sharing among members through the electric grid
  - Heating demand satisfied by each member independently
### Energy conversion and storage units:
  - Photovoltaic (PV)
  - Boiler (BOIL), fuelled by natural gas
  - Air-water heat pump (HP)
  - Thermal energy storage (TES): hot water tank
  - Electrical energy storage (EES): lithium battery
### How are these units distributed among the members
  - res (prosumer): PV, EES, HP, TES
  - ter (consumer): BOIL
  - com (consumer): BOIL
  - pub (prosumer): PV, EES, HP, TES

## Optimization models of the system
### The models are solved in each year of a training dataset (2005-2014)
### The model without uncertainty is based on 4 average seasonal days in each year
### Model with uncertainty:
  - solar irradiance and ambient temperature are the uncertain parameters
  - from 2 to 7 typical days for each season in each year
  - Example:
    - First case: 2 typical days for each season --> 8 scenarios in one year
    - Last case: 7 typical days for each season --> 28 scenarios in one year

## Instructions to use the software

### The software includes:
  - the script "Model_ML2.py" developed in Python, containing the optimization models with and without uncertainty, which are solved by the Gurobi optimizer
  - input data in the "Input_data" folder
  - results in the "Results" folder
  - this README.md file

### Insert the correct paths of your folders containing input data or results in the script "Model_ML2"
  - line 56: to import techno-economic data from PATH_data = "C:/.../Input_data/TechnoEconomicData.xlsx"
  - line 793: to import timeseries from PATH_clustering = "C:/.../Input_data/Timeseries/"+address
    - "address" is the city chosen in the section "AVERAGE SEASONAL DAYS AND STOCHASTIC SCENARIOS" (see the sections afterwards)
  - line 1053: to save results in PATH = "C:/.../Results/"+address
  - line 1070: to export results without uncertainty to PATH_export_det = "C:/.../Results/"+address+"/Deterministic_solution_"+address+".xlsx"
  - line 1105: to export results with uncertainty to PATH_export_sto = "C:/.../Results/"+address+"/Stochastic_solutions_"+address+".xlsx"
  - line 1133: to save plots with uncertainty in PATH_plots_sto = "C:/.../Results/"+address+"/Plots/Stochastic_scenarios"
  - line 1259: to save plots without uncertainty in PATH_plots_det = "C:/.../Results/"+address+"/Plots/Average_seasonal_days"

### In running the script "Model_ML2", the user can choose:
  - a city (Belluno, Padova, Palermo) as location for the weather parameters (solar irradiance and ambient temperature)
  - a year (from 2005 to 2014) and a number of typical days (from 2 to 7) for each season to obtain plots 

### The script "Model_ML2" is divided in sections that can be run independently
  - "MODEL": general comments on the models
  - "PYTHON PACKAGES": import Python packages
  - "TECHNO-ECONOMIC PARAMETERS": import techno-economic parameters of the conversion and storage units, price of natural gas and financial parameters
  - "DESIGN-OPERATION OPTIMIZATION MODEL": model with and without uncertainty built as a Python function
    - it includes the decision variables, the constraints, the objective function (total investment and operational costs actualized to one year) and the arrays of the results to be saved
    - the number of days "K" in one year is the input to the Python function
    - the optimal values of the objective function and some decision variables are the output of the Python function
  - "AVERAGE SEASONAL DAYS AND STOCHASTIC SCENARIOS":  
    - choose the city (variable "address") by writing "Belluno", "Padova", or "Palermo" in the console 
    - import the average seasonal days and the typical days of each season for each year of the training dataset
  - "SOLVING THE MODEL WITHOUT UNCERTAINTY": solve the model without uncertainty for 4 average seasonal days (K=4) for each year of the training dataset
  - "SOLVING THE MODEL WITH UNCERTAINTY": solve the model with uncertainty based on the same number of typical days for each season (2,3,4,5,6,7), resulting in different stochastic scenarios for each year of the training dataset (K=8,12,16,20,24,28)
  - "SAVE AND EXPORT RESULTS": 
    - save the main results with and without uncertainty in the npz files "Stochastic_solutions_+address.npz" and "Deterministic_solution_+address.npz", respectively
    - export the optimal total costs, investment costs and operational costs with and without uncertainty to the excel files "Stochastic_solutions_+address.xlsx" and "Deterministic_solution_+address.xlsx", respectively
  - "PLOTS": obtain the plots of solar irradiance, ambient temperature, electricity and heating demands, and grid sale price with and without uncertainty
    - choose a year (variable "Y") by writing a value in the range 2005-2014 and a number of typical days for each season (variable "C") by writing a value in the range 2-7 in the console
