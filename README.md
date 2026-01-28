# Optimization of a Multi-Energy System (MES) meeting the energy demands of a Renewable Energy Community (REC) under uncertainty

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

## Overview of software interface
The software consists of a main folder including three core Python scripts, each responsible for a distinct stage of the workflow. “main.py” constitutes the software interface. “Scenario_generation_function.py” contains two Python functions that perform data processing and clustering (“scenario_generation”), as well as post-processing (“post_clustering”). “Model_REC_function.py” encompasses two Python functions that carry out the clustering-based design-operation optimization (“clustering_design_operation”) and the full timeseries operation-only optimization (“fullts_operation”) of the REC. Specifically, the software is designed to:
1.  Import and pre-process raw input timeseries of weather data (e.g., solar irradiance and ambient temperature), electricity zonal market prices and user energy demand, as well as apply normalization and cleaning routines (through “scenario_generation” function).
2.	Generate 49 different sets of representative daily stochastic scenarios of PV capacity factor, ambient temperature and zonal electricity prices by applying the K-means clustering algorithm to a large training dataset (2005-2017) of historical timeseries (through “scenario_generation” function);
3.	Evaluate the clustering quality of the representative days using different metrics such as the Average Silhouette Value (ASV), Sum of Squared Errors (SSE), Davies-Bouldin Index (DBI) and Dunn Index (DI) (through “scenario_generation” function). In addition, plots of timeseries divided into clusters and representative scenarios are obtained (through “post_clustering” function);
4.	Formulate the REC model and solve the clustering-based design-operation optimization for each set of representative scenarios in the training dataset (through “clustering_design_operation” function). Each clustering-based solution consists of optimal investment decisions (i.e., technology sizing) and operational strategies (e.g., dispatch scheduling) subject to technical and economic constraints, including non-linearities due to the definition of the shared energy incentive;
5.	Assess the optimal unit sizes of the clustering-based design solutions in a testing dataset (2018-2023), using a full timeseries operation-only optimization that considers higher prices than in the training dataset (through “fullts_operation” function).
It is worth noting that “main.py” includes five sections that can be run independently. The first section imports the two Python files “Scenario_generation_function.py” and “Model_REC_function.py”, and enables the user to select input data by specifying the values of the following input parameters: 
-	“loc” (type: string), the location where renewable plants are installed (“Padova”, “Palermo” or “Trieste” by default);
-	“plot_ts” (type: string), deciding whether to create plots of the full timeseries of PV capacity factor, ambient temperature and zonal electricity prices (“yes” or “no”); 
-	“opt_K” (type: integer number), the “best” number of clusters for plots and the “best” representative days associated with the clustering-based design solution to be tested by the full timeseries operation-only optimization.

The above parameters are input to the Python functions implemented within “main.py”. In addition, the default input parameters are:
-	“year_start”, the first year of the available timeseries dataset (default: 2005);
-	“year_end”, the last year of the available timeseries dataset (default: 2023);
-	“seasons”, which decides whether to consider the division of the available dataset into seasons or not (0: no, 1: yes, default: 0);
-	“extreme_method”, the criterion used to include extreme zonal price days into each set of representative days ("A": adding criterion, "R": replacing criterion, default: “A”). According to the adding criterion, extreme days become the representative days of new clusters, after an iterative process that reassigns all elements to clusters. With the replacing criterion, the extreme days substitute the representative days of the clusters they belong to.

In addition, the main folder contains several subfolders for organizing input datasets, intermediate outputs, and final results. The main subfolders are:
-	“Input_data”, which contains the “TechnoEconomicData.xlsx” file with techno-economic parameters of technologies, as well as the “Demands” and “Prices” subfolders and “.csv” files for solar PV data and ambient temperature in defined locations;
-	“Plots” containing the “Timeseries” subfolder (divided into other subfolders for different locations) with plots of the full timeseries of PV capacity factor, ambient temperature and zonal electricity prices;
-	“Clustering” (divided into other subfolders for different locations) including clustering results, plots of metrics calculated to evaluate clustering quality, and plots of representative days of PV capacity factor, ambient temperature, zonal electricity prices and user energy demands;
-	“Results” (divided into other subfolders for different locations) containing results obtained from the clustering-based design-operation optimization and the full timeseries operation-only optimization.
The techno-economic parameters of technologies are investment and operational costs, efficiencies and lifetime. Input timeseries of solar PV capacity factor and ambient temperature are retrieved from PVGIS 5.3. “Demands” and “Prices” subfolders contain timeseries of various user energy demands, collected from an online database, and zonal electricity prices from the database of “Gestore Mercati Energetici” (GME), the Italian energy market operator, respectively. Among the intermediate outputs, the software (in particular “Scenario_generation_function.py”) produces several “.jpg” plots showing the timeseries categorized into a desired number of clusters (i.e., through “opt_K” input parameter) and the results of the clustering quality analysis. Moreover, “.npz” files are generated to save the clustering results, which are also used as input to “Model_REC_function.py”. The optimization results include optimal sizes for various technologies (e.g., PV, battery, heat pump) and total annualized investment and operational costs. These results are saved as “.npz”, “.xlsx” and “.jpg” files.

## Functions of the software
This section provides a practical guide to the Python functions implemented in the software, explaining how and when to use them, what inputs and outputs to expect, and some recommended best practices from a user’s operational perspective. 
The software is designed to work correctly when the input parameters “year_start”, “year_end”, “seasons” and “extreme_method” are set to their default values, as specified in the previous section. The first section of “main.py” must always be run when opening the software to import the two Python scripts and select the input parameters that define the location (“loc”), whether to plot the full timeseries (“plot_ts”), and the “best” number of clusters (“opt_K”). Note that the “best” number of clusters is not known until K-means clustering has been applied and its quality analyzed. Hence, for the first run of the software, it is suggested to fix an “opt_K” value between 4 and 20.
The second section generates the clustering outcomes, which are saved within the “Clustering” folder. Therefore, it is recommended to run the second section when the clustering results for the chosen location have not been obtained yet. Moreover, the second section also produces “.jpg” plots of the multi-year timeseries of zonal prices, PV capacity factor and ambient temperature, which are saved in the “Plots” folder. The third section is mainly used to produce the “.jpg” plots (within the “Clustering” folder) showing the price timeseries classified into a fixed “best” number of clusters, as well as the representative days for PV capacity factor, ambient temperature, zonal electricity prices and user energy demands. The fourth and fifth sections carry out the clustering-based design-operation optimization, which requires the results from the second section, and the full timeseries operation-only optimization, respectively. Optimization results are saved in the “Results” folder. As with the second section, it is recommended that the fourth section is run if the associated results are not yet available. Additionally, note that the full timeseries operation-only optimization is used to test the optimal sizes of the “best” clustering based-design solution. Consequently, the fifth section can only be run once the clustering-based solutions, saved as “.npz” files, have been obtained. The operational details of the four Python functions within the software are summarized in the following.
