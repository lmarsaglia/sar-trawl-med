# Shedding light on trawlers fishing activity in the Mediterranean Sea with remote sensing data 

In this repository all analyses and data to generate all results of Marsaglia et al. 2024 can be accessed.

The overall reccomendation is to open the `r` project file in order to have the different folders and notebooks to "speak" across them.

A description of each folder is the following: 

- **analysis**
  
  1. `STI_analysis.Rmd` is the notebook  running the Spatial Temporal Autocorrelation
  2. `preprocessing - aggragate all data to grid.Rmd` contains all the process to aggregate our data to a 0.2 decimal degrees grid
  3. `results analysis.Rmd` can be used to generate all figures in the paper except Figure 3 that is generated in the processing notebook and Figure 4 and 5 that are generated in the GAM script in the model folder.
  
- **data**: contains all data that are used in the analyses and created and read by the notebooks and scripts present in the repository.
  
- **figures**: folder with all paper figures and where script will sa

- **functions**: STI function used in notebook 1 lives here.

- **model**: GAM model scripts lives here

- **shapefiles**: contains all shapefiles used in the analyses



