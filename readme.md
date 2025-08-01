# Generating Vegetation Canopy Structure

## Overview
This repository contains the processing pipeline for upscaling vegetation canopy structure parameters (crown depth/width/ratio) from site-scale measurements to global 500m resolution dataset, supporting the manuscript:
> "A Global 500 m Resolution Tree Crown Structure Dataset"

## Code Organization

The code is organized into four main functional directories:

1. **Model Training** (`Model_Training/`)  
   Develops FLAML-based Random Forest models to predict crown structure with site-scale measurements.

2. **Feature Analysis** (`Feature_Analysis/`)  
   Performs SHAP analysis to evaluate and visualize feature importance.

3. **Canopy Height Processing** (`Canopy_Height_Processing/`)  
   Processes canopy height data from 1m to 500m resolution.

4. **Global Prediction** (`Global_Prediction/`)  
   Generates regional predictions at 500m and aggregates to global 0.5° scale.


## Directory Structure
```bash
├── Model_Training/                 
│   └── autoRF.py                   # Automated Random Forest training using FLAML
├── Feature_Analysis/               
│   ├── shap_calc.py                # SHAP value computation
│   └── plot_importance.py          # Feature ranking visualization
├── Canopy_Height_Processing/       
│   ├── tif2nc.py                   # GeoTIFF→NetCDF conversion (1m→30m)
│   └── aggr_CH_500m.F90            # Spatial aggregation (30m→500m)
└── Global_Prediction/              
    ├── predict_regional.py         # Generate 5°×5° regional crown structure predictions
    └── aggr_globe.F90              # aggregate regional predictions to global 0.5° grid
```
<br>

## Usage
Code is mostly written in Python 3.6+ with dependencies including sklearn, flaml, numpy, pandas, xarray and matplotlib. 

High-performance spatial aggregation is achieved through ​​Fortran 90+​​ implementations compiled with NetCDF library support.
