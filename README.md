# Generating Vegetation Canopy Structure

## Overview
This repository contains the processing pipeline for mapping vegetation canopy structure parameters (crown depth/width/ratio) from site-scale measurements to a global 500 m resolution dataset, supporting the manuscript:
> "A Global 500 m Resolution Tree Crown Structure Dataset"

## Code Organization

The code is organized into four main functional directories:

1. **Model Training** (`Model_Training/`)
   Develops FLAML-based Random Forest models to predict crown structure using site-scale measurements.

2. **Feature Analysis** (`Feature_Analysis/`)
   Performs SHAP analysis to evaluate and visualize feature importance.

3. **Canopy Height Processing** (`Canopy_Height_Processing/`)
   Processes canopy height data from 1 m to 500 m resolution.

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
│   ├── tif2nc.py                   # GeoTIFF → NetCDF conversion (1 m → 30 m)
│   └── aggr_CH_500m.F90            # Spatial aggregation (30 m → 500 m)
└── Global_Prediction/
    ├── predict_regional.py         # Generate 5°×5° regional crown structure predictions
    └── aggr_globe.F90              # Aggregate regional predictions to global 0.5° grid
```
<br>

## Usage
Code is primarily written in Python 3.6+ with dependencies including sklearn, flaml, numpy, pandas, xarray and matplotlib.

High-performance spatial aggregation is achieved through ​​Fortran 90+​​ implementations compiled with NetCDF library support.
