# -*- coding: utf-8 -*-
"""Process regional vegetation canopy data using random forest models."""

import numpy as np
import pandas as pd
import netCDF4 as nc4
import xarray as xr
import time
import pickle
from pathlib import Path
import os

start_time = time.time()
PFT_NAMES = ["MNE","BNE","BND","TBE","MBE","TBD","MBD","BBD"]
PFT_LABELS = {
    'BBD': 'BoBDT', 'BND': 'BoNDT', 'BNE': 'BoNET',
    'MBD': 'TeBDT', 'MBE': 'TeBET', 'MNE': 'TeNET',
    'TBD': 'TrBDT', 'TBE': 'TrBET'
}
BASE_PATH = Path('/tera10/yuanhua/xiangjy/factors_new/rawdata')

def load_height(reg):
    """Load canopy height metrics for a region.
    
    Args:
        reg: Region coordinates [lon0, lon1, lat0, lat1]
    
    Returns:
        Tuple of (mean height, 90th percentile height) arrays
    """
    fname = f'RG_{reg[0]:.0f}_{reg[1]:.0f}_{reg[2]:.0f}_{reg[3]:.0f}.CanopyHeight_500m.nc'
    with nc4.Dataset(BASE_PATH / 'CH_5x5' / fname) as ds:
        h_mean = ds['CH_mean'][:,:].flatten().reshape(-1, 1)
        h_90 = ds['CH_90'][:,:].flatten().reshape(-1, 1)
    return h_mean, h_90

def load_lai(reg): # 得到该种PFT的2001~2020气候态的月最大 PFT_LAI
    """Load and process 21-year average PFTLAI data for a region."""
    # Load base file for coordinates
    basename= 'RG_'+str(int(reg[0]))+'_'+str(int(reg[1]))+'_'+str(int(reg[2]))+'_'+str(int(reg[3]))
    BASE_PATH = Path('/tera10/yuanhua/xiangjy/factors_new/rawdata')
    with nc4.Dataset(BASE_PATH / 'PFTLAI_5x5' / f'{basename}.MOD2020.nc') as ds:
        lat=ds['lat'][:]
        lon=ds['lon'][:]
    # Initialize arrays for multi-year processing
    n_pfts = 8
    shape_4d = (21, n_pfts, 1200, 1200)
    shape_5d = (21, 12, n_pfts, 1200, 1200)
    all_years_lai = np.full(shape_5d, np.nan)
    # all_years_pct = np.full(shape_4d, np.nan)
    # Process each year (2000-2020)
    for iyear in range(2000, 2021): 
        with nc4.Dataset(BASE_PATH / 'PFTLAI_5x5' / f'{basename}.MOD{iyear}.nc') as ds:
            # pct = ds['PCT_PFT'][1:9,:,:].filled(fill_value=np.nan)
            lai_data = ds['MONTHLY_PFT_LAI'][:,1:9,:,:].filled(fill_value=np.nan)
        # all_years_pct[iyear-2000,:,:,:] = pct
        all_years_lai[iyear-2000,:,:,:,:] = lai_data

    # Calculate 21-year average of monthly max LAI
    yearly_max  = np.max(all_years_lai, axis=1) # Max over months
    avg_pftlai  = np.mean(yearly_max, axis=0) # Average over years
    return avg_pftlai, lat, lon
def load_bioclim(reg):
    """Load 20 bioclimatic variables for a region."""
    fname = f'RG_{reg[0]:.0f}_{reg[1]:.0f}_{reg[2]:.0f}_{reg[3]:.0f}_bio.nc'
    with nc4.Dataset(BASE_PATH / 'BioClim_5x5' / fname) as ds:
        bios = np.column_stack([ds[f'bio{i:02d}'][:].flatten() for i in range(1, 21)])
    return bios

def load_sw(reg):
    """Load and average shortwave radiation data."""
    fname = f'RG_{reg[0]:.0f}_{reg[1]:.0f}_{reg[2]:.0f}_{reg[3]:.0f}.SW.nc'
    with nc4.Dataset(BASE_PATH / 'SW_5x5' / fname) as ds:
        sw = np.array([ds['Band1'][t,:,:] for t in range(12)])
    return np.nanmean(sw, axis=0).reshape(-1, 1)

def load_wind(reg):
    """Load and get max wind speed."""
    fname = f'RG_{reg[0]:.0f}_{reg[1]:.0f}_{reg[2]:.0f}_{reg[3]:.0f}.Wind.nc'
    with nc4.Dataset(BASE_PATH / 'Wind_5x5' / fname) as ds:
        wind = np.array([ds['Band1'][t,:,:] for t in range(12)])
    return np.nanmax(wind, axis=0).reshape(-1, 1)

def load_pctt(reg):
    """Load tree cover percentage."""
    fname = f'RG_{reg[0]:.0f}_{reg[1]:.0f}_{reg[2]:.0f}_{reg[3]:.0f}.RAW2020.nc'
    with nc4.Dataset(BASE_PATH / 'PCTT_5x5' / fname) as ds:
        return ds['PCTT'][:].flatten().reshape(-1, 1)

def load_burn_area(reg):
    """Load burned area data."""
    fname = f'RG_{reg[0]:.0f}_{reg[1]:.0f}_{reg[2]:.0f}_{reg[3]:.0f}_burned_area.nc'
    with nc4.Dataset(BASE_PATH / 'BurnArea_5x5' / fname) as ds:
        return ds['burned_area'][0,:,:].flatten().reshape(-1, 1)

def process_pft(pft_idx, features, model_path):
    """Process predictions for a single PFT."""
    with open(model_path, 'rb') as f:
        model = pickle.load(f)
    
    # Handle missing values
    valid_mask = ~np.isnan(features).any(axis=1)
    pred = np.full(features.shape[0], -999.0)
    
    if np.any(valid_mask):
        pred[valid_mask] = model.predict(features[valid_mask])
    
    return pred.reshape(1200, 1200)

def process_region(ireg):
    """Process canopy structure predictions for a region."""
    # Load all input data
    h_mean, h_90 = load_height(ireg)
    lai, lat, lon = load_lai(ireg)
    bioclim = load_bioclim(ireg)
    sw = load_sw(ireg)
    wind = load_wind(ireg)
    pctt = load_pctt(ireg)
    burn = load_burn_area(ireg)

    # Initialize output arrays
    shape = (len(PFT_NAMES), 1200, 1200)
    cw_mean = np.full(shape, -999.0)
    cd_mean = np.full(shape, -999.0)
    cw_90 = np.full(shape, -999.0)
    cd_90 = np.full(shape, -999.0)

    # Process each PFT
    for i, pft in enumerate(PFT_NAMES):
        # Prepare features
        pft_lai = lai[i]
        pft_lai_flat = lai[i].reshape(-1, 1)
        features_mean = np.hstack([h_mean, pft_lai_flat, bioclim, sw, wind, pctt, burn])
        features_90 = np.hstack([h_90, pft_lai_flat, bioclim, sw, wind, pctt, burn])
        
        # Model paths
        model_dir = '/tera10/yuanhua/xiangjy/factors_new/train_new'
        cw_path = f'{model_dir}/CW/all/output/{PFT_LABELS[pft]}_autoRF.pickle'
        cd_path = f'{model_dir}/CD/all/output/{PFT_LABELS[pft]}_autoRF.pickle'
        
        # Make predictions
        cw_mean[i] = process_pft(i, features_mean, cw_path)
        cd_mean[i] = process_pft(i, features_mean, cd_path)
        cw_90[i] = process_pft(i, features_90, cw_path)
        cd_90[i] = process_pft(i, features_90, cd_path)

        # Handle LAI=0 grid cells (core modification)
        zero_mask = (pft_lai == 0)
        cw_mean[i][zero_mask] = -999.0
        cd_mean[i][zero_mask] = -999.0
        cw_90[i][zero_mask] = -999.0
        cd_90[i][zero_mask] = -999.0
    
    return cw_mean, cd_mean, cw_90, cd_90, lat, lon

def save_results(results, lat, lon, ireg, suffix):
    """Save canopy structure predictions to NetCDF with full metadata.
    
    Args:
        results: Tuple of (CW, CD, AR) prediction arrays
        lat: Latitude coordinates
        lon: Longitude coordinates
        ireg: Region coordinates [lon0, lon1, lat0, lat1]
        suffix: File suffix ('CHavg' or 'CH90')
    """
    CW, CD, AR = results[:3]    
    pftids = np.arange(1, 9, dtype=np.int32)
    pft_attrs = {
        'long_name': 'Index of Plant Functional Type, only including 8 tree types',
        'PFT_1': 'Temperate Needleleaf Evergreen Tree',
        'PFT_2': 'Boreal Needleleaf Evergreen Tree',
        'PFT_3': 'Boreal Needleleaf Deciduous Tree',
        'PFT_4': 'Tropical Broadleaf Evergreen Tree',
        'PFT_5': 'Temperate Broadleaf Evergreen Tree',
        'PFT_6': 'Tropical Broadleaf Deciduous Tree',
        'PFT_7': 'Temperate Broadleaf Deciduous Tree',
        'PFT_8': 'Boreal Broadleaf Deciduous Tree'
    }
    
    ds = xr.Dataset(
        data_vars={
            'CROWN_WIDTH': (
                ['PFT', 'lat', 'lon'], 
                CW.astype(np.float32),
                {'long_name': f'crown width predicted by the {suffix.replace("CH", "")} canopy-top height'}
            ),
            'CROWN_DEPTH': (
                ['PFT', 'lat', 'lon'],
                CD.astype(np.float32),
                {'long_name': f'crown depth predicted by the {suffix.replace("CH", "")} canopy-top height'}
            ),
            'ASPECT_RATIO': (
                ['PFT', 'lat', 'lon'],
                AR.astype(np.float32),
                {'long_name': f'crown depth_to_width ratio predicted by the {suffix.replace("CH", "")} canopy-top height'}
            )
        },
        coords={
            'PFT': (['PFT'], pftids, pft_attrs),
            'lat': (['lat'], lat.astype(np.float64), {
                'long_name': 'Latitude',
                'units': 'degrees_north'
            }),
            'lon': (['lon'], lon.astype(np.float64), {
                'long_name': 'Longitude', 
                'units': 'degrees_east'
            })
        }
    )

    # Encoding settings
    encoding = {
        var: {
            'dtype': 'float32',
            '_FillValue': -999.0,
            'zlib': True,
            'complevel': 6
        } for var in ds.data_vars
    }

    # Generate filename
    filename = f'./output/RG_{ireg[0]:.0f}_{ireg[1]:.0f}_{ireg[2]:.0f}_{ireg[3]:.0f}.CanopyStructure_15s_{suffix}.nc'

    # Save file
    ds.to_netcdf(filename, encoding=encoding)
    ds.close()

def main():
    """Main processing workflow."""
    regions = np.loadtxt('reg_5x5')
    os.makedirs('./output/', exist_ok=True)
    for ireg in regions:
        cw_mean, cd_mean, cw_90, cd_90, lat, lon = process_region(ireg)
        
        # Calculate aspect ratios
        ar_mean = np.where((cw_mean == -999.0) | (cd_mean == -999.0), -999.0, cd_mean/cw_mean)
        ar_90 = np.where((cw_90 == -999.0) | (cd_90 == -999.0), -999.0, cd_90/cw_90)
        
        # Save results
        save_results((cw_mean, cd_mean, ar_mean), lat, lon, ireg, 'CHavg')
        save_results((cw_90, cd_90, ar_90), lat, lon, ireg, 'CH90')
    
    print(f"Total execution time: {(time.time()-start_time)/3600:.2f} hours")

if __name__ == '__main__':
    main()