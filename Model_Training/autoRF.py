# -*- coding: utf-8 -*-
"""AutoML pipeline for crown width/depth training with FLAML"""
import numpy as np
import pandas as pd
import joblib
import time
import os
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from flaml import AutoML
from sklearn.metrics import (r2_score, 
                            mean_squared_error, 
                            mean_absolute_error)
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cloudpickle

# Global plotting settings
plt.rcParams.update({
    'axes.labelsize': 12,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'figure.dpi': 600,
    'figure.figsize': (5, 3.5),
    'savefig.bbox': 'tight'
})

# Constants
FEATURE_MAPPING = {
    'CD_m': 'crown_depth_m',
    'ETH_m': 'canopy height',
    'pft_lai': 'PFTLAI',
    'crown_width_m': 'crown_width_m',
    'crown_depth_m': 'crown_depth_m',
    'height_m': 'canopy height',
    'lai': 'PFTLAI',
    'bio01': 'Annual Mean T',
    'bio02': 'Mean Diurnal T Range',
    'bio03': 'Isothermality',
    'bio04': 'T Seasonality',
    'bio05': 'Max T of Warmest Month',
    'bio06': 'Min T of Coldest Month',
    'bio07': 'T Annual Range',
    'bio08': 'Mean T of Wettest Qtr',
    'bio09': 'Mean T of Driest Qtr',
    'bio10': 'Mean T of Warmest Qtr',
    'bio11': 'Mean T of Coldest Qtr',
    'bio12': 'Annual Prcp',
    'bio13': 'Prcp of Wettest Month',
    'bio14': 'Prcp of Driest Month',
    'bio15': 'Prcp Seasonality',
    'bio16': 'Prcp of Wettest Qtr',
    'bio17': 'Prcp of Driest Qtr',
    'bio18': 'Prcp of Warmest Qtr',
    'bio19': 'Prcp of Coldest Qtr',
    'bio20': 'Elevation',
    'wind_avg': 'Mean wind speed',
    'sw_avg': 'Mean shortwave radiation',
    'wind_max': 'Max wind speed',
    'sw_max': 'Max shortwave radiation',
    'pctt':'Tree cover',
    'burn': 'burned_area'
}

PFT_LIST = ['BBD', 'BND', 'BNE', 'MBD', 'MBE', 'MNE', 'TBD', 'TBE']

PFT_LABELS = {
    'BBD': 'BoBDT', 'BND': 'BoNDT', 'BNE': 'BoNET',
    'MBD': 'TeBDT', 'MBE': 'TeBET', 'MNE': 'TeNET',
    'TBD': 'TrBDT', 'TBE': 'TrBET'
}
TARGET_FULL = {'CW': 'crown_width_m', 'CD': 'crown_depth_m'}

def calculate_metrics(true, pred):
    """Calculate evaluation metrics"""
    return {
        'R2': r2_score(true, pred),
        'RMSE': np.sqrt(mean_squared_error(true, pred)),
        'MAE': mean_absolute_error(true, pred),
        'MBE': np.mean(pred - true)
    }

def plot_scatter(obs, pred, title, xlabel, ylabel, filename):
    """Generate scatter plot with density coloring"""
    obs, pred = np.asarray(obs), np.asarray(pred)
    
    # Calculate density and sort
    density = gaussian_kde(np.vstack([obs, pred]))(np.vstack([obs, pred]))
    idx = density.argsort()
    obs, pred, density = obs[idx], pred[idx], density[idx]
    
    # Determine axis limits
    max_val = max(obs.max(), pred.max())
    max_val = np.ceil(max_val) + 1
    
    # Create plot
    fig, ax = plt.subplots(figsize=(5, 3.5), dpi=600)
    sc = ax.scatter(obs, pred, c=density, s=10, cmap='Spectral_r')
    
    # Formatting
    ax.set(xlim=(0, max_val), ylim=(0, max_val),
           xlabel=xlabel, ylabel=ylabel, title=title)
    if max_val > 34:
        step = 6
    elif max_val > 24:
        step = 4
    else:
        step = 2
    ax.set_xticks(np.arange(0, max_val + 1, step))
    ax.set_yticks(np.arange(0, max_val + 1, step))
    ax.plot([0, max_val], [0, max_val], 'k--', lw=0.8)
    ax.set_aspect('equal', adjustable='box', anchor='C')

    # Add metrics
    ax.text(1.03, 1.01, f"N = {len(obs)}", transform=ax.transAxes, fontsize=10, va='bottom', ha='right')
    metrics = calculate_metrics(obs, pred)
    ax.text(0.05, 0.92, f"R² = {metrics['R2']:.2f}", transform=ax.transAxes)
    ax.text(0.05, 0.84, f"RMSE = {metrics['RMSE']:.2f}", transform=ax.transAxes)
    ax.text(0.05, 0.76, f"MBE = {metrics['MBE']:.2f}", transform=ax.transAxes)
   
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cbar = fig.colorbar(sc, cax=cax, label='density')

    plt.tight_layout()
    plt.savefig(filename, dpi=600, format="jpeg")
    plt.close()

def add_predictions(df, features, model, suffix=""):
    """Add model predictions to DataFrame"""
    df = df.copy()
    df[f'Predicted{suffix}'] = model.predict(df[features].values)
    return df

def train_model(X_train, y_train, pft_label, target, n_jobs=72):
    """Train AutoML model for given target"""
    automl = AutoML()
    setting = {
        "time_budget": 60,  # 30 minutes
        "metric": 'rmse', 
        "task": 'regression',
        "log_file_name": f'./output/{target}_{pft_label}_autoRF.log',
        "estimator_list": ['rf'], 
        "n_jobs": n_jobs,
        "n_splits": 10,
        "verbose": 2,
        "seed": 42,
        "early_stop": True
    }
    automl.fit(X_train, y_train, **setting)
    return automl

def process_target(pft, target, features, feature_names):
    """Process model training for single target"""
    print(f"\n==== Processing {PFT_LABELS[pft]} - {target} ====")
    start_time = time.time()
    
    # Load data
    train_df = pd.read_csv(f'/tera10/yuanhua/xiangjy/factors_new/add/mergedata/{target}_{pft}_train_new.csv')
    test_df = pd.read_csv(f'/tera10/yuanhua/xiangjy/factors_new/add/mergedata/{target}_{pft}_test_new.csv')
    train_df.rename(columns=FEATURE_MAPPING, inplace=True)
    test_df.rename(columns=FEATURE_MAPPING, inplace=True)
    
    X_train = train_df[feature_names].values
    y_train = train_df[TARGET_FULL[target]].values
    X_test = test_df[feature_names].values
    y_test = test_df[TARGET_FULL[target]].values
    
    print(f"Train: {X_train.shape[0]} samples, {X_train.shape[1]} features")
    print(f"Test: {X_test.shape[0]} samples")
    print(f"Target: {TARGET_FULL[target]}")
    print("Features:")
    features_per_line = 4
    for i in range(0, len(features), features_per_line):
        print("   ", ", ".join(features[i:i+features_per_line]))

    # Train model
    automl = train_model(X_train, y_train, PFT_LABELS[pft], target)
    
    # Save model
    with open(f'./output/{target}_{PFT_LABELS[pft]}_autoRF.pickle', 'wb') as f:
        cloudpickle.dump(automl, f)
    print(f"Saved model to ./output/{target}_{PFT_LABELS[pft]}_autoRF.pickle")
    
    # Evaluate
    y_pred_train = automl.predict(X_train)
    y_pred_test = automl.predict(X_test)
    
    # Save predictions
    full_df = pd.concat([
        add_predictions(train_df, feature_names, automl, "_train").assign(Dataset='train'),
        add_predictions(test_df, feature_names, automl, "_test").assign(Dataset='test')
    ])
    full_df.to_csv(f'./output/{target}_predictions_{PFT_LABELS[pft]}.csv', index=False)
    
    # Save metrics
    metrics = {
        'train': calculate_metrics(y_train, y_pred_train),
        'test': calculate_metrics(y_test, y_pred_test),
        'full': calculate_metrics(np.concatenate([y_train, y_test]), 
                                np.concatenate([y_pred_train, y_pred_test]))
    }
    pd.DataFrame(metrics).to_csv(f'./output/{target}_metrics_{PFT_LABELS[pft]}.csv')
    
    # Visualization
    plot_scatter(y_test, y_pred_test, 
                f'{PFT_LABELS[pft]} - {target}', 
                f'{target} observed (m)', 
                f'{target} predicted (m)', 
                f'{target}_scatter_test_{PFT_LABELS[pft]}.jpeg')
    
    print(f"Completed in {(time.time()-start_time)/60:.1f} minutes")

def run_automl_pipeline():
    """Main pipeline execution"""
    os.makedirs('./output', exist_ok=True)
    # Define features
    bio_features = [f'bio{num:02d}' for num in range(1, 21)]
    features = ['height_m', 'lai'] + bio_features + ['sw_avg','wind_max','pctt','burn']
    feature_names = [FEATURE_MAPPING.get(name, name) for name in features]
    
    # Process both CW and CD for each PFT
    for pft in PFT_LIST:
        for target in ['CW', 'CD']:
            process_target(pft, target, features, feature_names)

if __name__ == '__main__':
    run_automl_pipeline()