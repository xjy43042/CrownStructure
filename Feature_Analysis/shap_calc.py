# -*- coding: utf-8 -*-
"""SHAP analysis for crown width/depth prediction models"""
import numpy as np
import pandas as pd
import shap
import pickle
from joblib import Parallel, delayed
import time
from functools import partial
import matplotlib.pyplot as plt
import psutil
import os

# Global plotting settings
plt.rcParams.update({
    'axes.labelsize': 12,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'figure.dpi': 600,
    'figure.figsize': (6, 4.5),
    'savefig.bbox': 'tight'
})

# Feature mapping and PFT definitions
FEATURE_MAPPING = {
    'CD_m': 'crown_depth_m',
    'ETH_m': 'canopy height',
    'pft_lai': 'PFTLAI',
    'crown_width_m': 'crown_width_m',
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

PFT_LIST = ['MBD', 'MBE', 'MNE', 'TBD', 'TBE', 'BBD', 'BNE', 'BND']
PFT_LABELS = {
    'BBD': 'BoBDT', 'BND': 'BoNDT', 'BNE': 'BoNET',
    'MBD': 'TeBDT', 'MBE': 'TeBET', 'MNE': 'TeNET',
    'TBD': 'TrBDT', 'TBE': 'TrBET'
}

def automl_predict_wrapper(model, X):
    """Wrapper function for FLAML model prediction"""
    if hasattr(X, 'values'):
        X = X.values
    return model.predict(X)

def safe_shap_calculation(explainer, X, batch_size=1000):
    """Safe chunked SHAP value calculation with memory monitoring"""
    chunks = [X.iloc[i:i+batch_size] for i in range(0, len(X), batch_size)]
    results = []
    
    for chunk in chunks:
        if psutil.virtual_memory().percent > 85:
            raise MemoryError("Memory usage exceeds 85%, aborting")
        results.append(explainer.shap_values(chunk, check_additivity=False))
    
    return np.concatenate(results)

def process_target(pft, target):
    """Process SHAP analysis for single PFT under target (CW/CD)"""
    print(f"Processing {PFT_LABELS[pft]} - {target}")
    
    # 1. Data preparation
    bio_features = [f'bio{num:02d}' for num in range(1, 21) if num not in [3,4,5,6,10,11,13,14,16,17]] # removed highly correlated features
    features = ['height_m', 'lai'] + bio_features + ['sw_avg', 'wind_max', 'pctt', 'burn']
    feature_names = [FEATURE_MAPPING.get(name, name) for name in features]
    train_path = f'/tera10/yuanhua/xiangjy/factors_new/add/mergedata/{target}_{pft}_train_new.csv'
    train_df = pd.read_csv(train_path)
    train_df.rename(columns=FEATURE_MAPPING, inplace=True)
    X_train = train_df[feature_names]
    
    # 2. Load trained model
    model_path = f'/tera10/yuanhua/xiangjy/factors_new/train_new/{target}/rm1/output/{PFT_LABELS[pft]}_autoRF.pickle'
    with open(model_path, "rb") as f:
        model = pickle.load(f)
    
    # 3. Prepare background samples
    background = shap.kmeans(X_train.values, min(10, len(X_train)))
    
    # 4. Initialize explainer
    explainer = shap.KernelExplainer(
        partial(automl_predict_wrapper, model),
        background,
        feature_names=feature_names
    )
    
    # 5. Calculate SHAP values
    sample_size = min(5000, len(X_train)) 
    X_sample = X_train.iloc[:sample_size]
    shap_values = safe_shap_calculation(explainer, X_sample, batch_size=1000)
    
    # 6. Save results
    result = {
        'shap_values': shap_values,
        'X_train': X_sample,
        'feature_names': feature_names,
        'expected_value': explainer.expected_value,
        'pft': pft,
        'target': target
    }
    
    output_path = f'./output/{target}_shap_{PFT_LABELS[pft]}.pkl'
    with open(output_path, 'wb') as f:
        pickle.dump(result, f, protocol=pickle.HIGHEST_PROTOCOL)
    
    return pft, target, result

def main():
    """Main execution function"""
    start_time = time.time()
    os.makedirs('./output', exist_ok=True)
    
    # Process both CW and CD for each PFT
    tasks = [(pft, target) for pft in PFT_LIST for target in ['CW', 'CD']]
    
    with Parallel(n_jobs=4, prefer="processes") as parallel:
        results = parallel(delayed(process_target)(pft, target) for pft, target in tasks)
    
    # Print completion summary
    completion_msg = ["Completed processing:"]
    for pft, target, _ in results:
        completion_msg.append(f"- {PFT_LABELS[pft]} {target}")
    completion_msg.append(f"\nTotal time: {(time.time()-start_time)/60:.1f} minutes")
    
    print('\n'.join(completion_msg))
if __name__ == '__main__':
    main()