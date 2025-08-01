# -*- coding: utf-8 -*-
"""Simplified SHAP visualization for crown width/depth analysis"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import string
import os
from pylab import rcParams

# Global plotting settings
plt.rcParams.update({
    'axes.labelsize': 18,
    'font.size': 21,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'figure.dpi': 600,
    'savefig.bbox': 'tight',
    'axes.edgecolor': "#2A2929",
    'xtick.color': "#2A2929",
    'ytick.color': "#2A2929"
})

# Feature color mapping
FEATURE_COLORS = {
    'canopy height': '#89C15F',
    'PFTLAI': '#89C15F',
    'Annual Mean T': '#FF97B5',
    'Mean Diurnal T Range': '#FF97B5',
    'T Annual Range': '#FF97B5',
    'Mean T of Wettest Qtr': '#FF97B5',
    'Mean T of Driest Qtr': '#FF97B5',
    'Annual Prcp': '#7EB0DE',
    'Prcp Seasonality': '#7EB0DE',
    'Prcp of Warmest Qtr': '#7EB0DE',
    'Prcp of Coldest Qtr': '#7EB0DE',
    'Prcp of Direst Qtr': '#7EB0DE',
    'Elevation': 'brown',
    'Mean shortwave radiation': "#FFBC57",
    'Max wind speed': "#9E5CFA",
    'Tree cover':'#89C15F',
    'burned_area': "#B71500",
    'Others': "#ABAAAA"
}

# PFT definitions
PFT_LIST = ['BBD', 'BND', 'BNE', 'MBD', 'MBE', 'MNE', 'TBD', 'TBE']
PFT_LABELS = {
    'BBD': 'BoBDT', 'BND': 'BoNDT', 'BNE': 'BoNET',
    'MBD': 'TeBDT', 'MBE': 'TeBET', 'MNE': 'TeNET',
    'TBD': 'TrBDT', 'TBE': 'TrBET'
}

def plot_shap_importance(var, top_n=5):
    """Generate SHAP importance plots for given variable (CW/CD)"""
    os.makedirs('./fig', exist_ok=True)
    
    # Create 3x3 subplot grid
    fig, axes = plt.subplots(3, 3, figsize=(15, 9))
    axes = axes.flatten()
    letters = string.ascii_lowercase[:9]
    
    for i, pft in enumerate(PFT_LIST):
        print(f"Processing {PFT_LABELS[pft]} {var}")
        
        # Load SHAP data
        with open(f'/tera10/yuanhua/xiangjy/factors_new/SHAP/{var}/output/{var}_shap_{PFT_LABELS[pft]}.pkl', 'rb') as f:
            data = pickle.load(f)
        
        # Calculate global SHAP importance
        global_shap = np.abs(data['shap_values']).mean(axis=0)
        percentage_shap = (global_shap / global_shap.sum()) * 100
        
        # Create and sort importance DataFrame
        shap_df = pd.DataFrame({
            'feature': data['feature_names'],
            'importance': percentage_shap
        }).sort_values('importance', ascending=False)
        
        # Get top features and remaining
        top_features = shap_df.head(top_n)
        others_importance = shap_df.iloc[top_n:]['importance'].sum()
        plot_df = pd.concat([
            pd.DataFrame({'feature': ['Others'], 'importance': [others_importance]}),
            top_features.sort_values('importance', ascending=True)
        ])
        
        # Generate plot
        ax = axes[i]
        y_pos = np.arange(len(plot_df))
        colors = [FEATURE_COLORS.get(feat, '#ABAAAA') for feat in plot_df['feature']]
        
        bars = ax.barh(y_pos, plot_df['importance'], 
                      height=0.8, color=colors, edgecolor='k')
        
        # Format axes
        ax.set_yticks(y_pos)
        ax.set_yticklabels(plot_df['feature'], rotation=45, fontsize=16)
        ax.set_ylim(-0.5, len(plot_df)-0.5)
        ax.set_xlim(0, 100)
        ax.set_xticks(np.arange(0, 101, 20))
        
        # Add subplot label
        ax.set_title(f'({letters[i]}) {PFT_LABELS[pft]}', loc='left', fontsize=21)
        
        # Add value labels
        for j, (pos, row) in enumerate(zip(y_pos, plot_df.itertuples())):
            ha = 'right' if row.feature == 'canopy height' and row.importance > 40 else 'left'
            ax.text(row.importance + (1 if ha == 'left' else -1), 
                    pos, f"{row.importance:.1f}%", 
                    va='center', ha=ha, fontsize=18)
    
    # Clean up empty subplots
    for j in range(len(PFT_LIST), len(axes)):
        fig.delaxes(axes[j])
    
    # Add legend based on variable type
    legend_colors = {
        'CW': {
            'Intrinsic': '#89C15F',
            'Temperature': '#FF97B5',
            'Precipitation': '#7EB0DE',
            'Radiation': '#FFBC57',
            'Wind': '#9E5CFA',
            'Fire': "#B71500",
            'Others': "#ABAAAA"
        },
        'CD': {
            'Intrinsic': '#89C15F',
            'Temperature': '#FF97B5',
            'Precipitation': '#7EB0DE',
            'Radiation': '#FFBC57',
            'Wind': '#9E5CFA',
            'Fire': "#B71500",
            'Others': "#ABAAAA"
        }
    }
    
    handles = [plt.Rectangle((0,0),1,1, facecolor=color, label=label, edgecolor='k') 
               for label, color in legend_colors[var].items()]
    
    fig.legend(handles=handles, loc='lower right', 
              bbox_to_anchor=(0.9, 0.03), 
              fontsize=18)
    
    plt.subplots_adjust(wspace=0.82, hspace=0.5)
    plt.savefig(f"./fig/{var}_SHAP_Importance_top{top_n}.jpeg", 
               dpi=600, format="jpeg", bbox_inches='tight')
    plt.close()
    print(f"{var} visualization completed")

if __name__ == '__main__':
    # Process both crown width and depth
    for target_var in ['CW', 'CD']:
        plot_shap_importance(target_var, top_n=5)