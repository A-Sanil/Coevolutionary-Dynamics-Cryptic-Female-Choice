#!/usr/bin/env python3
"""
Plot RSC vs Mates comparison to verify the relationship is correct.
Since RSC is now the direct lambda for Poisson, mates should track RSC closely.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import glob
import os

# Find the most recent simulation results file
csv_files = glob.glob("parallel_sim_results_*.csv")
if not csv_files:
    raise FileNotFoundError("No simulation results CSV found")

# Get the most recent file
latest_file = max(csv_files, key=os.path.getmtime)
print(f"Using simulation file: {latest_file}")

# Read data
df = pd.read_csv(latest_file)

# Convert types
df['Generation'] = df['Generation'].astype(int)
df['Rep'] = df['Rep'].astype(str)

print(f"\n=== Data Summary ===")
print(f"Total generations: {df['Generation'].max()}")
print(f"Total replicates: {df['Rep'].nunique()}")
print(f"Mean RSC range: {df['MeanRSC'].min():.3f} to {df['MeanRSC'].max():.3f}")
print(f"Mean mates range: {df['MeanMates'].min():.3f} to {df['MeanMates'].max():.3f}")

# Create output directory
os.makedirs("Plots", exist_ok=True)

# Set style
plt.style.use('seaborn-v0_8-darkgrid')
fig_size = (10, 6)

# Plot 1: RSC and Mates evolution overlaid (all replicates)
fig, axes = plt.subplots(2, 5, figsize=(20, 8))
axes = axes.flatten()

reps = sorted(df['Rep'].unique())
for i, rep in enumerate(reps):
    rep_data = df[df['Rep'] == rep]
    ax = axes[i]
    
    ax.plot(rep_data['Generation'], rep_data['MeanRSC'], 
            label='Mean RSC (Lambda)', color='darkgreen', linewidth=1.5, alpha=0.8)
    ax.plot(rep_data['Generation'], rep_data['MeanMates'], 
            label='Mean Mates', color='darkblue', linewidth=1.5, 
            linestyle='--', alpha=0.8)
    
    ax.set_title(f'Replicate {rep}', fontsize=10, fontweight='bold')
    ax.set_xlabel('Generation', fontsize=9)
    ax.set_ylabel('Value', fontsize=9)
    ax.legend(fontsize=8, loc='best')
    ax.grid(True, alpha=0.3)

plt.suptitle('RSC (Lambda) vs Mean Mates Evolution by Replicate\n'
             'Solid: RSC | Dashed: Mean Mates', 
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('Plots/rsc_vs_mates_by_replicate.png', dpi=300, bbox_inches='tight')
plt.close()

# Plot 2: Summary across all replicates (mean ± SE)
summary = df.groupby('Generation').agg({
    'MeanRSC': ['mean', 'std'],
    'MeanMates': ['mean', 'std']
}).reset_index()

summary.columns = ['Generation', 'MeanRSC', 'StdRSC', 'MeanMates', 'StdMates']
n_reps = df['Rep'].nunique()
summary['SE_RSC'] = summary['StdRSC'] / np.sqrt(n_reps)
summary['SE_Mates'] = summary['StdMates'] / np.sqrt(n_reps)

fig, ax = plt.subplots(figsize=fig_size)
ax.plot(summary['Generation'], summary['MeanRSC'], 
        label='Mean RSC (Lambda)', color='darkgreen', linewidth=2)
ax.fill_between(summary['Generation'], 
                summary['MeanRSC'] - summary['SE_RSC'],
                summary['MeanRSC'] + summary['SE_RSC'],
                color='darkgreen', alpha=0.2)

ax.plot(summary['Generation'], summary['MeanMates'], 
        label='Mean Mates', color='darkblue', linewidth=2, linestyle='--')
ax.fill_between(summary['Generation'], 
                summary['MeanMates'] - summary['SE_Mates'],
                summary['MeanMates'] + summary['SE_Mates'],
                color='darkblue', alpha=0.2)

ax.set_xlabel('Generation', fontsize=12)
ax.set_ylabel('Value', fontsize=12)
ax.set_title('RSC (Lambda) vs Mean Mates: Summary Across All Replicates\n'
             'Mean ± SE | RSC should approximate Mean Mates since RSC = Poisson lambda',
             fontsize=14, fontweight='bold')
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('Plots/rsc_vs_mates_summary.png', dpi=300, bbox_inches='tight')
plt.close()

# Plot 3: Scatter plot - RSC vs Mates (should be approximately y=x)
fig, ax = plt.subplots(figsize=(10, 7))
ax.scatter(df['MeanRSC'], df['MeanMates'], alpha=0.4, s=20, color='steelblue')

# Linear fit
z = np.polyfit(df['MeanRSC'], df['MeanMates'], 1)
p = np.poly1d(z)
x_line = np.linspace(df['MeanRSC'].min(), df['MeanRSC'].max(), 100)
ax.plot(x_line, p(x_line), "r--", linewidth=2, label=f'Linear fit: y={z[0]:.3f}x+{z[1]:.3f}')

# y=x line (theoretical)
max_val = max(df['MeanRSC'].max(), df['MeanMates'].max())
ax.plot([0, max_val], [0, max_val], "k:", linewidth=2, label='y=x (theoretical)')

ax.set_xlabel('Mean RSC (Lambda)', fontsize=12)
ax.set_ylabel('Mean Mates per Female', fontsize=12)
ax.set_title('Relationship: RSC (Lambda) vs Mean Mates\n'
             'Red line: linear fit | Black dotted: y=x (theoretical perfect match)',
             fontsize=14, fontweight='bold')
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('Plots/rsc_vs_mates_scatter.png', dpi=300, bbox_inches='tight')
plt.close()

# Plot 4: Correlation over time
corr_data = df.groupby('Generation').apply(
    lambda x: x['MeanRSC'].corr(x['MeanMates']), include_groups=False
).reset_index()
corr_data.columns = ['Generation', 'Correlation']

fig, ax = plt.subplots(figsize=fig_size)
ax.plot(corr_data['Generation'], corr_data['Correlation'], 
        color='darkgreen', linewidth=2)
ax.axhline(y=1, color='r', linestyle='--', linewidth=1.5, label='Perfect correlation (y=1)')
ax.axhline(y=0.95, color='orange', linestyle=':', linewidth=1, label='High correlation (y=0.95)')

ax.set_xlabel('Generation', fontsize=12)
ax.set_ylabel('Correlation (RSC, Mates)', fontsize=12)
ax.set_title('Correlation between RSC and Mates over Time',
             fontsize=14, fontweight='bold')
ax.set_ylim(0, 1.05)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('Plots/rsc_vs_mates_correlation.png', dpi=300, bbox_inches='tight')
plt.close()

# Plot 5: Difference (RSC - Mates) over time
df['Difference'] = df['MeanRSC'] - df['MeanMates']
diff_summary = df.groupby('Generation').agg({
    'Difference': ['mean', 'std']
}).reset_index()
diff_summary.columns = ['Generation', 'MeanDiff', 'StdDiff']
diff_summary['SE_Diff'] = diff_summary['StdDiff'] / np.sqrt(n_reps)

fig, ax = plt.subplots(figsize=fig_size)
ax.plot(diff_summary['Generation'], diff_summary['MeanDiff'], 
        color='purple', linewidth=2)
ax.fill_between(diff_summary['Generation'], 
                diff_summary['MeanDiff'] - diff_summary['SE_Diff'],
                diff_summary['MeanDiff'] + diff_summary['SE_Diff'],
                color='purple', alpha=0.2)
ax.axhline(y=0, color='black', linestyle='--', linewidth=1, label='Zero difference')

ax.set_xlabel('Generation', fontsize=12)
ax.set_ylabel('Mean Difference (RSC - Mates)', fontsize=12)
ax.set_title('Difference: RSC - Mean Mates over Time\n'
             'Should be close to 0 | Positive = RSC > Mates | Negative = RSC < Mates',
             fontsize=14, fontweight='bold')
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('Plots/rsc_vs_mates_difference.png', dpi=300, bbox_inches='tight')
plt.close()

# Statistical analysis
print("\n=== Statistical Analysis ===")
overall_cor = df['MeanRSC'].corr(df['MeanMates'])
print(f"\nOverall correlation (all data): {overall_cor:.4f}")

# Linear model
from scipy import stats
slope, intercept, r_value, p_value, std_err = stats.linregress(df['MeanRSC'], df['MeanMates'])
print(f"\nLinear model: Mates ~ RSC")
print(f"Slope: {slope:.4f} (expected ≈ 1.0)")
print(f"Intercept: {intercept:.4f} (expected ≈ 0.0)")
print(f"R-squared: {r_value**2:.4f}")
print(f"P-value: {p_value:.2e}")

# Final generation statistics
final_gen = df['Generation'].max()
final_data = df[df['Generation'] == final_gen]
print(f"\n=== Final Generation Statistics ===")
print(f"Final generation: {final_gen}")
print(f"Mean RSC (final): {final_data['MeanRSC'].mean():.4f}")
print(f"Mean Mates (final): {final_data['MeanMates'].mean():.4f}")
print(f"Mean Difference (RSC - Mates): {(final_data['MeanRSC'] - final_data['MeanMates']).mean():.4f}")
final_cor = final_data['MeanRSC'].corr(final_data['MeanMates'])
print(f"Correlation (final gen): {final_cor:.4f}")

print("\n=== Plots Saved to Plots/ directory ===")
print("1. rsc_vs_mates_by_replicate.png - Side-by-side comparison by replicate")
print("2. rsc_vs_mates_summary.png - Summary with SE across replicates")
print("3. rsc_vs_mates_scatter.png - Scatter plot with y=x reference line")
print("4. rsc_vs_mates_correlation.png - Correlation over time")
print("5. rsc_vs_mates_difference.png - Difference (RSC - Mates) over time")

