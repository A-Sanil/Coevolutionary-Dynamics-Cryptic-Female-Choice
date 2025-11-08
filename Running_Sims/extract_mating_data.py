#!/usr/bin/env python3
"""
Extract only mating-related columns from simulation results CSV
"""

import pandas as pd
from pathlib import Path

# Read the full dataset
input_file = 'parallel_sim_results_2025-11-08_141821.csv'
df = pd.read_csv(input_file)

print(f"Original CSV shape: {df.shape}")
print(f"Original columns: {list(df.columns)}")

# Select only mating-related columns
mating_columns = [
    'Generation',
    'Rep',
    'MeanRSC',           # Risk of Sperm Competition (lambda for Poisson)
    'MeanMates',         # Mean number of mates per female
    'cor',               # Correlation between male and female traits (related to mating)
    'MeanCount',         # Mean sperm count (related to mating success)
    'MeanMale',          # Male trait (affects mating success)
    'MeanFemale',        # Female trait (affects mate choice)
]

# Check which columns exist
available_columns = [col for col in mating_columns if col in df.columns]
print(f"\nAvailable mating-related columns: {available_columns}")

# Create subset dataframe
df_mating = df[available_columns].copy()

# Add calculated columns if useful
if 'MeanRSC' in df_mating.columns and 'MeanMates' in df_mating.columns:
    df_mating['RSC_Mates_Diff'] = df_mating['MeanRSC'] - df_mating['MeanMates']
    df_mating['RSC_Mates_Ratio'] = df_mating['MeanMates'] / df_mating['MeanRSC']

# Sort by Rep and Generation for easier viewing
df_mating = df_mating.sort_values(['Rep', 'Generation']).reset_index(drop=True)

# Create output filename
output_file = 'parallel_sim_results_2025-11-08_141821_MATING_ONLY.csv'

# Save to CSV
df_mating.to_csv(output_file, index=False)

print(f"\nExtracted CSV shape: {df_mating.shape}")
print(f"Columns in extracted CSV: {list(df_mating.columns)}")
print(f"\nFirst few rows:")
print(df_mating.head(10))
print(f"\nLast few rows:")
print(df_mating.tail(10))

print(f"\n=== Summary Statistics ===")
print(f"Total rows: {len(df_mating)}")
print(f"Total replicates: {df_mating['Rep'].nunique()}")
print(f"Total generations: {df_mating['Generation'].max()}")

if 'MeanRSC' in df_mating.columns:
    print(f"\nMeanRSC - Min: {df_mating['MeanRSC'].min():.3f}, Max: {df_mating['MeanRSC'].max():.3f}, Mean: {df_mating['MeanRSC'].mean():.3f}")

if 'MeanMates' in df_mating.columns:
    print(f"MeanMates - Min: {df_mating['MeanMates'].min():.3f}, Max: {df_mating['MeanMates'].max():.3f}, Mean: {df_mating['MeanMates'].mean():.3f}")

if 'RSC_Mates_Diff' in df_mating.columns:
    print(f"RSC-Mates Difference - Min: {df_mating['RSC_Mates_Diff'].min():.3f}, Max: {df_mating['RSC_Mates_Diff'].max():.3f}, Mean: {df_mating['RSC_Mates_Diff'].mean():.3f}")

print(f"\n✓ Saved to: {output_file}")

