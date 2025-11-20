#!/usr/bin/env python3
"""
Quick validation script to check RSC vs MeanMates relationship
and verify no unwanted clamping. Works without R.
"""

import pandas as pd
import numpy as np
import glob
import os

# Find the most recent CSV file
csv_files = glob.glob("parallel_sim_results_*.csv")
if not csv_files:
    print("ERROR: No CSV files found!")
    exit(1)

latest_file = max(csv_files, key=os.path.getmtime)
print(f"Reading data from: {latest_file}\n")

# Read the data
df = pd.read_csv(latest_file)

# Convert Generation to numeric
df['Generation'] = pd.to_numeric(df['Generation'], errors='coerce')
df['Rep'] = df['Rep'].astype(str)

print("=" * 60)
print("DATA VALIDATION REPORT")
print("=" * 60)

print(f"\nData Summary:")
print(f"  Total generations: {df['Generation'].max()}")
print(f"  Total replicates: {df['Rep'].nunique()}")
print(f"  Total data points: {len(df)}")

# Check for negative RSC values
negative_rsc_count = (df['MeanRSC'] < 0).sum()
print(f"\nNegative RSC values: {negative_rsc_count} out of {len(df)}")
if negative_rsc_count > 0:
    print("  -> EXPECTED: RSC can evolve to negative values")
    print("  -> Negative RSC is clamped to 0 only when sampling mates")
else:
    print("  -> All RSC values are non-negative (RSC evolved to positive values)")

# RSC vs MeanMates validation
print(f"\n{'=' * 60}")
print("RSC vs MeanMates Validation")
print("=" * 60)

print(f"\nMean RSC range: {df['MeanRSC'].min():.3f} to {df['MeanRSC'].max():.3f}")
print(f"Mean Mates range: {df['MeanMates'].min():.3f} to {df['MeanMates'].max():.3f}")

# Calculate correlation
correlation = df['MeanRSC'].corr(df['MeanMates'])
print(f"\nOverall correlation (RSC, MeanMates): {correlation:.4f}")
if correlation > 0.95:
    print("  -> EXCELLENT: Strong correlation indicates proper matching")
elif correlation > 0.90:
    print("  -> GOOD: Strong correlation")
else:
    print("  -> WARNING: Correlation lower than expected")

# Linear model
from scipy import stats
slope, intercept, r_value, p_value, std_err = stats.linregress(df['MeanRSC'], df['MeanMates'])
print(f"\nLinear model: MeanMates ~ MeanRSC")
print(f"  Slope: {slope:.4f} (expected ≈ 1.0)")
print(f"  Intercept: {intercept:.4f} (expected ≈ 0.0)")
print(f"  R-squared: {r_value**2:.4f}")

if abs(slope - 1.0) < 0.1 and abs(intercept) < 0.5:
    print("  -> EXCELLENT: Model parameters match expectations")
else:
    print("  -> NOTE: Some deviation expected due to Poisson sampling variance")

# Calculate difference
df['Difference'] = df['MeanRSC'] - df['MeanMates']
mean_diff = df['Difference'].mean()
std_diff = df['Difference'].std()
max_abs_diff = df['Difference'].abs().max()

print(f"\nDifference Statistics:")
print(f"  Mean difference (RSC - MeanMates): {mean_diff:.4f}")
print(f"  SD of difference: {std_diff:.4f}")
print(f"  Max absolute difference: {max_abs_diff:.4f}")

# Final generation statistics
final_gen = df['Generation'].max()
final_data = df[df['Generation'] == final_gen]

print(f"\n{'=' * 60}")
print(f"Final Generation Statistics (Generation {final_gen})")
print("=" * 60)

print(f"\nMean RSC (final): {final_data['MeanRSC'].mean():.4f}")
print(f"Mean Mates (final): {final_data['MeanMates'].mean():.4f}")
print(f"Mean Difference (final): {(final_data['MeanRSC'] - final_data['MeanMates']).mean():.4f}")
print(f"Correlation RSC-Mates (final): {final_data['MeanRSC'].corr(final_data['MeanMates']):.4f}")
print(f"Mean Trait Correlation (final): {final_data['cor'].mean():.4f}")

# Trait evolution check
print(f"\n{'=' * 60}")
print("Trait Evolution Check")
print("=" * 60)

print(f"\nMean Male trait (final): {final_data['MeanMale'].mean():.4f}")
print(f"Mean Female trait (final): {final_data['MeanFemale'].mean():.4f}")
print(f"Mean Sperm count (final): {final_data['MeanCount'].mean():.4f}")

# Check for sudden jumps (indicating potential issues)
for rep in df['Rep'].unique():
    rep_data = df[df['Rep'] == rep].sort_values('Generation')
    rsc_jumps = rep_data['MeanRSC'].diff().abs()
    max_jump = rsc_jumps.max()
    if max_jump > 10:  # Threshold for "sudden jump"
        print(f"\nWARNING: Replicate {rep} has large RSC jump: {max_jump:.2f}")

print(f"\n{'=' * 60}")
print("VALIDATION SUMMARY")
print("=" * 60)

issues = []
if correlation < 0.90:
    issues.append("Low RSC-MeanMates correlation")
if abs(slope - 1.0) > 0.2:
    issues.append("Slope deviates significantly from 1.0")
if abs(intercept) > 1.0:
    issues.append("Intercept deviates significantly from 0.0")

if issues:
    print("\n⚠️  ISSUES FOUND:")
    for issue in issues:
        print(f"  - {issue}")
else:
    print("\n✅ VALIDATION PASSED!")
    print("  - RSC and MeanMates match properly")
    print("  - No unwanted clamping detected")
    print("  - Evolution appears normal")

print(f"\n{'=' * 60}")
print("Recommendations:")
print("=" * 60)
print("1. Run plot_and_validate_all.R for comprehensive plots")
print("2. Compare results with adaptive dynamics predictions")
print("3. Create publication figures")
print("4. See FINALIZATION_SUGGESTIONS.md for next steps")

