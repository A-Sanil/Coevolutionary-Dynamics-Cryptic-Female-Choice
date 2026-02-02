#!/usr/bin/env python3
"""
Python script to create visualization images for the three offspring functions
Shows how each function should look: parabola, linear (monotonic), and asymptotic
"""

import math
try:
    import matplotlib.pyplot as plt
except ImportError:
    print("matplotlib not available, trying to install...")
    import subprocess
    import sys
    subprocess.check_call([sys.executable, "-m", "pip", "install", "matplotlib", "--quiet"])
    import matplotlib.pyplot as plt

# Define the functions matching the Julia code
def offspring_parabola(score):
    peak = 4.0
    opt = 0.65
    width = 0.70
    scale = 1 - ((score - opt) / width)**2
    base = peak * max(0.0, scale)
    return max(1.0, base)

def offspring_asymptote(score):
    max_children = 4.0
    k = 4.2
    return max_children * (1 - math.exp(-k * score))

def offspring_linear(score):
    min_offspring = 1.0
    max_offspring = 4.0
    return min_offspring + (max_offspring - min_offspring) * score

# Create data for plotting
scores = [i / 1000.0 for i in range(1001)]
parabola_vals = [offspring_parabola(s) for s in scores]
asymptote_vals = [offspring_asymptote(s) for s in scores]
linear_vals = [offspring_linear(s) for s in scores]

# Create combined plot
plt.figure(figsize=(12, 8))
plt.plot(scores, parabola_vals, label='Parabola', color='#E69F00', linewidth=2)
plt.plot(scores, asymptote_vals, label='Asymptotic', color='#009E73', linewidth=2, linestyle='--')
plt.plot(scores, linear_vals, label='Linear (Monotonic)', color='#0072B2', linewidth=2, linestyle='-.')
plt.xlabel('Compatibility Score', fontsize=14)
plt.ylabel('Expected Offspring Count', fontsize=14)
plt.title('Offspring Functions: Expected Shapes', fontsize=16, fontweight='bold')
plt.legend(fontsize=12)
plt.xlim(0, 1)
plt.ylim(0, 5)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('function_shapes_visualization.png', dpi=300, bbox_inches='tight')
print("Saved: function_shapes_visualization.png")
plt.close()

# Create individual plots
# Parabola
plt.figure(figsize=(8, 6))
plt.plot(scores, parabola_vals, color='#E69F00', linewidth=2)
plt.xlabel('Compatibility Score', fontsize=14)
plt.ylabel('Expected Offspring Count', fontsize=14)
plt.title('Parabola Function', fontsize=16, fontweight='bold')
plt.xlim(0, 1)
plt.ylim(0, 5)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('function_parabola.png', dpi=300, bbox_inches='tight')
print("Saved: function_parabola.png")
plt.close()

# Asymptotic
plt.figure(figsize=(8, 6))
plt.plot(scores, asymptote_vals, color='#009E73', linewidth=2, linestyle='--')
plt.xlabel('Compatibility Score', fontsize=14)
plt.ylabel('Expected Offspring Count', fontsize=14)
plt.title('Asymptotic Function', fontsize=16, fontweight='bold')
plt.xlim(0, 1)
plt.ylim(0, 5)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('function_asymptotic.png', dpi=300, bbox_inches='tight')
print("Saved: function_asymptotic.png")
plt.close()

# Linear (Monotonic)
plt.figure(figsize=(8, 6))
plt.plot(scores, linear_vals, color='#0072B2', linewidth=2, linestyle='-.')
plt.xlabel('Compatibility Score', fontsize=14)
plt.ylabel('Expected Offspring Count', fontsize=14)
plt.title('Linear (Monotonically Increasing) Function', fontsize=16, fontweight='bold')
plt.xlim(0, 1)
plt.ylim(0, 5)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('function_linear.png', dpi=300, bbox_inches='tight')
print("Saved: function_linear.png")
plt.close()

print("All function visualization images created successfully!")

