import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the latest CFD data
df = pd.read_csv('flow_field_1600.csv')

# Basic pressure statistics
print("=== PRESSURE ANALYSIS ===")
print(f"Min pressure: {df['p'].min():.2f} Pa")
print(f"Max pressure: {df['p'].max():.2f} Pa") 
print(f"Mean pressure: {df['p'].mean():.2f} Pa")
print(f"Std deviation: {df['p'].std():.2f} Pa")
print(f"Pressure range: {df['p'].max() - df['p'].min():.2f} Pa")

# Check for pressure anomalies
print("\n=== PRESSURE ANOMALIES ===")
# Check corner values
corners = [
    (0.0, 0.0),    # bottom-left
    (0.99, 0.0),   # bottom-right  
    (0.0, 0.99),   # top-left
    (0.99, 0.99)   # top-right
]

for x_corner, y_corner in corners:
    corner_data = df[(np.abs(df['x'] - x_corner) < 0.005) & (np.abs(df['y'] - y_corner) < 0.005)]
    if not corner_data.empty:
        p_corner = corner_data['p'].iloc[0]
        print(f"Corner ({x_corner}, {y_corner}): p = {p_corner:.2f} Pa")

# Look for pressure gradients near boundaries
print("\n=== BOUNDARY PRESSURE GRADIENTS ===")
# Bottom wall pressure variation
bottom_wall = df[df['y'] == 0.0]
print(f"Bottom wall pressure range: {bottom_wall['p'].min():.2f} - {bottom_wall['p'].max():.2f} Pa")

# Top wall pressure variation  
top_wall = df[df['y'] == 0.99]
print(f"Top wall pressure range: {top_wall['p'].min():.2f} - {top_wall['p'].max():.2f} Pa")

# Left wall pressure variation
left_wall = df[df['x'] == 0.0]
print(f"Left wall pressure range: {left_wall['p'].min():.2f} - {left_wall['p'].max():.2f} Pa")

# Right wall pressure variation
right_wall = df[df['x'] == 0.99]
print(f"Right wall pressure range: {right_wall['p'].min():.2f} - {right_wall['p'].max():.2f} Pa")

# Check for pressure variations
print("\n=== PRESSURE VARIATIONS ===")
pressure_range = df['p'].max() - df['p'].min()
reference_pressure = df['p'].mean()
relative_variation = pressure_range / reference_pressure * 100

print(f"Relative pressure variation: {relative_variation:.4f}%")
if relative_variation > 1.0:
    print("WARNING: Large pressure variations detected (may indicate instability)")
elif relative_variation < 0.001:
    print("WARNING: Pressure field too uniform (weird)")
else:
    print("Pressure variation reasonable")

# Check pressure-velocity coupling issues
print("\n=== PRESSURE-VELOCITY COUPLING ===")
high_vel_data = df[df['velocity_magnitude'] > 0.01]
if not high_vel_data.empty:
    print(f"Pressure in high velocity regions: {high_vel_data['p'].mean():.2f} ± {high_vel_data['p'].std():.2f} Pa")

low_vel_data = df[df['velocity_magnitude'] < 0.001]
if not low_vel_data.empty:
    print(f"Pressure in low velocity regions: {low_vel_data['p'].mean():.2f} ± {low_vel_data['p'].std():.2f} Pa")