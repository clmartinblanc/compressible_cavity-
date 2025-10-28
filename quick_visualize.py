#!/usr/bin/env python3
"""
Quick CFD Visualization - Simple script for immediate results

This script creates a quick visualization of the latest CFD results.
Run this after your CFD simulation to see immediate results.

Requirements: matplotlib, numpy, pandas
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
import os

def quick_plot():
    """Create a quick 4-panel visualization of the latest results"""
    
    # Find the latest data file
    files = glob.glob("flow_field_*.csv")
    if not files:
        print("No CFD data files found! Run the simulation first.")
        return
    
    latest_file = max(files, key=lambda f: int(f.split('_')[-1].split('.')[0]))
    print(f"Visualizing: {latest_file}")
    
    # Load data
    df = pd.read_csv(latest_file)
    
    # Get grid dimensions
    x_unique = sorted(df['x'].unique())
    y_unique = sorted(df['y'].unique())
    nx, ny = len(x_unique), len(y_unique)
    
    # Create meshgrid
    X, Y = np.meshgrid(x_unique, y_unique)
    
    # Reshape variables
    vel_mag = df['velocity_magnitude'].values.reshape(ny, nx)
    pressure = df['p'].values.reshape(ny, nx)
    temperature = df['T'].values.reshape(ny, nx)
    vorticity = df['vorticity'].values.reshape(ny, nx)
    u = df['u'].values.reshape(ny, nx)
    v = df['v'].values.reshape(ny, nx)
    
    # Create 2x2 subplot
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('CFD Results - Latest Timestep', fontsize=16, fontweight='bold')
    
    # Velocity magnitude
    cs1 = axes[0,0].contourf(X, Y, vel_mag, levels=20, cmap='plasma')
    axes[0,0].set_title('Velocity Magnitude (m/s)')
    axes[0,0].set_xlabel('x (m)')
    axes[0,0].set_ylabel('y (m)')
    axes[0,0].set_aspect('equal')
    plt.colorbar(cs1, ax=axes[0,0])
    
    # Pressure
    cs2 = axes[0,1].contourf(X, Y, pressure, levels=20, cmap='viridis')
    axes[0,1].set_title('Pressure (Pa)')
    axes[0,1].set_xlabel('x (m)')
    axes[0,1].set_ylabel('y (m)')
    axes[0,1].set_aspect('equal')
    plt.colorbar(cs2, ax=axes[0,1])
    
    # Temperature
    cs3 = axes[1,0].contourf(X, Y, temperature, levels=20, cmap='coolwarm')
    axes[1,0].set_title('Temperature (K)')
    axes[1,0].set_xlabel('x (m)')
    axes[1,0].set_ylabel('y (m)')
    axes[1,0].set_aspect('equal')
    plt.colorbar(cs3, ax=axes[1,0])
    
    # Vorticity with streamlines (handle NaN/inf values)
    vorticity_clean = np.nan_to_num(vorticity, nan=0.0, posinf=0.0, neginf=0.0)
    vort_max = np.abs(vorticity_clean).max()
    
    if vort_max > 0:
        levels_vort = np.linspace(-vort_max, vort_max, 21)
        cs4 = axes[1,1].contourf(X, Y, vorticity_clean, levels=levels_vort, cmap='RdBu_r')
    else:
        cs4 = axes[1,1].contourf(X, Y, vorticity_clean, levels=20, cmap='RdBu_r')
    
    # Add streamlines
    axes[1,1].streamplot(X, Y, u, v, density=1.5, color='black', linewidth=0.8)
    
    axes[1,1].set_title('Vorticity (1/s) + Streamlines')
    axes[1,1].set_xlabel('x (m)')
    axes[1,1].set_ylabel('y (m)')
    axes[1,1].set_aspect('equal')
    plt.colorbar(cs4, ax=axes[1,1])
    
    # Add boundary walls to all plots
    for ax in axes.flat:
        # Domain boundaries
        ax.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], 'k-', linewidth=2)
        # Moving top wall
        ax.plot([0, 1], [1, 1], 'r-', linewidth=3, label='Moving Wall')
    
    plt.tight_layout()
    
    # Save the plot
    output_file = 'quick_cfd_visualization.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"Quick visualization saved as: {output_file}")
    
    # Print some statistics
    print(f"\n=== Flow Statistics ===")
    print(f"Max velocity: {vel_mag.max():.4f} m/s")
    print(f"Min pressure: {pressure.min():.2f} Pa")
    print(f"Max pressure: {pressure.max():.2f} Pa")
    print(f"Min temperature: {temperature.min():.2f} K")
    print(f"Max temperature: {temperature.max():.2f} K")
    print(f"Max |vorticity|: {np.abs(vorticity_clean).max():.4f} 1/s")

if __name__ == "__main__":
    quick_plot()