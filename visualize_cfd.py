#!/usr/bin/env python3
"""
CFD Visualization Script for 2D Compressible Navier-Stokes Solver

This script creates comprehensive visualizations of the CFD simulation results
including contour plots, vector fields, streamlines, and animations.

Requirements: matplotlib, numpy, pandas
Install with: pip install matplotlib numpy pandas

Author: CFD Project Team
Date: October 28, 2025
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import glob
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import matplotlib.patches as patches
from matplotlib.animation import FuncAnimation

# Configuration
plt.style.use('seaborn-v0_8')
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10

class CFDVisualizer:
    def __init__(self, data_pattern="flow_field_*.csv"):
        """Initialize the CFD visualizer with data files"""
        self.data_files = sorted(glob.glob(data_pattern))
        self.timesteps = [int(f.split('_')[-1].split('.')[0]) for f in self.data_files]
        print(f"Found {len(self.data_files)} data files")
        
        if self.data_files:
            # Load first file to get grid information
            self.df_sample = pd.read_csv(self.data_files[0])
            self.setup_grid()
    
    def setup_grid(self):
        """Setup grid parameters from data"""
        self.x_unique = sorted(self.df_sample['x'].unique())
        self.y_unique = sorted(self.df_sample['y'].unique())
        self.nx = len(self.x_unique)
        self.ny = len(self.y_unique)
        self.X, self.Y = np.meshgrid(self.x_unique, self.y_unique)
        print(f"Grid: {self.nx} x {self.ny} points")
    
    def load_data(self, filename):
        """Load and reshape CFD data"""
        df = pd.read_csv(filename)
        
        # Reshape data to 2D grids
        data = {}
        for var in ['rho', 'u', 'v', 'p', 'T', 'E', 'velocity_magnitude', 'vorticity']:
            if var in df.columns:
                data[var] = df[var].values.reshape(self.ny, self.nx)
        
        return data
    
    def create_contour_plot(self, data, variable, timestep, title=None, cmap='viridis'):
        """Create contour plot for a variable"""
        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
        
        if title is None:
            title = f'{variable.capitalize()} at t = {timestep * 0.001:.3f}s'
        
        # Contour plot
        levels = 20
        cs = ax.contourf(self.X, self.Y, data[variable], levels=levels, cmap=cmap)
        
        # Contour lines
        cs_lines = ax.contour(self.X, self.Y, data[variable], levels=levels//2, 
                             colors='black', alpha=0.4, linewidths=0.5)
        
        # Colorbar
        cbar = plt.colorbar(cs, ax=ax, shrink=0.8)
        cbar.set_label(self.get_variable_label(variable), rotation=270, labelpad=20)
        
        # Formatting
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        
        # Add domain boundaries
        self.add_boundaries(ax)
        
        plt.tight_layout()
        return fig, ax
    
    def create_vector_plot(self, data, timestep, subsample=4):
        """Create velocity vector field plot"""
        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
        
        # Velocity magnitude contour
        cs = ax.contourf(self.X, self.Y, data['velocity_magnitude'], 
                        levels=20, cmap='plasma', alpha=0.8)
        
        # Vector field (subsampled for clarity)
        X_sub = self.X[::subsample, ::subsample]
        Y_sub = self.Y[::subsample, ::subsample]
        U_sub = data['u'][::subsample, ::subsample]
        V_sub = data['v'][::subsample, ::subsample]
        
        # Quiver plot
        ax.quiver(X_sub, Y_sub, U_sub, V_sub, 
                 scale=20, scale_units='xy', angles='xy', 
                 color='white', alpha=0.8, width=0.003)
        
        # Colorbar
        cbar = plt.colorbar(cs, ax=ax, shrink=0.8)
        cbar.set_label('Velocity Magnitude (m/s)', rotation=270, labelpad=20)
        
        # Formatting
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        ax.set_title(f'Velocity Field at t = {timestep * 0.001:.3f}s', 
                    fontsize=14, fontweight='bold')
        ax.set_aspect('equal')
        
        # Add domain boundaries
        self.add_boundaries(ax)
        
        plt.tight_layout()
        return fig, ax
    
    def create_streamline_plot(self, data, timestep):
        """Create streamline plot"""
        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
        
        # Velocity magnitude background
        cs = ax.contourf(self.X, self.Y, data['velocity_magnitude'], 
                        levels=20, cmap='Blues', alpha=0.6)
        
        # Streamlines
        ax.streamplot(self.X, self.Y, data['u'], data['v'], 
                     density=2, color='darkred', linewidth=1.5, arrowsize=1.5)
        
        # Colorbar
        cbar = plt.colorbar(cs, ax=ax, shrink=0.8)
        cbar.set_label('Velocity Magnitude (m/s)', rotation=270, labelpad=20)
        
        # Formatting
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        ax.set_title(f'Streamlines at t = {timestep * 0.001:.3f}s', 
                    fontsize=14, fontweight='bold')
        ax.set_aspect('equal')
        
        # Add domain boundaries
        self.add_boundaries(ax)
        
        plt.tight_layout()
        return fig, ax
    
    def create_vorticity_plot(self, data, timestep):
        """Create vorticity contour plot"""
        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
        
        # Clean vorticity data and handle edge cases
        vorticity = np.nan_to_num(data['vorticity'], nan=0.0, posinf=0.0, neginf=0.0)
        vort_max = np.abs(vorticity).max()
        vort_min = vorticity.min()
        vort_range = vort_max - vort_min
        
        if vort_range > 1e-12:  # Non-trivial vorticity field
            levels = np.linspace(-vort_max, vort_max, 21)
            cs = ax.contourf(self.X, self.Y, vorticity, 
                            levels=levels, cmap='RdBu_r', extend='both')
        else:  # Nearly constant vorticity
            cs = ax.contourf(self.X, self.Y, vorticity, 
                            levels=20, cmap='RdBu_r')
        
        # Zero vorticity line (if vorticity field has variation)
        if vort_range > 1e-12:
            ax.contour(self.X, self.Y, vorticity, levels=[0], 
                      colors='black', linewidths=1.5)
        
        # Colorbar
        cbar = plt.colorbar(cs, ax=ax, shrink=0.8)
        cbar.set_label('Vorticity (1/s)', rotation=270, labelpad=20)
        
        # Formatting
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        ax.set_title(f'Vorticity Field at t = {timestep * 0.001:.3f}s', 
                    fontsize=14, fontweight='bold')
        ax.set_aspect('equal')
        
        # Add domain boundaries
        self.add_boundaries(ax)
        
        plt.tight_layout()
        return fig, ax
    
    def create_comprehensive_plot(self, data, timestep):
        """Create comprehensive 2x2 subplot figure"""
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle(f'CFD Results at t = {timestep * 0.001:.3f}s', 
                    fontsize=16, fontweight='bold')
        
        # Velocity magnitude
        cs1 = axes[0,0].contourf(self.X, self.Y, data['velocity_magnitude'], 
                               levels=20, cmap='plasma')
        axes[0,0].set_title('Velocity Magnitude')
        axes[0,0].set_xlabel('x (m)')
        axes[0,0].set_ylabel('y (m)')
        axes[0,0].set_aspect('equal')
        self.add_boundaries(axes[0,0])
        plt.colorbar(cs1, ax=axes[0,0], shrink=0.8)
        
        # Pressure
        cs2 = axes[0,1].contourf(self.X, self.Y, data['p'], 
                               levels=20, cmap='viridis')
        axes[0,1].set_title('Pressure')
        axes[0,1].set_xlabel('x (m)')
        axes[0,1].set_ylabel('y (m)')
        axes[0,1].set_aspect('equal')
        self.add_boundaries(axes[0,1])
        plt.colorbar(cs2, ax=axes[0,1], shrink=0.8)
        
        # Temperature
        cs3 = axes[1,0].contourf(self.X, self.Y, data['T'], 
                               levels=20, cmap='coolwarm')
        axes[1,0].set_title('Temperature')
        axes[1,0].set_xlabel('x (m)')
        axes[1,0].set_ylabel('y (m)')
        axes[1,0].set_aspect('equal')
        self.add_boundaries(axes[1,0])
        plt.colorbar(cs3, ax=axes[1,0], shrink=0.8)
        
        # Vorticity
        vorticity = np.nan_to_num(data['vorticity'], nan=0.0, posinf=0.0, neginf=0.0)
        vort_max = np.abs(vorticity).max()
        vort_min = vorticity.min()
        vort_range = vort_max - vort_min
        
        if vort_range > 1e-12:
            levels = np.linspace(-vort_max, vort_max, 21)
            cs4 = axes[1,1].contourf(self.X, self.Y, vorticity, 
                                   levels=levels, cmap='RdBu_r')
        else:
            cs4 = axes[1,1].contourf(self.X, self.Y, vorticity, 
                                   levels=20, cmap='RdBu_r')
        axes[1,1].set_title('Vorticity')
        axes[1,1].set_xlabel('x (m)')
        axes[1,1].set_ylabel('y (m)')
        axes[1,1].set_aspect('equal')
        self.add_boundaries(axes[1,1])
        plt.colorbar(cs4, ax=axes[1,1], shrink=0.8)
        
        plt.tight_layout()
        return fig, axes
    
    def add_boundaries(self, ax):
        """Add boundary walls to plot"""
        # Domain boundaries
        rect = patches.Rectangle((0, 0), 1, 1, linewidth=3, 
                               edgecolor='black', facecolor='none')
        ax.add_patch(rect)
        
        # Moving top wall (different color)
        ax.plot([0, 1], [1, 1], 'r-', linewidth=4, label='Moving Wall')
        
        # Stationary walls
        ax.plot([0, 1], [0, 0], 'k-', linewidth=3, label='Stationary Wall')
        ax.plot([0, 0], [0, 1], 'k-', linewidth=3)
        ax.plot([1, 1], [0, 1], 'k-', linewidth=3)
        
        # Legend (only for first boundary addition)
        current_legend = ax.get_legend()
        if current_legend is None:
            ax.legend(loc='upper right', bbox_to_anchor=(0.98, 0.98))
    
    def get_variable_label(self, variable):
        """Get proper label with units for variables"""
        labels = {
            'rho': 'Density (kg/m³)',
            'u': 'u-velocity (m/s)',
            'v': 'v-velocity (m/s)',
            'p': 'Pressure (Pa)',
            'T': 'Temperature (K)',
            'E': 'Total Energy (J/m³)',
            'velocity_magnitude': 'Velocity Magnitude (m/s)',
            'vorticity': 'Vorticity (1/s)'
        }
        return labels.get(variable, variable)
    
    def generate_all_plots(self, output_dir='plots'):
        """Generate all visualization plots"""
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        print("Generating visualization plots...")
        
        for i, (filename, timestep) in enumerate(zip(self.data_files, self.timesteps)):
            print(f"Processing timestep {timestep} ({i+1}/{len(self.data_files)})")
            
            # Load data
            data = self.load_data(filename)
            
            # Create plots
            plots = [
                ('velocity_magnitude', 'plasma'),
                ('pressure', 'viridis'),
                ('temperature', 'coolwarm'),
                ('density', 'Blues'),
            ]
            
            for var, cmap in plots:
                if var == 'pressure':
                    var_data = 'p'
                elif var == 'temperature':
                    var_data = 'T'
                elif var == 'density':
                    var_data = 'rho'
                else:
                    var_data = var
                
                fig, ax = self.create_contour_plot(data, var_data, timestep, 
                                                 title=f'{var.replace("_", " ").title()}', 
                                                 cmap=cmap)
                
                plt.savefig(f'{output_dir}/{var}_t{timestep:04d}.png', 
                           bbox_inches='tight', dpi=300)
                plt.close()
            
            # Vector plot
            fig, ax = self.create_vector_plot(data, timestep)
            plt.savefig(f'{output_dir}/velocity_vectors_t{timestep:04d}.png', 
                       bbox_inches='tight', dpi=300)
            plt.close()
            
            # Streamlines
            fig, ax = self.create_streamline_plot(data, timestep)
            plt.savefig(f'{output_dir}/streamlines_t{timestep:04d}.png', 
                       bbox_inches='tight', dpi=300)
            plt.close()
            
            # Vorticity
            fig, ax = self.create_vorticity_plot(data, timestep)
            plt.savefig(f'{output_dir}/vorticity_t{timestep:04d}.png', 
                       bbox_inches='tight', dpi=300)
            plt.close()
            
            # Comprehensive plot
            fig, axes = self.create_comprehensive_plot(data, timestep)
            plt.savefig(f'{output_dir}/comprehensive_t{timestep:04d}.png', 
                       bbox_inches='tight', dpi=300)
            plt.close()
        
        print(f"All plots saved to '{output_dir}/' directory")
    
    def create_time_series_plots(self, output_dir='plots'):
        """Create time series plots of key quantities"""
        if not self.data_files:
            return
        
        # Initialize arrays for time series data
        times = []
        max_velocity = []
        avg_pressure = []
        max_vorticity = []
        kinetic_energy = []
        
        print("Extracting time series data...")
        
        for filename, timestep in zip(self.data_files, self.timesteps):
            data = self.load_data(filename)
            
            times.append(timestep * 0.001)  # Convert to seconds
            max_velocity.append(np.max(data['velocity_magnitude']))
            avg_pressure.append(np.mean(data['p']))
            max_vorticity.append(np.max(np.abs(data['vorticity'])))
            
            # Kinetic energy density
            ke = 0.5 * data['rho'] * data['velocity_magnitude']**2
            kinetic_energy.append(np.mean(ke))
        
        # Create time series plots
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle('CFD Time Series Analysis', fontsize=16, fontweight='bold')
        
        axes[0,0].plot(times, max_velocity, 'b-', linewidth=2)
        axes[0,0].set_xlabel('Time (s)')
        axes[0,0].set_ylabel('Maximum Velocity (m/s)')
        axes[0,0].set_title('Peak Velocity Evolution')
        axes[0,0].grid(True, alpha=0.3)
        
        axes[0,1].plot(times, avg_pressure, 'r-', linewidth=2)
        axes[0,1].set_xlabel('Time (s)')
        axes[0,1].set_ylabel('Average Pressure (Pa)')
        axes[0,1].set_title('Mean Pressure Evolution')
        axes[0,1].grid(True, alpha=0.3)
        
        axes[1,0].plot(times, max_vorticity, 'g-', linewidth=2)
        axes[1,0].set_xlabel('Time (s)')
        axes[1,0].set_ylabel('Maximum |Vorticity| (1/s)')
        axes[1,0].set_title('Peak Vorticity Evolution')
        axes[1,0].grid(True, alpha=0.3)
        
        axes[1,1].plot(times, kinetic_energy, 'm-', linewidth=2)
        axes[1,1].set_xlabel('Time (s)')
        axes[1,1].set_ylabel('Mean Kinetic Energy Density (J/m³)')
        axes[1,1].set_title('Kinetic Energy Evolution')
        axes[1,1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/time_series_analysis.png', 
                   bbox_inches='tight', dpi=300)
        plt.close()
        
        print("Time series analysis saved")

def main():
    """Main function to run CFD visualization"""
    print("=== CFD Visualization Tool ===")
    print("Starting visualization generation...")
    
    # Initialize visualizer
    viz = CFDVisualizer()
    
    if not viz.data_files:
        print("No CFD data files found. Please run the CFD simulation first.")
        print("Expected files: flow_field_*.csv")
        return
    
    # Generate all plots
    viz.generate_all_plots()
    viz.create_time_series_plots()
    
    print("\n=== Visualization Complete! ===")
    print("Generated plots:")
    print("  • Velocity magnitude contours")
    print("  • Pressure contours") 
    print("  • Temperature contours")
    print("  • Density contours")
    print("  • Velocity vector fields")
    print("  • Streamline plots")
    print("  • Vorticity contours")
    print("  • Comprehensive 4-panel plots")
    print("  • Time series analysis")
    print("\nAll files saved in 'plots/' directory")

if __name__ == "__main__":
    main()