# CFD Visualization Guide

This folder includes visualizations for the 2D Navier-Stokes CFD simulation results. The visualization system generates plots showing velocity fields, pressure distributions, temperature maps, vorticity patterns, and time evolution.

## Generated Visualization Files

### Data Files (Automatic Export from .cpp)
- **CSV Format**: `flow_field_XXXX.csv` - Complete flow field data with computed quantities
- **VTK Format**: `flow_field_XXXX.vtk` - 3D visualization format

### Plot Types Generated

#### 1. **Quick Visualization**
- **File**: `quick_cfd_visualization.png`
- **Content**: 2×2 panel showing velocity magnitude, pressure, temperature, and vorticity with streamlines

#### 2. **Individual Contour Plots**
- **Velocity Magnitude**: `velocity_magnitude_tXXXX.png`
- **Pressure**: `pressure_tXXXX.png` 
- **Temperature**: `temperature_tXXXX.png`
- **Density**: `density_tXXXX.png`
- **Features**: High-resolution contours with colorbar, domain boundaries, and proper scaling

#### 3. **Vector Field Plots**
- **File**: `velocity_vectors_tXXXX.png`
- **Content**: Velocity magnitude contours with overlaid velocity vector field
- **Features**: Subsampled arrows for clarity, magnitude-based coloring

#### 4. **Streamline Plots**  
- **File**: `streamlines_tXXXX.png`
- **Content**: Flow streamlines over velocity magnitude background
- **Features**: Shows flow patterns and circulation zones

#### 5. **Vorticity Analysis**
- **File**: `vorticity_tXXXX.png`
- **Content**: Vorticity contours with zero-vorticity lines
- **Features**: Red/blue diverging colormap for positive/negative rotation

#### 6. **Comprehensive Analysis**
- **File**: `comprehensive_tXXXX.png`
- **Content**: 2×2 panel with velocity, pressure, temperature, and vorticity
- **Features**: Complete flow state at each time step

#### 7. **Time Series Analysis**
- **File**: `time_series_analysis.png`
- **Content**: Evolution of peak velocity, mean pressure, max vorticity, and kinetic energy
- **Features**: Shows temporal dynamics and flow development

## Visualization Scripts

### 1. Quick Visualization (`quick_visualize.py`)
```bash
python3 quick_visualize.py
```
- **Purpose**: Immediate visualization of latest results
- **Output**: Single 4-panel plot + flow statistics
- **Runtime**: ~5 seconds

### 2. Comprehensive Visualization (`visualize_cfd.py`)
```bash
python3 visualize_cfd.py
```
- **Purpose**: Generate complete visualization series
- **Output**: All plot types for all time steps + time series analysis

