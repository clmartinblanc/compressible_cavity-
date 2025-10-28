/*Code compressible flow solve: 
- Using finite difference method to solve 1D compressible flow equations (Euler equations)
- Flow inside a square sides of length L with left, bottom and right walls are stationary.
- Top wall moves with a sinusoidal velocity profile in time.
- Initial conditions: flow is stationary with inside the box.
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <algorithm>

/* Reference parameters */
const double L = 1.0;        // Length of the square domain
const double T = 1.0;        // Total time of simulation
const int nx = 100;      // Number of grid points in x-direction
const int ny = 100;      // Number of grid points in y-direction
const int nt = 2000;     // Number of time steps
const double gamma = 1.4;   // Specific heat ratio
const double dx = L / nx;   // Grid spacing in x-direction
const double dy = L / ny;   // Grid spacing in y-direction
double dt = T / nt;   // Time step size
const double CFL = 0.1;   // CFL number for stability

const double omega = 2 * M_PI / T; // Angular frequency
const double nu = 1.0e-6;              // Kinematic viscosity
const double mu = 1.8e-5;               // Dynamic viscosity
const double U_w = 1.0;                // Characteristic velocity (max velocity of top wall)
const double L_char = L * sqrt(omega / (2 * nu)); // Characteristic length scale based on frequency
const double Re = U_w * L_char / nu; // Reynolds number
const double Ma = 0.025;                  // Mach number U_w / a (a: speed of sound)
const double Pr = 0.71;                   // Prandtl number rho * cp * mu / k
const double rho0 = 1.225;              // Reference density (kg/m³)
const double cp = 1005.0;               // Specific heat at constant pressure (J/kg·K)
const double lambda = mu * cp / Pr;     // Thermal conductivity
const double k = lambda / (rho0 * cp);  // Thermal diffusivity


/* Initialize flow variables as 2D arrays */
std::vector<std::vector<double>> rho(ny, std::vector<double>(nx, rho0));      // Density 
std::vector<std::vector<double>> u(ny, std::vector<double>(nx, 0.0));         // Velocity in x-direction
std::vector<std::vector<double>> v(ny, std::vector<double>(nx, 0.0));         // Velocity in y-direction
std::vector<std::vector<double>> p(ny, std::vector<double>(nx, 101325.0));    // Pressure (Pa)
std::vector<std::vector<double>> Temp(ny, std::vector<double>(nx, 300.0));    // Temperature 300 K
std::vector<std::vector<double>> E(ny, std::vector<double>(nx));              // Total energy per unit volume
std::vector<std::vector<double>> a(ny, std::vector<double>(nx, 347.2));       // Speed of sound at 300 K

/*Sinusoidal profile at the top wall*/
double topWallVelocity(double t) {
    return U_w * sin(omega * t);
}

/* Flux calculation functions for 2D Navier-Stokes equations */

// Calculate convective flux in x-direction
void calculateConvectiveFluxX(double rho, double u, double v, double p, double E, 
                             double* F_rho, double* F_rhou, double* F_rhov, double* F_E) {
    *F_rho = rho * u;
    *F_rhou = rho * u * u + p;
    *F_rhov = rho * u * v;
    *F_E = (E + p) * u;
}

// Calculate convective flux in y-direction
void calculateConvectiveFluxY(double rho, double u, double v, double p, double E,
                             double* G_rho, double* G_rhou, double* G_rhov, double* G_E) {
    *G_rho = rho * v;
    *G_rhou = rho * u * v;
    *G_rhov = rho * v * v + p;
    *G_E = (E + p) * v;
}

// Calculate viscous flux in x-direction
void calculateViscousFluxX(double rho, double u, double v, double T,
                          double dudx, double dudy, double dvdx, double dvdy, double dTdx,
                          double* Fv_rho, double* Fv_rhou, double* Fv_rhov, double* Fv_E) {
    double mu_eff = mu; // Could be temperature-dependent
    double k_eff = lambda;
    
    // Viscous stress tensor components
    double tau_xx = (2.0/3.0) * mu_eff * (2.0 * dudx - dvdy);
    double tau_xy = mu_eff * (dudy + dvdx);
    
    *Fv_rho = 0.0;
    *Fv_rhou = tau_xx;
    *Fv_rhov = tau_xy;
    *Fv_E = u * tau_xx + v * tau_xy + k_eff * dTdx;
}

// Calculate viscous flux in y-direction
void calculateViscousFluxY(double rho, double u, double v, double T,
                          double dudx, double dudy, double dvdx, double dvdy, double dTdy,
                          double* Gv_rho, double* Gv_rhou, double* Gv_rhov, double* Gv_E) {
    double mu_eff = mu;
    double k_eff = lambda;
    
    // Viscous stress tensor components
    double tau_yy = (2.0/3.0) * mu_eff * (2.0 * dvdy - dudx);
    double tau_yx = mu_eff * (dudy + dvdx);
    
    *Gv_rho = 0.0;
    *Gv_rhou = tau_yx;
    *Gv_rhov = tau_yy;
    *Gv_E = u * tau_yx + v * tau_yy + k_eff * dTdy;
}

// Equation of state for ideal gas: p = ρRT
double calculatePressure(double rho, double T) {
    const double R = 287.0; // Specific gas constant for air (J/kg·K)
    return rho * R * T;
}

// Calculate total energy per unit volume: E = ρ(e + 1/2(u² + v²))
double calculateTotalEnergy(double rho, double u, double v, double T) {
    double e = cp * T / gamma; // Internal energy per unit mass
    return rho * (e + 0.5 * (u * u + v * v));
}

/* Visualization and Data Export Functions */

// Export flow field data to CSV format
void exportToCSV(const std::vector<std::vector<double>>& rho,
                 const std::vector<std::vector<double>>& u,
                 const std::vector<std::vector<double>>& v,
                 const std::vector<std::vector<double>>& p,
                 const std::vector<std::vector<double>>& Temp,
                 const std::vector<std::vector<double>>& E,
                 int timestep) {
    
    std::string filename = "flow_field_" + std::to_string(timestep) + ".csv";
    std::ofstream file(filename);
    
    // Header
    file << "x,y,rho,u,v,p,T,E,velocity_magnitude,vorticity\n";
    
    // Data points
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            double x = i * dx;
            double y = j * dy;
            double vel_mag = sqrt(u[j][i]*u[j][i] + v[j][i]*v[j][i]);
            
            // Calculate vorticity (∂v/∂x - ∂u/∂y)
            double vorticity = 0.0;
            if (i > 0 && i < nx-1 && j > 0 && j < ny-1) {
                double dvdx = (v[j][i+1] - v[j][i-1]) / (2.0 * dx);
                double dudy = (u[j+1][i] - u[j-1][i]) / (2.0 * dy);
                vorticity = dvdx - dudy;
            }
            
            file << x << "," << y << "," << rho[j][i] << "," << u[j][i] << "," 
                 << v[j][i] << "," << p[j][i] << "," << Temp[j][i] << "," 
                 << E[j][i] << "," << vel_mag << "," << vorticity << "\n";
        }
    }
    
    file.close();
    std::cout << "Exported data to " << filename << std::endl;
}

// Export to VTK format for ParaView visualization
void exportToVTK(const std::vector<std::vector<double>>& rho,
                 const std::vector<std::vector<double>>& u,
                 const std::vector<std::vector<double>>& v,
                 const std::vector<std::vector<double>>& p,
                 const std::vector<std::vector<double>>& Temp,
                 int timestep) {
    
    std::string filename = "flow_field_" + std::to_string(timestep) + ".vtk";
    std::ofstream file(filename);
    
    // VTK Header
    file << "# vtk DataFile Version 3.0\n";
    file << "2D CFD Flow Field Data\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_GRID\n";
    file << "DIMENSIONS " << nx << " " << ny << " 1\n";
    file << "POINTS " << nx*ny << " float\n";
    
    // Grid points
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            file << i*dx << " " << j*dy << " 0.0\n";
        }
    }
    
    // Point data
    file << "POINT_DATA " << nx*ny << "\n";
    
    // Velocity vector field
    file << "VECTORS velocity float\n";
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            file << u[j][i] << " " << v[j][i] << " 0.0\n";
        }
    }
    
    // Scalar fields
    file << "SCALARS pressure float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            file << p[j][i] << "\n";
        }
    }
    
    file << "SCALARS density float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            file << rho[j][i] << "\n";
        }
    }
    
    file << "SCALARS temperature float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            file << Temp[j][i] << "\n";
        }
    }
    
    file.close();
    std::cout << "Exported VTK data to " << filename << std::endl;
}

int main() {
    std::cout << "Starting 2D CFD simulation..." << std::endl;
    
    // Initialize total energy field
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            E[j][i] = calculateTotalEnergy(rho[j][i], u[j][i], v[j][i], Temp[j][i]);
        }
    }
    
    // Calculate stable time step based on CFL condition
    double max_speed = sqrt(gamma * calculatePressure(rho0, 300.0) / rho0); // Speed of sound
    double dt_cfl = CFL * std::min(dx, dy) / max_speed;
    dt = std::min(dt, dt_cfl);
    std::cout << "Using time step: " << dt << " seconds" << std::endl;
    
    /* Time-stepping loop */
    for (int n = 0; n < nt; n++) {
        double t = n * dt; // Current time
        double u_top = topWallVelocity(t); // Velocity at the top wall
        
        // Create temporary arrays for new values
        std::vector<std::vector<double>> rho_new = rho;
        std::vector<std::vector<double>> u_new = u;
        std::vector<std::vector<double>> v_new = v;
        std::vector<std::vector<double>> p_new = p;
        std::vector<std::vector<double>> Temp_new = Temp;
        std::vector<std::vector<double>> E_new = E;
        
        // Update interior points using flux-based Navier-Stokes equations
        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                
                // Calculate gradients needed for viscous fluxes
                double dudx = (u[j][i+1] - u[j][i-1]) / (2.0 * dx);
                double dudy = (u[j+1][i] - u[j-1][i]) / (2.0 * dy);
                double dvdx = (v[j][i+1] - v[j][i-1]) / (2.0 * dx);
                double dvdy = (v[j+1][i] - v[j-1][i]) / (2.0 * dy);
                double dTdx = (Temp[j][i+1] - Temp[j][i-1]) / (2.0 * dx);
                double dTdy = (Temp[j+1][i] - Temp[j-1][i]) / (2.0 * dy);
                
                // Calculate fluxes at cell faces (i+1/2, i-1/2, j+1/2, j-1/2)
                
                // Convective fluxes at i+1/2 face (average states)
                double rho_ip = 0.5 * (rho[j][i] + rho[j][i+1]);
                double u_ip = 0.5 * (u[j][i] + u[j][i+1]);
                double v_ip = 0.5 * (v[j][i] + v[j][i+1]);
                double p_ip = 0.5 * (p[j][i] + p[j][i+1]);
                double E_ip = 0.5 * (E[j][i] + E[j][i+1]);
                
                double F_rho_ip, F_rhou_ip, F_rhov_ip, F_E_ip;
                calculateConvectiveFluxX(rho_ip, u_ip, v_ip, p_ip, E_ip, 
                                       &F_rho_ip, &F_rhou_ip, &F_rhov_ip, &F_E_ip);
                
                // Convective fluxes at i-1/2 face
                double rho_im = 0.5 * (rho[j][i-1] + rho[j][i]);
                double u_im = 0.5 * (u[j][i-1] + u[j][i]);
                double v_im = 0.5 * (v[j][i-1] + v[j][i]);
                double p_im = 0.5 * (p[j][i-1] + p[j][i]);
                double E_im = 0.5 * (E[j][i-1] + E[j][i]);
                
                double F_rho_im, F_rhou_im, F_rhov_im, F_E_im;
                calculateConvectiveFluxX(rho_im, u_im, v_im, p_im, E_im,
                                       &F_rho_im, &F_rhou_im, &F_rhov_im, &F_E_im);
                
                // Convective fluxes at j+1/2 face
                double rho_jp = 0.5 * (rho[j][i] + rho[j+1][i]);
                double u_jp = 0.5 * (u[j][i] + u[j+1][i]);
                double v_jp = 0.5 * (v[j][i] + v[j+1][i]);
                double p_jp = 0.5 * (p[j][i] + p[j+1][i]);
                double E_jp = 0.5 * (E[j][i] + E[j+1][i]);
                
                double G_rho_jp, G_rhou_jp, G_rhov_jp, G_E_jp;
                calculateConvectiveFluxY(rho_jp, u_jp, v_jp, p_jp, E_jp,
                                       &G_rho_jp, &G_rhou_jp, &G_rhov_jp, &G_E_jp);
                
                // Convective fluxes at j-1/2 face
                double rho_jm = 0.5 * (rho[j-1][i] + rho[j][i]);
                double u_jm = 0.5 * (u[j-1][i] + u[j][i]);
                double v_jm = 0.5 * (v[j-1][i] + v[j][i]);
                double p_jm = 0.5 * (p[j-1][i] + p[j][i]);
                double E_jm = 0.5 * (E[j-1][i] + E[j][i]);
                
                double G_rho_jm, G_rhou_jm, G_rhov_jm, G_E_jm;
                calculateConvectiveFluxY(rho_jm, u_jm, v_jm, p_jm, E_jm,
                                       &G_rho_jm, &G_rhou_jm, &G_rhov_jm, &G_E_jm);
                
                // Viscous fluxes
                double Fv_rho, Fv_rhou, Fv_rhov, Fv_E;
                calculateViscousFluxX(rho[j][i], u[j][i], v[j][i], Temp[j][i],
                                    dudx, dudy, dvdx, dvdy, dTdx,
                                    &Fv_rho, &Fv_rhou, &Fv_rhov, &Fv_E);
                
                double Gv_rho, Gv_rhou, Gv_rhov, Gv_E;
                calculateViscousFluxY(rho[j][i], u[j][i], v[j][i], Temp[j][i],
                                    dudx, dudy, dvdx, dvdy, dTdy,
                                    &Gv_rho, &Gv_rhou, &Gv_rhov, &Gv_E);
                
                // Update conservative variables using finite volume method
                // ∂U/∂t + ∂F/∂x + ∂G/∂y = ∂Fv/∂x + ∂Gv/∂y
                
                rho_new[j][i] = rho[j][i] - dt * (
                    (F_rho_ip - F_rho_im) / dx + (G_rho_jp - G_rho_jm) / dy
                );
                
                // Update momentum in x-direction: ∂(ρu)/∂t + ∂F_rhou/∂x + ∂G_rhou/∂y = ∂Fv_rhou/∂x + ∂Gv_rhou/∂y
                double rhou_old = rho[j][i] * u[j][i];
                double rhou_new = rhou_old - dt * (
                    (F_rhou_ip - F_rhou_im) / dx + (G_rhou_jp - G_rhou_jm) / dy -
                    (Fv_rhou * 2.0 / dx + Gv_rhou * 2.0 / dy)  // Simple viscous flux discretization
                );
                
                // Update momentum in y-direction
                double rhov_old = rho[j][i] * v[j][i];
                double rhov_new = rhov_old - dt * (
                    (F_rhov_ip - F_rhov_im) / dx + (G_rhov_jp - G_rhov_jm) / dy -
                    (Fv_rhov * 2.0 / dx + Gv_rhov * 2.0 / dy)
                );
                
                // Update energy
                E_new[j][i] = E[j][i] - dt * (
                    (F_E_ip - F_E_im) / dx + (G_E_jp - G_E_jm) / dy -
                    (Fv_E * 2.0 / dx + Gv_E * 2.0 / dy)
                );
                
                // Calculate new velocities
                u_new[j][i] = rhou_new / rho_new[j][i];
                v_new[j][i] = rhov_new / rho_new[j][i];
                
                // Update temperature and pressure using equation of state
                double e_new = E_new[j][i] / rho_new[j][i] - 0.5 * (u_new[j][i] * u_new[j][i] + v_new[j][i] * v_new[j][i]);
                Temp_new[j][i] = e_new * gamma / cp;
                
                // SIMPLE-based pressure-velocity coupling for realistic pressure field
                
                // Velocity divergence (continuity violation)
                double div_u = (u[j][i+1] - u[j][i-1]) / (2.0 * dx) + (v[j+1][i] - v[j-1][i]) / (2.0 * dy);
                
                // Pressure gradients from momentum balance
                double dp_dx = (p[j][i+1] - p[j][i-1]) / (2.0 * dx);
                double dp_dy = (p[j+1][i] - p[j-1][i]) / (2.0 * dy);
                
                // Dynamic pressure from Bernoulli equation
                double vel_mag_sq = u_new[j][i]*u_new[j][i] + v_new[j][i]*v_new[j][i];
                double dynamic_pressure = 0.5 * rho_new[j][i] * vel_mag_sq;
                
                // Pressure Poisson source term (strong momentum-pressure coupling)
                double momentum_source = -rho_new[j][i] * div_u / dt;
                
                // Pressure correction based on neighboring cells (Poisson solver step)
                double neighbor_avg = 0.25 * (p[j][i+1] + p[j][i-1] + p[j+1][i] + p[j-1][i]);
                double poisson_correction = neighbor_avg + momentum_source * dx * dx * 0.25;
                
                // Combine effects with proper weighting
                double base_pressure = calculatePressure(rho_new[j][i], Temp_new[j][i]);
                double pressure_correction = 0.7 * poisson_correction + 0.3 * base_pressure;
                
                // Add dynamic pressure effects
                p_new[j][i] = pressure_correction + dynamic_pressure * 0.1;
            }
        }
        
        // Apply boundary conditions
        for (int i = 0; i < nx; i++) {
            // Bottom wall (j = 0) - stationary, no-slip
            u_new[0][i] = 0.0;
            v_new[0][i] = 0.0;
            // Adiabatic wall: zero temperature gradient
            Temp_new[0][i] = Temp_new[1][i];
            // Calculate pressure from equation of state at boundary
            p_new[0][i] = calculatePressure(rho_new[0][i], Temp_new[0][i]);
            rho_new[0][i] = rho_new[1][i]; // Density from adjacent cell
            E_new[0][i] = calculateTotalEnergy(rho_new[0][i], u_new[0][i], v_new[0][i], Temp_new[0][i]);
            
            // Top wall (j = ny-1) - moving with sinusoidal profile, no-slip in y
            u_new[ny-1][i] = u_top;
            v_new[ny-1][i] = 0.0;
            Temp_new[ny-1][i] = Temp_new[ny-2][i];  // Adiabatic  
            // Calculate pressure from equation of state at moving wall
            p_new[ny-1][i] = calculatePressure(rho_new[ny-1][i], Temp_new[ny-1][i]);
            rho_new[ny-1][i] = rho_new[ny-2][i]; // Density from adjacent cell
            E_new[ny-1][i] = calculateTotalEnergy(rho_new[ny-1][i], u_new[ny-1][i], v_new[ny-1][i], Temp_new[ny-1][i]);
        }
        
        for (int j = 0; j < ny; j++) {
            // Left wall (i = 0) - stationary, no-slip
            u_new[j][0] = 0.0;
            v_new[j][0] = 0.0;
            Temp_new[j][0] = Temp_new[j][1];  // Adiabatic
            // Calculate pressure from equation of state at left wall
            p_new[j][0] = calculatePressure(rho_new[j][0], Temp_new[j][0]);
            rho_new[j][0] = rho_new[j][1]; // Density from adjacent cell
            E_new[j][0] = calculateTotalEnergy(rho_new[j][0], u_new[j][0], v_new[j][0], Temp_new[j][0]);
            
            // Right wall (i = nx-1) - stationary, no-slip
            u_new[j][nx-1] = 0.0;
            v_new[j][nx-1] = 0.0;
            Temp_new[j][nx-1] = Temp_new[j][nx-2];  // Adiabatic
            // Calculate pressure from equation of state at right wall  
            p_new[j][nx-1] = calculatePressure(rho_new[j][nx-1], Temp_new[j][nx-1]);
            rho_new[j][nx-1] = rho_new[j][nx-2]; // Density from adjacent cell
            E_new[j][nx-1] = calculateTotalEnergy(rho_new[j][nx-1], u_new[j][nx-1], v_new[j][nx-1], Temp_new[j][nx-1]);
        }
        
        // Update flow variables for next time step
        rho = rho_new;
        u = u_new;
        v = v_new;
        p = p_new;
        Temp = Temp_new;
        E = E_new;
        
        // Print progress and export data occasionally
        if (n % 200 == 0) {
            std::cout << "Time step " << n << " / " << nt << " completed" << std::endl;
            
            // Export data for visualization every 400 time steps
            if (n % 400 == 0) {
                exportToCSV(rho, u, v, p, Temp, E, n);
                exportToVTK(rho, u, v, p, Temp, n);
            }
        }
    }
    
    std::cout << "Simulation completed!" << std::endl;
    return 0;
}
