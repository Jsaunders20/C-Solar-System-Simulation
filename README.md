# Orbital Simulation and Visualization Toolkit

## Introduction
C++ toolkit for simulating celestial dynamics using Runge-Kutta 4th Order (RK4) integration. Models single/multi-body gravitational interactions and visualizes planetary motion.

## Core Components
1. RK4 Numerical Integration
2. Gravitational Force Calculations 
3. Data Logging and Visualization

## Code Files

### 1. `invelfindsave.cpp`
- Calculates single-body orbital parameters
- Implements gravitational force calculations
- Uses RK4 integration
- Outputs to `invelfind.dat`
- Generates Gnuplot visualizations

### 2. `invelfind.cpp`
- Enhanced version of `invelfindsave.cpp`
- Dynamically adjusts orbital parameters
- Matches desired orbital periods
- Outputs radius, period, and accuracy metrics

### 3. `solarsystem.cpp`
- Multi-body gravitational simulation
- Tracks 5 planets by default
- Reads from `mass.dat` and `ssinitial.dat`
- Outputs to `solarsystem.dat`
- Generates postscript plots

## Methodology

### RK4 Integration
Computes future states using intermediate steps for improved accuracy.

### Gravitational Dynamics
Newton's law: $\vec{F} = -\frac{G M m}{r^3} \vec{r}$

### Files
**Input:**
- `mass.dat`: Planet masses
- `ssinitial.dat`: Initial conditions

**Output:**
- `invelfind.dat`: Orbital trajectories
- `solarsystem.dat`: Simulation data
- `invelfind.ps`, `solarsystem.ps`: Visualizations

## Usage

### Compilation
```bash
g++ invelfindsave.cpp -o invelfindsave
g++ invelfind.cpp -o invelfind
g++ solarsystem.cpp -o solarsystem
```

### Running Simulations

#### Orbit Finding
```bash
./invelfindsave A P e target_T tau
```
**Parameters:**
- `A`: Aphelion
- `P`: Perihelion
- `e`: Orbital eccentricity
- `target_T`: Target period
- `tau`: Time step

#### Solar System
```bash
./solarsystem T
```
**Parameter:**
- `T`: Duration (years)

### Visualization
Programs automatically generate Gnuplot plots

## Future Work
- Relativistic effects
- Additional visualization formats
- User-friendly parameter input
