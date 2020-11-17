# AE_LE_crack
An analytic element for a crack in a linear elastic half-space. A program which calculated the stress field for given number of cracks in a uniform stress field. The input files for the plot properties and the physical constants are provided through the MATLAB program `bin_write.m` and the program `bin_read.m` can be used to plot the results calculated by `AE_LE_master.cpp`.

At the moment the programs does not sync automatically and must be manually run in both MATLAB and the preferred C++ editor. Also, the function for the stresses are coded in the C++ file and changing these, or appending new functions, requires changes in the main code. 

The program includes the Cartesian stress field (<img src="https://latex.codecogs.com/gif.latex?\sigma_{11}"/> , <img src="https://latex.codecogs.com/gif.latex?\sigma_{22}"/>  & <img src="https://latex.codecogs.com/gif.latex?\sigma_{12}"/> ), the principal stress field (<img src="https://latex.codecogs.com/gif.latex?\sigma_{1}"/> , <img src="https://latex.codecogs.com/gif.latex?\sigma_{2}"/>  & <img src="https://latex.codecogs.com/gif.latex?\theta_{p}"/> ) and the principal stress trajectories. Also, it computes the displacement filed <img src="https://latex.codecogs.com/gif.latex?w"/> and the displacement trajectories. For a given resolution and coordinates.

A log file, `log.txt`, is created with all the variable data and time stamp.

This program has been developed using *Microsoft Visual Stuido* and only the `.cpp` file is included in the repository. The solution also uses the Eigen library (Guennebaud & Jacob, 2010).

## Instructions
To generate plot for stresses filed (here only for gravity at the moment) follow the following procedure:
1. Run `bin_write.m`
   - Generates the input data and plot data in `input_data.bin` and `plot_data.bin`.
2. Run `AE_LE_master.cpp`
   - Solves for the coefficients and calcualtes the stress fields and provide them in `data.bin` and the dimensions in `dim_data.bin`.
3. Run `bin_read.m`
   - Reads the data and plots the stress fields.
   
## Included
The program inludes several analytic element for linear elasticity. Currently the following anlaytic elements are included:
- uniform stress and
- cracks.

## Input data
This list contains all the definitions of the user input data defined in `bin_wirte.m`. These are the model properties:
- `nu` the Poisson's ratio (only used in `bin_wirte.m`)
- `kappa` the bulk modulus (is calculated by default with `kappa = 3-4*nu`)
- `sigma_11inf` the uniform stress state
- `nc` number of cracks (calculated in `bin_wirte.m`)
- `m` number of coefficients to sovle
- `z1` starting coordinates for cracks
- `z2` end coordinates for cracks
- `L` length of cracks (calculated in `bin_wirte.m`)
- `mu` angel of cracks (calculated in `bin_wirte.m`)
- `beta` coefficients for cracks

These are the plotting properties:
- `xfrom` the starting value for x-axis
- `xto` the end value for the x-axis
- `yfrom` the starting value for the y-axis
- `yto` the end value for the y-axis
- `Nx` the number of grid points in the x-diraction
- `Ny` the number of grid points in the y-direction
- `Ntraj` the number of steps for the stress trajectories
- `lvs_traj` the number of stress trajectories
- `xtraj` the vector containing the start and end point for the line hwere the trajectores are evenly spaced, for <img src="https://latex.codecogs.com/gif.latex?\sigma_{1}"/>
- `ytraj` the vector containing the start and end point for the line hwere the trajectores are evenly spaced, for <img src="https://latex.codecogs.com/gif.latex?\sigma_{2}"/>
- `Nw` the number of grid points in the x- and y-direction for the displacements

## References
Guennebaud, G., Jacob, B., et al. (2010). Eigen v3. http://eigen.tuxfamily.org.
