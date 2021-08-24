# AE_LE_crack
An analytic element for intersecting cracks in a linear elastic half-space. A program which calculated the stress field for given number of cracks, which may intersect, in a uniform stress field. The input files for the plot properties and the physical constants are provided through the MATLAB program `run_ae_master.m` and then calls the program `AE_LE_master.cpp` to solve and plot the system.

The program includes the Cartesian stress field (<img src="https://latex.codecogs.com/gif.latex?\sigma_{11}"/> , <img src="https://latex.codecogs.com/gif.latex?\sigma_{22}"/>  & <img src="https://latex.codecogs.com/gif.latex?\sigma_{12}"/> ), the principal stress field (<img src="https://latex.codecogs.com/gif.latex?\sigma_{1}"/> , <img src="https://latex.codecogs.com/gif.latex?\sigma_{2}"/>  & <img src="https://latex.codecogs.com/gif.latex?\theta_{p}"/> ) and the principal stress trajectories. Also, it computes the displacement filed <img src="https://latex.codecogs.com/gif.latex?w"/> and the displacement trajectories. For a given resolution and coordinates.

This program has been developed using *Microsoft Visual Stuido* and only the `.cpp` file and `.exe` is included in the repository. The solution also uses the Eigen library (Guennebaud & Jacob, 2010).

## Instructions
The plots are generated using the MATLAB program `run_ae_master.m`. The script calls the C++ program which solves the system and plots the results. To run teh program simply run the MATLAB script.

## Input data
This list contains all the definitions of the user input data defined in `run_ae_master.m`. These are the model properties:
- `nu` the Poisson's ratio
- `sigma_11inf` the uniform stress state
- `z1` starting coordinates for cracks
- `z2` end coordinates for cracks
- `p` vector for pressures of the cracks
- `ma` number of coefficients to sovle
- `Na` number of interation points along the cracks

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

## Citations
The program has been used in the following paper:
- Comming soon.

## References
Guennebaud, G., Jacob, B., et al. (2010). Eigen v3. http://eigen.tuxfamily.org.
