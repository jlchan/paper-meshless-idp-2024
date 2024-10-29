# paper-meshless-idp-2024
Reproducibility repo for first-order invariant domain-preserving meshless methods

`.mat` files for more resolved operators are also available [here](https://drive.google.com/drive/folders/10aj7Ek_eW2ce-M97Pt8a47SvBFAV6pRl?usp=drive_link).

Each file contains the following set of matrices:
points: Matrix that gives the entire set of points in the domain and classifies it as a boundary point or not
H_min: Optimized H matrix
E_x, E_y: The E matrices used in SBP
Q_x, Q_y: The Q matrices used in SBP
normal_x, normal_y: The normal vectors at each boundary point 

matrices_k_circle refer to Ω1. The first number is n_x (number of evenly spaced points on diameter). The second number is n_boudnary (number of evenly spaced points on boudnary)

matrices_sharingan_k_circle refer to Ω2. The first number is n_x (number of evenly spaced points on diameter). The second number is n_boudnary (number of evenly spaced points on exterior boudnary). The third number is n_int (number of evenly spaced points on each circular interior boudnary).

Only the smaller sets of matrices are uploaded on this GitHub repository due to size constraints. Refer to this googledrive for all relevant sets: https://drive.google.com/drive/folders/10aj7Ek_eW2ce-M97Pt8a47SvBFAV6pRl?usp=sharing

This folder contains the various files for solving different types of non-linear conservation laws using the meshfree SBP operators.
Contains the following files:
1. "advection.jl": Solves advection equation and gives the convergence rates
2. "burgers.jl": Solves burgers equation and gives the convergence rates
3. "euler.jl": Solves the 2D incompressible euler equations and contains functions for visualization/plotting
4. "euler_convergence.jl": Specifically outputs the convergence rates using the meshfree SBP operators
5. "operator_setup.jl": function that converts the .mat files provided 
