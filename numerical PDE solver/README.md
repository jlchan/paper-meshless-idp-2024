This folder contains the various files for solving different types of non-linear conservation laws using the meshfree SBP operators.
Contains the following files:
1. "advection.jl": Solves advection equation and gives the convergence rates
2. "burgers.jl": Solves burgers equation and gives the convergence rates
3. "euler.jl": Solves the 2D incompressible euler equations and contains functions for visualization/plotting
4. "euler_convergence.jl": Specifically outputs the convergence rates using the meshfree SBP operators
5. "operator_setup.jl": function that converts the .mat files provided 
