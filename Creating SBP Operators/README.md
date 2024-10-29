This folder contains the files to create the SBP operators for the two point clouds described in the paper (Ω1 and Ω2).

Steps to follow to generate SBP opeartors:
1. Select 'main_function_simple.m' for Ω1 and 'main_function_sharingan.m' for Ω2.
2. Set desired parameters in main_function and run it to create all the SBP operator matrices (Q, E, S, H etc...).The parameters listed in the files are the ones presented the paper. Only adjusting n_x, n_y, n_boundary, and n_boundary_supp (for 'main_function_sharingan.m' only) is needed to reproduce results at various densities. 
3. Run differetiation_operator_experiments to run the tests that were shown in the paper.

File descriptions:
1. differentiation_operator_experiments: Contains the various numerical experiments that were performed on the SBP operators
2. main_function_simple: calls matrix_generator to create all the SBP matrices for Ω1
3. main_function_sharingan: calls matrix_generator to create all the SBP matrices for Ω2
4. matrix_generator: contains numerous functions that create the SBP matrices
5. point_generator: function that generates a point cloud in Ω1 or Ω2. 
