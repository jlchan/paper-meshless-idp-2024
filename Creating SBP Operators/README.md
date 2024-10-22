This folder contains the files to create the SBP operators for the two point clouds described in the paper (Ω1 and Ω2).

Steps to follow:
1. Set your parameters in main_function and run it to create all the SBP operator matrices (Q, E, S, H etc...)
2. Run differetiation_operator_experiments to run the tests that were shown in the paper.

File descriptions:
1. differetiation_operator_experiments: Contains the various numerical experiments that were performed on the SBP operators
2. main_function: simply calls matrix_generator to create all the SBP matrices
3. matrix_generator: contains numerous functions that create the SBP matrices
4. point_generator: funciton that generates a point cloud in Ω1 or Ω2. 
