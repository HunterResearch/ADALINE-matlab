# ADALINE-matlab
Code for ADALINE. See https://link.springer.com/article/10.1007/s10107-021-01667-6

To run the code for ADALINE, please go to servo.m first. 
There are three potential oracles to choose from: Bus scheduling, discrete quadratic, and dynamic newsvendor.
There are two algorithms to choose from, ADALINE and R-SPLINE. 

The servo_parfor.m code creates a parallel pool to do the final runs for the paper. This code creates datafiles that are then called by PlotQuantilesOracleName.m files. These datafiles initially live in the main folder; to keep the files from overwriting each other, manually copy them into an appropriately named subfolder after the runs are complete. So, the procedure is to run servo_parfor.m to get .mat files and then run PlotQuantilesOracleName.m.

Credit to Eric Applegate for creating the initial template that generates the plots.
