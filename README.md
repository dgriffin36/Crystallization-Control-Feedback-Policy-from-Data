# Crystallization Control: Feedback Policy from Data

This code has been posted in association with the manuscript titled "Data-Driven Modeling and Dynamic Programming Applied to Batch Cooling Crystallization" by D. J. Griffin, M. A. Grover, Y. Kawajiri, and R. W. Rousseau. The main function, ObtainPolicy.m, is a MATLAB function that takes as input experimental data on crystallization dynamics and outputs a feedback control policy for reaching a target crystall mass and chord count. To run this function, the m-files in the subfunctions folder must be in the same path. The output policy can be viewed as a time-varying colormap using viewPolicy.m; the output policy can be interpreted using inputFromPolicy.m. 

## ObtainPolicy
The process has two important steps: 1. data-driven modeling, and 2. dynamic programming. In the first step a convex optimization problem is posed. This problem is solved with CVX: Software for Disciplined Convex Programming [1],[2]. For this function to run, CVX must be installed and in the appropriate path. For commercial use with non-free solvers, such as MATLAB, please obtain the appropriate license: (http://cvxr.com/cvx/licensing/).

###--------------------------------------Inputs----------------------------------------
There are a number of inputs. The inputs specify: the training data, the length of the time steps, the target position in mass-count space, the batch run time, adjustable parameters in the optimization formulation, and the space-input discretization to use.
   
REQUIRED

  Xtr     - Contains training data 'positions'.
  
  Utr     - Contains training data inputs.
  
  dXtr    - Contains training data 'movements'. 
  
  dt      - The measurement time interval in minutes.
  
  Dt      - The time interval over which the input is held constant.
  
  xTarget - The target position [count, mass].
  
  N       - The number of time intervals in the control run.

OPTIONAL

  rho     - Parameter that weights the input-effort cost (default: 25).
  
  gamma 	- Parameter that weights the running distance-to-target (default: 5).
  
  Grid    - Structure containing the grid spacing for: count (Grid.c), mass (Grid.m), and supersaturation (Grid.s).
 
###--------------------------------------Note------------------------------------------
Depending on input data and discretization, the function will require substantial computation time (on the order of 30 minutes for the example data). The modeling step takes the longest. Prompts and visuals have been added as progress checks. These should be commented-out to run the function unattended.

###----------------------------------Bibliography-------------------------------------
[1] Michael Grant and Stephen Boyd. CVX: Matlab software for disciplined convex programming, version 2.0 beta. http://cvxr.com/cvx, September 2013.

[2] Michael Grant and Stephen Boyd. Graph implementations for nonsmooth convex programs, Recent Advances in Learning and Control (a tribute to M. Vidyasagar), V. Blondel, S. Boyd, and H. Kimura, editors, pages 95-110, Lecture Notes in Control and Information Sciences, Springer, 2008. http://stanford.edu/~boyd/graph_dcp.html.

###-----------------------------------Copyright---------------------------------------
Copyright 2015, Daniel Griffin.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details: <http://www.gnu.org/licenses/>.

