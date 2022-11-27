# Selected Code in Numerical Analysis
Author: Tobias Timofeyev

## Content:
### 1. Approximating Bessel Functions of the First Kind:
  - `Bessel_Project_Presentation.pptx` A slideshow that provides context and setup for the project
  -  `Tobias_Euler_2ndOrder.m` Function that can be called to approximate a DE using the euler approximation. 
      Takes as input the function, f, initial position, y0, initial change in position, v0, total x-interval, and step size, h.
  -  `Tobias_cRK_2ndOrder.m` Function that can be called to approximate a DE using the Classical Runge-Kutta (rk4) approximation. 
      Takes as input the function, f, initial position, y0, initial change in position, v0, total x-interval, and step size, h.
  - `Bessel_Approx.m` Main file that approximates a few Bessel functions of the first kind
 
### 2. Selected Numerical Differential Equations Code:
   #### Functions for numerical approximation methods
  -  `Tobias_Euler_2ndOrder.m` Function that can be called to approximate a DE using the euler approximation. 
      Takes as input the function, f, initial position, y0, initial change in position, v0, total x-interval, and step size, h.
  -  `Tobias_cRK_2ndOrder.m` Function that can be called to approximate a DE using the Classical Runge-Kutta (rk4) approximation. 
      Takes as input the function, f, initial position, y0, initial change in position, v0, total x-interval, and step size, h.
  -  `Tobias_Center_2ndOrder.m` Function that can be called to approximate a DE using the Central difference approximation. 
      Takes as input the function, f, initial position, y0, initial change in position, v0, total x-interval, and step size, h.
  -  `Tobias_ModEuler_2ndOrder.m` Function that can be called to approximate a DE using the Modified Euler approximation. 
      Takes as input the function, f, initial position, y0, initial change in position, v0, total x-interval, and step size, h.
  -  `Tobias_ModEuler_Matrix.m` Function that can be called to approximate a DE using the Modified Euler approximation. 
      This differs syntactically from the previous in that it accepts and modifies a matrix instead of distinct vectors. 
      Takes as input the function, f, initial condition (vector), y0, , total x-interval, and step size, h.
  -  `Tobias_SympEuler_2ndOrder.m` Function that can be called to approximate a DE using the Symplectic Euler approximation. 
      Takes as input the function, f, initial position, y0, initial change in position, v0, total x-interval, and step size, h.
  -  `Tobias_Verlet1_2ndOrder.m` Function that can be called to approximate a DE using the Verlet approximation. 
      Takes as input the function, f, initial position, y0, initial change in position, v0, total x-interval, and step size, h.
  -  `Tobias_ModEuler_2ndOrder.m` Function that can be called to approximate a DE using the Modified Euler approximation. 
      Takes as input the function, f, initial position, y0, initial change in position, v0, total x-interval, and step size, h.
      
  #### Code Exemplifying the use of numerical methods as well as issues surrounding their use
  -   `12_Von_Neuman_1DHE_Dominating_Error.m` Graph of numerical scheme for the heat equation where von neuman analysis would tell us that that the error will dominate the solution.
  -   `13_Crank_Nicholson_1DHE_HomogBCs.m` Demonstrates Crank Nicholson Scheme on a 1 dimensional heat equation with homogeneous boundary conditions.
  -   `13_Von_Neumann_AmpFact_Analysis.m` When doing von neuman analysis on a numerical scheme, an amplification factor, rho, is found, which represents the growth of error in the problem. 
       This code was used in the analysis of this amplification factor  for different parameter values in a model problem.
  -   `14_Crank_Nicholson_1DHE_MixedBCs.m` Demonstrates Crank Nicholson Scheme on a 1 dimensional heat equation with mixed boundary conditions.
  -   `14_Crank_Nicholson_1DHE_VariableCoeff.m` Demonstrates Crank Nicholson Scheme on a 1 dimensional heat equation where the model equation has variable coefficients.
  -   `14_Semi_Implicit_NonlinearPDE.m` Demonstrates Semi Implicit scheme for a nonlinear partial differential equation.
  -   `15_Peaceman_Rachforth_2DHE_NonHomogBCs.m` Demonstrates Peaceman Rachforth scheme on a 2 dimensional heat equation with nonhomogeneous boundary conditions.
  -   `15_Peaceman_Rachforth_2DHeat_Equation.m` Demonstrates Peaceman Rachforth scheme on a 2 dimensional heat equation with nonhomogeneous boundary conditions.
  -   `15_Simple_Explicit_2DHE_NonHomogBCs.m` Demonstrates a simple explicit scheme for the 2 dimensional heat equation with nonhomogeneous boundary conditions.
  -   `15_Simple_Explicit_2DHeat_Equation.m` Demonstrates a simple explicit scheme for the 2 dimensional heat equation with homogeneous boundary conditions.
  -   `4_Adams_Bashforth_Stability_Region.m` Plots the stability region of the Adams Bashforth method.
  -   `4_Mod_Euler_Stability_Region.m`  Plots the stability region of the Modified Euler method.
  -   `5_Comparison_HarmOscillator_Friction_Symptlectic.m` Plots 2nd order Simple Euler and Symplectic Euler to compare their performance in approximating long term behavior of a harmonic oscillator with a friction term. 
      The result is that the Symplectic Euler method much more accurately describes the weakening of the oscillatory behavior over time, as we would expect in the corresponding physical system.
  -   `5_Comparison_HarmOscillator_Symplectic.m` Compares Verlet, Modified Euler, Classical Runge-Kutta, and Central Difference methods, along with matlab's built-in ode45 function, in their ability to approximate a simple harmonic oscillator over a long time interval.
  -   `7_BVP_Nonlinear_Shooting.m` Demonstrates a nonlinear boundary value problem solved with the shooting method.
  -   `7_BVP_Shooting_Method.m` Demonstrates a second order boundary value problem solved with the shooting method.
  -   `7_BVP_Third_order_Shooting.m` Demonstrates a third order boundary value problem solved with the shooting method.
  -   `8_BVP_Finite_Element_Method.m` Demonstrates a boundary value problem solved with a finite difference method.
  -   `9_BVP_Finite_Element_Colocation.m` Uses Finite Element Collocation method to solve a boundary value problem.



