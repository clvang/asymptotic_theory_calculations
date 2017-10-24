# asymptotic_theory_calculations

## Purpose:
- This program calculates the diffusivities by solving the asymptotic theory equation from Aharon and Shaw using Bisection method.

- The program reads the experimental paramters from the file "fcprops.txt" and returns the uncertainty in the diffusivity calculations using the Taylor Series Method.

- The main Fortran program is contained in the file "asymTheoryDiffusivity.f90".  The make file for the main program is "makefile"

- To compile envoke "make -f makefile" in the command prompt.

- You can change the experimental data input variables in "fcprops.txt".

- The script "mc_analysis.R" calculates uncertainties in D using the Monte Carlo Method.  Random numbers are generated for the independent variables in the asymptotic theory equation and the roots for D are calculated using the shared Fortran object "mc_uncertainty.so".  The make file for "mc_uncertainty.so" is "Rf90_sharedObjects.make".  To compile the shared object, infolve "make -f Rf90_sharedObjects.make" in the command prompt.
