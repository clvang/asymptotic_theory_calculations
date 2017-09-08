# asymptotic_theory_calculations

## Purpose:
- This program calculates the diffusivities by solving the asymptotic theory equation from Aharon and Shaw using Bisection method.

- The program also returns the uncertainty in the diffusivity calculations using the Taylor Series Method.

- The main program is contained in the file "asymTheoryDiffusivity.f90".  The make file for the main program is "makefile"

- To compile envoke "make -f makefile" in the command prompt.

- You can change the experimental data input variables in "fcprops.txt".

- The script "mc_analysis.R" calculates uncertainties in D using the Monte Carlo Method.  It is an R script, and relies on built in R function, as well as the shared object "mc_uncertainty.so".  "Rf90_sharedObjects.make" is the make file for the shared object "mc_uncertainty.so".  This object is required to run the MC analysis in the script "mc_analysis.R".

- "normal.f90" is an external library used to generate random variables for a normal distribution with specified mean and standard variation.
