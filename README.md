# singlecellagingdynamics

Stability Analysis of Single-Cell Aging Models
This repository contains MATLAB and Python for theoretical and numerical stability analysis of fixed points in models of single-cell aging. The models investigate regulatory mechanisms involving SIR2, HAP4, and heme dynamics.

Repository Structure
1. SIR2-hap4 Folder (2D Model)
	•	Implements the 2D model describing interactions between SIR2 and HAP4.
	•	Includes:
	◦	bistabilitynew.ipynb – generates phase portraits for WT and mutant strains.
	◦	fp_sir2_hap.m – computes fixed points for the 2D system.
	◦	limitnew.m – analyzes limit cycles in the 2D system.

2. SIR2-HAP4-heme Folder (3D Model)
	•	Implements the 3D model describing interactions between SIR2, HAP4, and heme.
	•	Focused on wild-type (WT) strain parameters.
	•	Includes:
	◦	Phase portraits for the 3D model.
	◦	fpsir_hap_heme.m – computes fixed points for the 3D system.
	◦	stabilitytest3d_eig_numerical.m – compares theoretical (Jacobian + eigenvalues) vs numericalstability analysis.

Supporting Files (3D Model Only)
	•	hillfunction.m Defines Hill functions used for parameter modeling in the 3D system.
	•	hungarian.m Implements the Hungarian algorithm for parameter fitting in the 3D system.

How the Analysis Works
	1	Fixed Point Calculation
	◦	2D and 3D models compute equilibrium points for WT and mutant strains.
	2	Stability Analysis
	◦	Theoretical: Using Jacobian matrices and eigenvalue analysis.
	◦	Numerical: Perturbing points around fixed points and checking convergence behavior.
	3	Visualization
	◦	Phase portraits and limit cycles are generated for both 2D and 3D systems to visualize trajectories and stability.

Usage:
	4	For 2D model:
	◦	Navigate to SIR2-hap4 folder.
	◦	Run bistabilitynew.ipynb for phase portraits.
	◦	Use fp_sir2_hap.m for fixed-point computation.
	◦	Use limitnew.m for limit cycle analysis.
	5	For 3D model:
	◦	Navigate to SIR2-HAP4-heme folder.
	◦	Run stabilitytest3d_eig_numerical.m for stability comparison.
	◦	Use fpsir_hap_heme.m to compute fixed points.

Dependencies
	•	MATLAB core functions (few toolboxes required).
	•	hillfunction.m and hungarian.m (supporting files for 3D model only).

