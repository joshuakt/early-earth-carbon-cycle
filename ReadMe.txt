Version 1.0

This set of python scripts runs our geological carbon cycle model for the last 4.0 Ga and plots selected outputs alongside proxy data from the literature. The model is described in  J. Krissansen-Totton, G. Arney, and D. C. Catling (2018) "Constraining the climate and ocean pH of the early Earth with a geological carbon cycle model", Proceedings of the National Academy of Sciences. With the appropriate choice of parameters, this code can reproduce Fig. 3, 4, and 5 in the main text, in addition to many supplementary figures.

As a matter of courtesy, we request that people using this code please cite Krissansen-Totton et al. (2018). In the interest of an "open source" approach, we also request that authors who use and modify the code, please send a copy of papers and modified code to the lead author (joshkt@uw.edu)

REQUIREMENTS: Python, including numpy, pylab, and scipy modules.

HOW TO RUN CODE:
(1) Put all the python scripts in the same directory, and ensure python is working in this directory.
(2) Open Main_code_parallel.py and check desired parameter ranges, options, and number of iterations (default parameters reproduce Fig. 3 in main text)
(3) Run Main_code.py. Code will output model confidence interval alongside proxy data. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EXPLANATION OF CODE STRUCTURE:

%% Main_code_parallel.py
This script provides the shell to repeatedly call the forward model and plot the output distributions. The "Options" section of the code contains various parameters that can be modified by the user:

- it_num: number of forward model calls used to build output distributions. The larger this is, the longer the code will take to run.
- Parallelize: determines whether multiple threads should be used (faster) or whether all iterations should be computer using the same thread (slower).
- Carbon_chem: determines whether carbon equilibrium chemistry constants should be temperature dependent (slower but more accurate) or fixed (faster but less accurate).
- methane: controls whether methane should be added to the Proterozoic and Archean atmospheres

Ranges for uncertain parameters can also be modified by the user. Modifying anything else in 'Main_code_parallel.py' may produce errors in the code. Once parameter ranges have been defined, the script calls the forward model ('Forward_Model') once and uses the dimensions of the outputs to define an output array to store all future outputs ('all_output'). The forward model is then called many times (equal to 'it_num') and parameter ranges are randomly sampled for each forward model call. Forward model calls resulting in errors or non-physical outputs are discarded (the code may still print error messages from these discarded outputs). The remainder of the script calculates 95% confidence intervals for model outputs based on the distribution of outputs, and plots these alongside proxy data and other literature estimates (see below).

%% model_functions_reformate.py
This script contains the following functions which, taken together, define and solve the forward model:

% Forward_Model - Given parameter inputs, Forward_Model calculates the initial (modern) conditions for the carbon cycle e.g. equilibrium ocean chemistry and modern fluxes. Proportionality constants for carbon cycle functions are also calculated from initial conditions. Next, the ODE solver is called, and the system of equations describing the carbon cycle are solved. The ODE solver only returns DIC and ALK as a function of time for both the ocean and the pore space. These outputs are fed back into carbon_cycle to obtain the time evolution for carbon cycle fluxes and ocean chemistry. Selected outputs are returned to Main_code_parallel.py. 

% system_of_equations - Contains the ODEs that describe the time evolution of the carbon cycle (equation S1 in manuscript). The function takes the current state of the carbon cycle (calculated using the carbon_cycle function), and returns the time derivatives of DIC and ALK in the ocean and the pore space. This function is fed into the ODE solver to compute the time evolution of DIC and ALK for the ocean and pore space.

% lf - Calculates continental land fraction as a function of time using equation S6.

% carbon_cycle - This function computes equilibrium chemistry and carbon cycle fluxes from DIC and ALK for both the ocean and the pore space. First, carbon_cycle calculates equilibrium chemistry for the ocean and the pore space using equations S12-S16 in the manuscript and heatflow and outgassing parameters from equations S8-S10. Next, it calculates surface and pore space temperatures from pCO2 using our climate model and surface-deep ocean relationship (equations S20, S4). Note that the option 'Carbon_chem' in Main_code_parallel.py determines whether this procedure is done once (fast, less accurate), or several times such that the correct temprature is converged upon for calculating carbon chemistry equilibrium constants (slow, more accuate). Finally, carbon cycle fluxes are calculated from this information. Both carbon cycle fluxes and equilibrium chemistry are returned to system_of_equations or Forward_Model.

%% thermodynamic_variables.py
Contains the function Sol_prod that calculates carbonate solubility product as a function of temperature from Pilson, M. E. (1998) "An Introduction to the Chemistry of the Sea", Prentice-Hall, Inc. This function reproduces table G.1 in Pilson for salinity=35. thermodynamic_variables.py also contains the temperature-dependent functions for carbon equilibrium constants, equil_cont.

%% clim_model.py
This contains the the parameterized climate models used in this study. clim_fun_CO2_only is the polynomial fit to the CO2-only climate model (equation S28), clim_fun_lowCH4 is the polynomial fit to the 100 ppm methane model (equation S29), and clim_fun_highCH4 is the polynomial fit to the 1% methane model (equation S30).

%% n_out_functions_no_land.py
This function is only needed for the no Archean land sensitivity tests. The nominal model rarely fails when sampling the default parameter ranges. However, for the no Archean land sensitivity test, a larger number of forward model calls result in model failure. Furthermore, model failure is slighly biased towards high Archean heatflow values. To correct for this bias, the n_out parameter is drawn from a distribution skewed towards high Archean heatflow values. This compensates for model failure and the resultant final distribution for n_out is approximately uniform, as it should be. The function no_land_no_methane should be used for the zero methane, no Archean land case (Fig. 4), whereas the function no_land_with_methane should be used for the high methane, no Archean land sensitivity test (Fig. S7).

%% Other files included in the code:
- Halevy_Bachan_high.txt and Halevy_Bachan_low.txt are ocean pH estimates from Halevy and Bachan than are plotted along our model outputs for compairson
- options_array.npy is an array created by the code for storing user choices about parallelization, methane, iteration number, and carbon chemistry precision.

END EXPLANATION OF CODE STRUCTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

-------------------------
Contact e-mail: joshkt@uw.edu

