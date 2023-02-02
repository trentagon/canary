
    _____          _   _          _______     __
   / ____|   /\   | \ | |   /\   |  __ \ \   / /   \\
  | |       /  \  |  \| |  /  \  | |__) \ \_/ /    (o>
  | |      / /\ \ |   ` | / /\ \ |  _  / \   /  \\_//)
  | |____ / ____ \| |\  |/ ____ \| | \ \  | |    \_/_)
   \_____/_/    \_\_| \_/_/    \_\_|  \_\ |_|     _|_

A Mars Atmospheric Evolution Model that tracks [CA]rbon, [N]itrogen, and [AR]gon, [Y]ay!                                       
Developed by Trent Thomas and Renyu Hu 

This set of Matlab codes runs our model for the evolution of carbon, nitrogen, and argon in the Martian atmosphere over the last 3.8 billion years - A.K.A. “CANARY”. A full model description is found in:

 “Constraints on the Size and Composition of the Ancient Martian Atmosphere from Coupled CO2-N2-Ar Isotopic Evolution Models”, by Trent B. Thomas, Renyu Hu, and Daniel Y. Lo; 2023, Planetary Science Journal. 

As a matter of courtesy, we request that people using this code please cite Thomas et al. (2023). In the interest of an "open source" approach, we also request that authors who use and modify the code, please send a copy of papers and modified code to the authors (renyu.hu@jpl.nasa.gov and tbthomas@uw.edu)

REQUIREMENTS: Matlab (this code was written in version 2018a)

HOW TO RUN CODE:
(1) Download canary and make sure everything is in the same directory.
(2) Open “canary_fly.m”
	—> Set the path to the directory containing CANARY.
	—> Choose the input file (the default file recreates the case in figure 5).
	—> Modify parameters of your choice.
(3) Run “canary_fly.m”. Code will output plots of the evolution and print results in the Matlab console. 


%% canary_fly.m
This is the top-level script to run the model. Set the path, choose an input file, run the script.

%% inputFile_fig5.m
This is the main input file that sets relevant parameters and values used in the model. It recreates figure 5 in the main text. Use this file to explore how results change with parameters.

%% Everything else
Model functions are contained in the “functions” subdirectory. External data used to run the model is contained in the subdirectory “input_data”.