# TuTraSt

Matlab code for TuTraST analysis as presented in [1]. 

TuTraSt is a novel algorithm to predict self-diffusion of a mobile guest particle in a crystalline material. It detects the energies at which diffusion paths are formed, allowing for easy identification of diffusive systems, and furthermore partitions the potential energy field into energy basins and transitions states. This TUnnel and TRAnsition STate search algorithm permits a transition state theory based analysis for fast prediction of the diffusion coefficients with an automated multiscale modeling approach.

TuTraSt_main.m executes the analysis reading grid.cube and input.param as input.

grid.cube is a the potential energy grid file in cube-format. Please check that the formatting of your file matches the example file given. 

input.param contains the user defined parameters for the TuTraSt analysis and kMC simualtions. The ordering in this file must be as defined. 

For questions on technical issues, usage, maintenance and development please contact Amber Mace (amber.mace@epfl.ch). We are also happy discuss any ideas new applications or potential areas of development, so free to contact us for any related matters. 

Please cite!
[1] Amber Mace, Senja Barthel and Bernd Smit, J. Chem. Theory Comput. 2019,  15, 4, 2127-2141
