%% This file compiles all c-mex functions. You will need a proper installation of a compiler.
% If you have troubles compiling, please install a correct compiler by
% following Matlab's instructions. If you struggle with openmp flags, you
% can remove CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp",
% at the cost of not exploiting all the cores of your machine.

disp('Compiling prox_lse_2D_c.cpp')
mex prox_lse_2D.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"

disp('Compiling lambertw2_c.cpp')
mex lambertw2_c.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"

disp('Compiling lambertw_c.cpp')
mex lambertw_c.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"

disp('Success!')