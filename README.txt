%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This toolbox was developed to support the companion paper:
% "Multiview Attenuation Estimation and Correction" 
% by V. Debarnot, J. Kahn and P. Weiss, January 2017.
% Please cite this paper in case you found this toolbox useful for your own research.
%
% The toolbox provides various tools to:
% i) recover attenuation from two attenuated 1D or 2D signals
% ii) correct attenuation given two attenuated images and an attenuation map
% iii) compute the prox of logsumexp (mex function in c++ with openmp support)
% iv) compute the LambertW function efficiently (mex function in c++ with openmp support)
% v) various examples of how to use the toolbox
%
% Developers: Valentin Debarnot and Pierre Weiss, 2016.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

To start, open matlab and type
>> compile_all
this will compile all mex files on your system. By default, only the linux64 architecture should be working.

Then, various demos are available
>> demo_first
>> demo_warmstart
>> demo_correct_attenuation
>> demo_alternate_minimization

Main functions: 

- Warm_Start_SDMM : computes the warm start solution
- Min_alpha_SDMM : evaluates attenuation only knowing density
- Min_beta_SDMM : evaluates density only knowing attenuation
- lambertw2_c.cpp : solves x exp(x)=exp(alpha) (openmp support)
- prox_lse_2D.cpp : finds the proximal operator of the logsumexp function (openmp support)