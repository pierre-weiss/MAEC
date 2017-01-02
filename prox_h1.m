% This function returns the solution of 
% min_z gamma* ( uz + exp(-z)beta ) + 1/2 (z-z0)^2
function z = prox_h1(gamma,u,beta,z0)

a=z0-gamma.*u;
%z=a+lambertW(gamma.*beta.*exp(-(a)));
%z=a+lambertw_c(gamma.*beta.*exp(-(a)));

z=a+lambertw2_c(log(gamma.*beta)-a);


