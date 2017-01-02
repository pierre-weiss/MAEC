% function [x,CF] = prox_f1(y,u,u1,gamma,nit)
%
% Find the solution of
%
%           arg min_x  < exp(-x) , u >
%               + < u1 , x > +1/(2gamma)*||x-y||_2^2
%
% INPUT:
%
% - y : variable of the proximal operator of the objective function
% - u : the original image
% - u2 : an attenuated images
% - gamma : proximal parmeter
% - nit : number of iteration
%
% OUTPUT :
%x
% - x : the value wich reach the minimum.
% - CF : cost function wrt iterations
%
% Developer: Valentin Debarnot, June 27 2016

function x = prox_f1(y,u,uj,gamma,table_prox,v_discr,c_discr)

v=gamma*u;
c= - (gamma*uj-y);
x_approx = find_prox_f(table_prox,v,c,v_discr,c_discr);

%% Dichotomie classique
%     %Assume x lives in [a,b]
% a=x_approx-0.4;
% b=x_approx+0.4;
% %     %/!\ v=gamma*u
% %     %    c=gamma*u1-y
% for i=1:45
%     m=(a+b)/2;
%     f1pa=-exp(-a).*(gamma*u)+a-c;
%     f1pm=-exp(-m).*(gamma*u)+m-c;
%     I1=f1pa.*f1pm<0;
%     I2=f1pa.*f1pm>=0;
%     a(I2)=m(I2);
%     b(I1)=m(I1);
% end
% x=m;

%% Jonas' method
% x=x_approx;
% fgrad=-gamma*exp(-x).*u+uj*gamma+(x-y);
% norm_grad=norm(fgrad);
% nit=1;
% while norm_grad > 1e-10 && nit<10
%     x=gamma*exp(-x).*u+c;
%     fgrad=-gamma*exp(-x).*u+uj*gamma+(x-y);
%     norm_grad=norm(fgrad);
%     nit=nit+1;
% end

%% Newton method
iter=1;
nit=10;
x=x_approx;
while (iter<nit)
    f1p=-exp(-x).*v+x-c;
    f1pp=exp(-x).*v+1;
    %CF(iter)=sum(sum(exp(-x).*v+gamma*uj.*x))+sum(sum((x-y).^2));
    %Iterations
    x=x-f1p./f1pp;
    iter=iter+1;
end