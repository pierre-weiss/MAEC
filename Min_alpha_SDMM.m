% function [x,CF] = Min_alpha_SDMM(a_ini,beta,u1,u2,lambda,gamma,c3,c4,nit)
%
% This function finds the minimizer of:
%
%            min_a f1(a) + f2(a) + f3(a) + f4(a)
%  
% With :
%
%   - f1(a) = < exp(-A1 a) , beta > + < u1 , A1 a>
%   - f2(a) = < exp(-A2 a) , beta > + < u2 , A2 a>
%   - f3(a) = lambda * || nabla a ||_{1,eps}
%   - f4(a) = indicator of the positive orthant
%
% Given two signals u1 and u2, and an estimation of the signal u, this
% function returns an estimate of the attenuation map. The procedure
% wors with parallel proximal algorithm (SDMM).
%
%
% INPUT:
%
% - a_ini : the initial attenuation map
% - beta : original signal
% - u1, u2 : two attenuated signals
% - lambda : total variation regularisation parameter
% - gamma : proximal parameter for the SDMM
% - c3,c4 : parameters to tune the algorithm
% - nit : number of iterations
%
% OUTPUT :
%
% - x : attenuation function
% - CF : cost function wrt iterations
%
% Developer: Valentin Debarnot & Pierre Weiss, July 19 2016

function [x,CF] = Min_alpha_SDMM(a_ini,beta,u1,u2,lambda,gamma,c3,c4,nit)
x=a_ini;

y1=M1(a_ini); %M1
y2=M2(a_ini); %M2
y3_1=drond1(a_ini);
y3_2=drond2(a_ini); %TV
y4=a_ini; %positivity

z1=zeros(size(a_ini));
z2=zeros(size(a_ini));
z3_1=zeros(size(a_ini));
z3_2=zeros(size(a_ini));
z4=zeros(size(beta));

tol=1e-1;
iter_max=500;
CF=zeros(nit,1);

figure(100);colormap gray;
for k=1:nit 
    %Compute x
    tmp=M1T(y1-z1) + M2T(y2-z2) + c3*(drond1T(y3_1-z3_1)+drond2T(y3_2-z3_2)) + c4*(y4-z4);
    x=conjugate_gradient(tmp,x,c3,c4,tol,iter_max);
    
    %Compute s_i
    s1=M1(x);
    s2=M2(x);
    s3_1=c3*drond1(x);
    s3_2=c3*drond2(x);
    s4=c4*x;
        
    %Compute y_i
    y1 = prox_h1(gamma,u1,beta,s1+z1);
    y2 = prox_h1(gamma,u2,beta,s2+z2);
    
    %[y1,~]=gradient_descent_proxh1(gamma,u1,beta,s1+z1);
    %[y2,~]=gradient_descent_proxh1(gamma,u2,beta,s2+z2);
    
    y3_1=s3_1+z3_1;
    y3_2=s3_2+z3_2;
    
    ng=sqrt(y3_1.^2+y3_2.^2);
    sh=max(abs(ng)-lambda*gamma/c3,0);
    y3_1=y3_1./(max(ng,1e-16)).*sh;
    y3_2=y3_2./(max(ng,1e-16)).*sh;
    
    y4=max(s4+z4,0);
    
    %Compute z_i
    z1=z1+s1-y1;
    z2=z2+s2-y2;
    z3_1=z3_1+s3_1-y3_1;
    z3_2=z3_2+s3_2-y3_2;
    z4=z4+s4-y4;
    
    % Compute cost function
    g1=drond1(x);g2=drond2(x);
    ng=sqrt(g1.^2+g2.^2);
    
    CF(k)=sum(exp(-s1(:)).*beta(:)+u1(:).*s1(:))...
        +sum(exp(-s2(:)).*beta(:)+u2(:).*s2(:))...
        +lambda*(sum(ng(:).^2));
    
    %Plot
    if (mod(k,1000)==0)
        figure(100);colormap gray;imagesc(x);title(sprintf('ATTENUATION: Iter %i/%i -- CF:%1.2e',k,nit,CF(k)));colorbar;
        drawnow
    end
end
