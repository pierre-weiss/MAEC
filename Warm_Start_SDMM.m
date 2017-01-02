% function [u,a,CF] = Minimisation_u_estim_SDMM(u1,u2,dynamics,nit,alpha,gamma,ainit)
%
% Solve : 
%
%   min_a>0  sum_j <u_j, M_j a + log(exp(-M1a)+exp(-M2a))> + lambda || \nabla a||.
%
% Given two 2D images u1 and u2, this function returns an attenuation map
% a and an estimate of the non attenuated image u. The procedure works with
% a simultaneaous direction methoh of multipliers.
%
% Input : 
% - u1, u2 : two attenuated images
% - dynamics : maximum value of grey-levels
% - nit : number of iteration wanted
% - lambda : regularization parameter
% - gamma : proximal parameter
% - c3,c4 : constants
% - ainit : initialization for the attenuation map
%
% Output : 
% u : original image
% a : attenuation map reconstructed
%
% Developer: Valentin Debarnot & Pierre Weiss, December 2016

function [u,a,CF] = Warm_Start_SDMM(u1,u2,dynamics,nit,lambda,gamma,c3,c4,ainit)

if nargin<=8
    a=zeros(size(u1));
else
    a=ainit;
end
x=a;

y1=M1(a);
y2=M2(a);
y3_1=drond1(a);
y3_2=drond2(a); %TV
y4=a; %positivity

z1=zeros(size(a));
z2=zeros(size(a));
z3_1=zeros(size(a));
z3_2=zeros(size(a));
z4=zeros(size(a));


tol=1e-4*length(u1)^2;
iter_max=500;

CF=zeros(nit,1);

for k=1:nit
    %Compute x
    x=conjugate_gradient_initialization(c3*(drond1T(y3_1-z3_1)+drond2T(y3_2-z3_2))+c4*(y4-z4)+M1T(y1-z1)+M2T(y2-z2),x,c3,c4,tol,iter_max);
    
    %Compute s_i
    s1=M1(x);
    s2=M2(x);
    s3_1=c3*drond1(x);
    s3_2=c3*drond2(x);
    s4=c4*x;
    
    %CF
    e1=exp(-s1);e2=exp(-s2);
    logz1z2=log(e1(:)+e2(:));
    g1=drond1(x);g2=drond2(x);
    ng=sqrt(g1.^2+g2.^2);
    CF(k)=sum(u1(:).*(s1(:)+logz1z2)) + sum(u2(:).*(s2(:)+logz1z2)) + lambda*sum(ng(:));
  
    %Compute y_i
    [y1,y2] = prox_g1(s1+z1,s2+z2,u1,u2,gamma);    
    y3_1=s3_1+z3_1;
    y3_2=s3_2+z3_2;
    
    ng=sqrt(y3_1.^2+y3_2.^2);
    sh=max(abs(ng)-lambda*gamma/c3,0);
    y3_1=y3_1./(max(ng,1e-16)).*sh;
    y3_2=y3_2./(max(ng,1e-16)).*sh;
    
    y4=min(max(s4+z4,0),dynamics);
    
    %Compute z_i
    z1=z1+s1-y1;
    z2=z2+s2-y2;
    z3_1=z3_1+s3_1-y3_1;
    z3_2=z3_2+s3_2-y3_2;
    z4=z4+s4-y4;
    
    %Plot
    if (mod(k,100)==0)
        figure(100);colormap gray;imagesc(x);title(sprintf('ATTENUATION: Iter %i/%i -- CF:%1.2e',k,nit,CF(k)));colorbar;
        u=(u1+u2)./(e1+e2);
        figure(101);colormap gray;imagesc(u);title(sprintf('IMAGE: Iter %i/%i -- CF:%1.2e',k,nit,CF(k)));colorbar;
        drawnow
    end
end
u=(u1+u2)./(e1+e2);
a=x;