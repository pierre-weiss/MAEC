% function [beta,CF] = Min_beta_SDMM(alpha,beta_ini,u1,u2,lambda,nit,c1,c2)
%
% This function finds the minimizer of:
%
% min_u < exp(-A1 alpha) , beta > - < u1 , log(beta) > + < exp(-A2 alpha ) , beta  > - < u2 ,
% log(beta ) >  + lambda * || nabla beta ||_{1,eps}
%
% Given two 2D images u1 and u2, and an attenuation map alpha, this function,
% returns an estimate of the non attenuated image beta . The procedure works
% with the SDMM algorithm. 
%
% INPUT:
%
% - alpha : attenuation map.
% - beta : starting value of the signals u.
% - u1, u2 : two attenuated images
% - lambda : total variation regularisation parameter
% - nit : number of iterations
% - c1,c2 : parameter to tune the algorithm
%
% OUTPUT : 
%
% - beta : denoised image
% - CF : cost function
%
% Developer: Pierre Weiss, December 28 2016

function [beta,CF] = Min_beta_SDMM(alpha,beta_ini,u1,u2,lambda,nit,c1,c2)

%Cost function
CF=zeros(nit,1);
gamma=1;

% Definition of linear operators
der1=zeros(size(u1));
der1(1,1)=1;der1(end,1)=-1;
der2=zeros(size(u1));
der2(1,1)=1;der2(1,end)=-1;
Der1=fft2(der1);Der2=fft2(der2);
d1=@(x) ifft2(Der1.*fft2(x));
d2=@(x) ifft2(Der2.*fft2(x));
d1T=@(x) ifft2(conj(Der1).*fft2(x));
d2T=@(x) ifft2(conj(Der2).*fft2(x));
denom=c1^2+c2^2*(abs(Der1).^2+abs(Der2).^2);

%Initialisation (should be modified to take beta_ini in account)
y1=zeros(size(u2));
y21=zeros(size(u2));
y22=zeros(size(u2));

z1=zeros(size(u2));
z21=zeros(size(u2));
z22=zeros(size(u2));

%Useful notations
a=exp(-M1(alpha))+exp(-M2(alpha));
u=u1+u2;

%% SDMM algorithm
figure(101);
for i=1:nit
    %% Step 1
    tmp=c1*(y1-z1)+c2*(d1T(y21-z21)+d2T(y22-z22));
    x=ifft2(fft2(tmp)./denom);
        
    %% Step 2
    s1=c1*x;
    tmp=s1+z1;
    y1=(-(gamma*a/c1-tmp)+sqrt((gamma*a/c1-tmp).^2+4*gamma*u))/2;
    z1=z1+s1-y1;
    
    s21=c2*d1(x);s22=c2*d2(x);
    tmp1=s21+z21;tmp2=s22+z22;
    ng=sqrt(tmp1.^2+tmp2.^2);
    ng(ng==0)=1;
    y21=tmp1./ng.*max(ng-gamma*lambda/c2,0);
    y22=tmp2./ng.*max(ng-gamma*lambda/c2,0);
    z21=z21+s21-y21;
    z22=z22+s22-y22;
    
    %% Step 3 Cost function
    x(x<1e-16)=1e-16;
    CF(i)=sum(a(:).*x(:)-u(:).*log(x(:)))+lambda*sum(ng(:));
    
    %% Display
    if (mod(i,10000)==0)
        figure(101);colormap gray;imagesc(x);title(sprintf('IMAGE: Iter %i/%i -- CF:%1.2e',i,nit,CF(i)));colorbar;
        drawnow
        %pause;
    end
end
beta=x;
close(101); 
