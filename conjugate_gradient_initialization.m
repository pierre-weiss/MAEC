% function x = conjugated_gradient_initialization(b,x_ini,c3,c4,tol)
%
% 1D solver
%
% Compute the solution of linear system with matrix:
%        M1TM1+M2TM2 + nablaT nabla + I,
%
% INPUT:
%
% - b : second member
% - x_ini : first guess for x
% - tol : tolerance from error
%
% OUTPUT :
%
% - x approximation of the linear system
% 
% Developer: Valentin Debarnot, August 1 2016

function x = conjugate_gradient_initialization(b,x_ini,c3,c4,tol,iter_max)
n1=size(x_ini,1);
n2=size(x_ini,2);
x=x_ini(:);
Ax=(M1T(M1(x_ini))+M2T(M2(x_ini))+c3^2*(drond1T(drond1(x_ini))+drond2T(drond2(x_ini)))+c4^2*x_ini);
r=b(:)-Ax(:);
p=r;

i=1;
while norm(r)>tol && i<iter_max
    p=reshape(p,n1,n2);
    Apk=(M1T(M1(p))+M2T(M2(p))+c3^2*(drond1T(drond1(p))+drond2T(drond2(p)))+c4^2.*p);
    p=p(:);
    s1 = (r'*r) ./ (p'*Apk(:));
    x=x+s1*p;
    r_=r;
    r=r-s1*Apk(:);
    if norm(r)>tol
        s2=(r'*r)./(r_'*r_);
        p=r+s2*p;
    end
%     nr(i)=norm(r);
    i=i+1;
end
% if i>1
%     plot(nr);
% end
x=reshape(x,n1,n2);
% sprintf('Convergence in %i iterations with residual norm equal to : %d',i,norm(r))


