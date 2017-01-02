function [x1,x2] = prox_g1(z1,z2,u1,u2,gamma)

a=gamma*(u1+u2);

y1=gamma*u1-z1;
y2=gamma*u2-z2;
[x1,x2]=prox_lse_2D(a,y1,y2);
x1=-x1;
x2=-x2;

