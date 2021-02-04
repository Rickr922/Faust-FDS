import("stdfaust.lib");

/*
u_l^n+1 = 2(1-lambda^2)u_l^n - u_l^n-1 + lambda^2(u_l+1^n + u_l-1^n)
*/

k=1/ma.SR;
c=344;
h=c*k;
lambda=c*k/h;

A = 2*(1-lambda^2);
B = lambda^2;
C = -1;
 
r=1; t=1;

midCoeff = B,A,B;


midCoeffDel = 0,C,0;

/*
u_l^n+1 = 2(1-lambda^2)u_l^n - u_l^n-1 + 2*lambda^2(u_l+1^n)*/

D = 2*lambda^2;
leftCoeff = 0,A,D;

scheme(points) = leftCoeff,midCoeffDel,par(i,points-1,midCoeff,midCoeffDel);


play = button("Play") : ba.impulsify;
inPoint = hslider("Input Point", floor(nPoints/2),0,nPoints-1,0.01);
outPoint = hslider("Output Point", floor(nPoints/2),0,nPoints-1,0.01);

process = play<:fd.linInterp1D(nPoints,inPoint):fd.model1D(nPoints,r,t,scheme(nPoints)):fd.linInterp1DOut(nPoints,outPoint);

nPoints = 100;
