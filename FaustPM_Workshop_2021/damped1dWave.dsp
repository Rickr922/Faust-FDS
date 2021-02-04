/*
NOTE: There is an error both in the slides and in the code used in the video: the coefficients C1 and C2 are wrong, here you can find the correct ones.
There is not much perceptual difference, but it's error in the discretization of the PDE.
*/

import("stdfaust.lib");

k=1/ma.SR;
c=344;
h=c*k;
lambda=c*k/h;

play = button("Play") : ba.impulsify;
inPoint = hslider("Input Point", floor(nPoints/2),0,nPoints-1,0.01);
outPoint = hslider("Output Point", floor(nPoints/2),0,nPoints-1,0.01);

process = play<:fd.linInterp1D(nPoints,inPoint):fd.model1D(nPoints,r,t,scheme(nPoints)):fd.linInterp1DOut(nPoints,outPoint)<:_,_;

/*
u_l^n+1 = 2(1-lambda^2)/C1 u_l^n +C2/C1 u_l^n-1 + lambda^2/C1 (u_l+1^n + u_l-1^n)
C1 = sigma_0 k + 1;
C2 = sigma_0 k - 1;*/

sigma0 = 5;

r=1; t=1;

A = 2*(1-lambda^2)/C1;
B = lambda^2/C1;
C = C2/C1;

C1 = (sigma0*k) + 1;
C2 = (sigma0*k) - 1;

midCoeff = B,A,B;
midCoeffDel = 0,C,0;

scheme(points) = par(i,points, midCoeff, midCoeffDel);

nPoints=100;