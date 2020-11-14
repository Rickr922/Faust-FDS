import("stdfaust.lib");
import("fds.lib");

//--------------------------------Model Settings-----------------------------//
k = 1/ma.SR;
K = 20;
s0 = 2;
s1 = 0.05;
c=344;

coeff = c^2*k^2+4*s1*k;
h = sqrt((coeff+sqrt(coeff*coeff+16*K*K*k*k)));

nPointsX = 30;
nPointsY = 10;

lambda = c*k/h;

//----------------------------------Equations--------------------------------//
mu=K*K*k*k/(h^4);
den = 1+s0*k;
A = 2*(1-10*mu-2*lambda*lambda-4*s1*k*k)/den;
B = (s0*k+4*k*k-1)/den;
C = (8*mu + lambda*lambda + 2*s1*k*k)/den;
D = -2*mu/den;
E = -mu/den;
F = -2*s1*k*k/den;

midCoeff = 0,0,E,0,0,
           0,D,C,D,0,
           E,C,A,C,E,
           0,D,C,D,0,
           0,0,E,0,0;

midCoeffDelay1 = 0,0,0,0,0,
                 0,0,F,0,0,
                 0,F,B,F,0,
                 0,0,F,0,0,
                 0,0,0,0,0;

r=2;
t=1;

coefficients = midCoeff,midCoeffDelay1;

scheme(pointsX,pointsY) = par (i, pointsX,
                                par(j,pointsY, coefficients));

//----------------------------------Controls---------------------------------//
inPointX=hslider("input point x", floor(nPointsX/2),0,nPointsX-1,0.01);
inPointY=hslider("input point y", floor(nPointsY/2),0,nPointsY-1,0.01);
outPointX=hslider("output point x",floor(nPointsX/2),0,nPointsX-1,0.01);
outPointY=hslider("output point y",floor(nPointsY/2),0,nPointsY-1,0.01);
hit = button("play");

forceModel = hit:ba.impulsify;

process =
    forceModel<:linInterp2D(nPointsX,nPointsY,inPointX,inPointY):
        model2D(nPointsX,nPointsY,r,t,scheme(nPointsX,nPointsY)):
            linInterp2DOut(nPointsX,nPointsY,outPointX,outPointY);
