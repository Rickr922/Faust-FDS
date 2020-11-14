import("stdfaust.lib");
import("fds.lib");
//--------------------------------Model Settings-----------------------------//
nPoints = 90;

k = 1/ma.SR;
c = 344;
h = c*k;
s0 = 500;
lambda = c*k/h;

//----------------------------------Equations--------------------------------//
C1 = 1+(2*s0*k*k/h);
C2 = (2*s0*k*k/h)-1;
A = 2*(1-lambda*lambda)/C1;
B = lambda*lambda/C1;
C = C2/C1;

midCoeff = B,A,B;
midCoeffDel = 0,C,0;

r=1;
t=1;

scheme(points) = par(i,points,midCoeff,midCoeffDel);

/*
//west=left, east=right
schemeMidPoint(fIn,u_w,u,u_e) =
        u*2*(1-lambda*lambda)/C1+u'*C2/C1+ lambda*lambda*(u_w+u_e)/C1 + fIn;

schemeFreePointWest(fIn,u_w,u,u_e) =
        u*2*(1-lambda*lambda)/C1+u'*C2/C1+ lambda*lambda*2*u_e/C1 + fIn;

schemeFreePointEast(fIn,u_w,u,u_e) =
        u*2*(1-lambda*lambda)/C1+u'*C2/C1+ lambda*lambda*2*u_w/C1 + fIn;*/

//----------------------------------Controls---------------------------------//
play = button("hit");
stop = button("Stop");
inPoint=hslider("input point", floor(nPoints/2),0,nPoints-1,0.01);
outPoint=hslider("output point",floor(nPoints/2),0,nPoints-1,0.01):si.smoo;

forceModel = play:ba.impulsify;

process = forceModel<:linInterp1D(nPoints,inPoint):
  model1D(nPoints,r,t,scheme(nPoints)):
  linInterp1DOut(nPoints,outPoint)<:_,_;
