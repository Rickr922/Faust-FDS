import("stdfaust.lib");

//--------------------------------Model Settings-----------------------------//
nPoints = 6;

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

//----------------------------------Library-------------------------------//
schemePoint1D(R,T,fIn) = coeffs,neighbors<:
    sum(t,T+1,
        sum(i,nNeighbors,
            ba.selector(int(i+t*nNeighbors),nNeighbors*(T+1),coeffs)*
            ba.selector(i,nNeighbors,neighbors)@(t)))
                + fIn
with
{
    nNeighbors = (2*R+1);
    neighbors = si.bus(nNeighbors);
    coeffs=si.bus(nNeighbors*(T+1));
};

buildScheme1D(points,R,T) =
    par (x, points,schemePoint1D(R,T));

//----------------------------------Interpolation---------------------------------//
linInterp1D(points,point) = par(i,points,_*select2(
        i==int(point), select2(i==int(point+1),0,fraction),(1-fraction)))
with
{
    fraction = ma.frac(point);
};

//----------------------------------Force---------------------------------//
forceModel = play:ba.impulsify;

stairsForce(points,point) = ba.selectoutn(points,point);

linInterp1DOut(points,point) = linInterp1D(points,point):>_;

route1D(points, R, T) = route(points*2+points*nCoeffs, points*nInputs,
                                par(x, nPoints, connections(x)))
with
{
    connections(x) =  par(k,nCoeffs,x*nCoeffs+k+1,C(x,k+1)),
                      P(x) + points, C(x,0),
                      par(i, nNeighbors, P(x),C(x-R+i,nInputs-1-i));

    P(x) = x+1 + nCoeffs*points;
    C(x,count) = (1 + count + (x*nInputs)) * (x>=0) * (x<points);

    nNeighbors = 2*R+1;
    nCoeffs = nNeighbors*(T+1);
    nInputs = nNeighbors+1+nCoeffs;
};

model1D(points,R,T,scheme) =
    (route1D(points,R,T,scheme) : buildScheme1D(points,R,T)) ~ si.bus(points);

process = forceModel<:linInterp1D(nPoints,inPoint):
  model1D(nPoints,r,t,scheme(nPoints)):
  linInterp1DOut(nPoints,outPoint)<:_,_;
