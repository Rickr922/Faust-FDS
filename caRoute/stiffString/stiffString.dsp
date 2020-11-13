import("stdfaust.lib");

nPoints = 10;
L = 0.1;                  // String length [m]
//nPoints=int(L/h);

k = 1/ma.SR;
//Stability condition
coeff = c^2*k^2 + 4*sigma1*k;
h =sqrt((coeff + sqrt((coeff)^2 + 16*k^2*K^2))/2);

T = 150;                 // Tension [N]
//T = hslider("Tension",150,10,1000,0.1);
radius = 3.5560e-04;    // Radius (0.016 gauge) [m]
rho = 8.05*10^3;        // Density [kg/m^3];
Area = ma.PI*radius^2;        // Area of string section
I = (ma.PI*radius^4)/ 4;   // Moment of Inertia
Emod = 174e4;              // Young modulus [Pa]
K = sqrt(Emod*I/rho/Area);    // Stiffness parameter
c = sqrt(T/rho/Area);      // Wave speed
sigma1 = 0.01;         // Frequency dependent damping
sigma0 = 0.0005;

//----------------------------------Equations--------------------------------//
den = 1+sigma0*k;
A = (2*h^4-2*c^2*k^2*h^2-4*sigma1*k*h^2+6*K^2*k^2)/den/h^4;
B = (sigma0*k*h^2-h^2+4*sigma1*k)/den/h^2;
C = (c^2*k^2*h^2+2*sigma1*k*h^2-4*K^2*k^2)/den/h^4;
D = -2*sigma1*k/den/h^2;
E = K^2*k^2/den/h^4;

midCoeff = E,C,A,C,E;
midCoeffDel = 0,D,B,D,0;

r=2;
t=1;

scheme(points) = par(i,points,midCoeff,midCoeffDel);


//----------------------------------Controls---------------------------------//
play = button("hit");
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
