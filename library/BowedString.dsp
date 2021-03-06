import("stdfaust.lib");
import("fds.lib");

//----------------------------------String Settings---------------------------//
// Generic string

//nPoints=int(Length/h);
nPoints = 100;

k = 1/ma.SR;
//Stability condition
coeff = c^2*k^2 + 4*sigma1*k;
h =sqrt((coeff + sqrt((coeff)^2 + 16*k^2*K^2))/2);

T = 150;                  // Tension [N]
radius = 3.5560e-04;      // Radius (0.016 gauge) [m]
rho = 8.05*10^3;          // Density [kg/m^3];
Area = ma.PI*radius^2;    // Area of string section
I = (ma.PI*radius^4)/ 4;  // Moment of Inertia
Emod = 174e4;             // Young modulus [Pa]
K = sqrt(Emod*I/rho/Area);// Stiffness parameter
c = sqrt(T/rho/Area);     // Wave speed
sigma1 = 0.01;            // Frequency dependent damping
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
inPoint=hslider("input point", floor(nPoints/2),0,nPoints-1,1);
outPoint=hslider("output point",floor(nPoints/2),0,nPoints-1,0.01):si.smoo;

//----------------------------------Force---------------------------------//
Vb = hslider("bow vel", 0,-10,10,0.01); //bow velocity [m/s]
Fb = 1000000; //[m/s^2]
J = Fb*k^2/den/h;
alpha=0.0001;

//----------------------------------Process---------------------------------//
//TODO: lin interp in input causes 0 output at .5 due to opposite phase
process =
    (stairsInterp1D(nPoints,inPoint):>bow(J,alpha,k,Vb)<:linInterp1D(nPoints,inPoint):
  model1D(nPoints,r,t,scheme(nPoints)))~si.bus(nPoints):linInterp1DOut(nPoints,outPoint)
    <:_,_;
