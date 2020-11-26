import("stdfaust.lib");
import("fds.lib");

nPoints = 90;
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
play = button("hit"):ba.impulsify;
inPoint=hslider("input point", floor(nPoints/2),0,nPoints-1,0.01);
outPoint=hslider("output point",floor(nPoints/2),0,nPoints-1,0.01):si.smoo;

//----------------------------------Force---------------------------------//
JCoeff = (k^2/den/rho/Area);

KH = 300;
mH = 1;
omega0SqrH = KH/mH;
sigma0H = 14;
alpha = 2.5;
offset=1;
forceScaling = 0.009;

gain = 500;
process = (linInterp1D(nPoints,inPoint):>
    nlHammer(omega0SqrH,sigma0H,10000,alpha,k,offset,play*forceScaling)*(JCoeff)<:
        linInterp1D(nPoints,inPoint):
            model1D(nPoints,r,t,scheme(nPoints)))~si.bus(nPoints):
                linInterp1DOut(nPoints,outPoint)*gain<:_,_;
