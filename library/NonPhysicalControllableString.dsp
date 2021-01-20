import("stdfaust.lib");

/*DISCLAMER:
    This is a creative "non physical" finite difference scheme physical model
    of a string, intended to show how changing the different physical 
    parameters has an impact on the sounding characteristics of the string. 
    
    My scientific consciousness forces me to say that to make things 
    physically correct, the number of string points should change according 
    to the variations of each parameter. However, this cannot be done at run 
    time, so physically correct models can become a bit boring.

    Have fun!
*/

//----------------------------------String Settings---------------------------//
//nPoints=int(Length/h);
nPoints = 200;

k = 1/ma.SR;
//Stability condition
coeff = c^2*k^2 + 4*sigma1*k;
h = sqrt((coeff + sqrt((coeff)^2 + 16*k^2*K^2))/2);

T = hslider("String Tension (N)", 150,20,1000,0.1);                     // Tension [N]
radius = hslider("String Radius (m)", 3.6e-04,2e-5,1e-3,0.00001);       // Radius (0.016 gauge) [m]
rho = hslider("String Material Density (kg/m^3)", 8.05*10^3,1e1,1e6,1); // Density [kg/m^3];
Emod = hslider("String Young Modulus (Pa)",174e4,1e2,1e8,1);            // Young modulus [Pa]
Area = ma.PI*radius^2;                                                  // Area of string section
I = (ma.PI*radius^4)/ 4;                                                // Moment of Inertia
K = sqrt(Emod*I/rho/Area);                                              // Stiffness parameter
c = sqrt(T/rho/Area);                                                   // Wave speed
sigma1 = 0.01;                                                          // Frequency dependent damping
sigma0 = 0.0005;                                                        // Frequency independent damping

//----------------------------------Equations--------------------------------//
den = 1+sigma0*k;
A = (2*h^4-2*c^2*k^2*h^2-4*sigma1*k*h^2+6*K^2*k^2)/den/h^4;
B = (sigma0*k*h^2-h^2+4*sigma1*k)/den/h^2;
C = (c^2*k^2*h^2+2*sigma1*k*h^2-4*K^2*k^2)/den/h^4;
D = -2*sigma1*k/den/h^2;
E = K^2*k^2/den/h^4;

midCoeff = E,C,A,C,E;
midCoeffDel = 0,D,B,D,0;

r = 2;
t = 1;

scheme(points) = par(i,points,midCoeff,midCoeffDel);

//----------------------------------Controls---------------------------------//
play = button("Play");
inPoint = hslider("Input Point",floor(nPoints/2),0,nPoints-1,0.01);
outPoint = hslider("Output Point",floor(nPoints/2),0,nPoints-1,0.01):si.smoo;

//----------------------------------Force---------------------------------//
forceModel = play:ba.impulsify;

//----------------------------------Process---------------------------------//
process = forceModel<:fd.linInterp1D(nPoints,inPoint):
  fd.model1D(nPoints,r,t,scheme(nPoints)):
  fd.linInterp1DOut(nPoints,outPoint)<:_,_;
