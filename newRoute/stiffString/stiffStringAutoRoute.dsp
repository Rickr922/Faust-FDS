import("stdfaust.lib");

//--------------------------------Model Settings-----------------------------//
nPoints = 100;
L = 0.1;                  // String length [m]
//nPoints=int(L/h);

k = 1/48000;//1/ma.SR;
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
nInputs = inputs(schemeMidPoint);

den = 1+sigma0*k;
A = (2*h^4-2*c^2*k^2*h^2-4*sigma1*k*h^2+6*K^2*k^2)/den/h^4;
B = (sigma0*k*h^2-h^2+4*sigma1*k)/den/h^2;
C = (c^2*k^2*h^2+2*sigma1*k*h^2-4*K^2*k^2)/den/h^4;
D = -2*sigma1*k/den/h^2;
E = K^2*k^2/den/h^4;

//west=left, east=right
schemeMidPoint(u_ww,u_w,u,u_e,u_ee,fIn) =
        A*u + B*u' + C*(u_e+u_w) + D*(u_e'+u_w') + E*(u_ee+u_ww) + fIn;

schemeFixedPoint(u_ww,u_w,u,u_e,u_ee,fIn) = 0;

//----------------------------------Controls---------------------------------//
hit = button("hit"):ba.impulsify;
stop = button("Stop");
inPoint=hslider("input point", floor(nPoints/2),0,nPoints-1,0.01);
outPoint=hslider("output point",floor(nPoints/2),0,nPoints-1,0.01):si.smoo;

//----------------------------------Force---------------------------------//
forceModel = hit;
linInterp1DForce(i,inPoint,force) = force*select2(
    i==int(inPoint), select2(i==int(inPoint+1),0,fraction),(1-fraction))
with
{
    fraction = ma.frac(inPoint);
};

linInterp1D(selectedPoint,nPoints) = par(i,nPoints,_*select2(
        i==int(selectedPoint), select2(i==int(selectedPoint+1),0,fraction),(1-fraction)))
with
{
    fraction = ma.frac(outPoint);
};
//----------------------------------Output-------------------------------//
linInterp1DOut(outPoint,nPoints) = linInterp1D(outPoint,nPoints):>_;

//----------------------------------Build Model-------------------------------//
buildScheme(nPoints) = par (i, nPoints, schemeMidPoint);

  /*schemeFixedPoint,
    par (i, nPoints-2, schemeMidPoint),
    schemeFreePointEast;*/

routing(nPoints,nInputs) =   route(nPoints+nPoints, nPoints*nInputs,
    par(x, nPoints, connections(x)))
with
{
    connections(x) = par(i, nInputs-1, U(x),I(x-spaceDepX+i,nInputs-2-i)),
                                       U(x)+nPoints, I(x,nInputs-1);
    U(x) = x+1;
    I(x,j)  = (1 + j + (x*nInputs)) * (x>=0) * (x<nPoints);
    spaceDepX = (nInputs - 2)/2;
};

model(nPoints) = (routing(nPoints,nInputs) : buildScheme(nPoints)) ~
    (par(i, nPoints, _*(stop==0)));

process = forceModel<:linInterp1D(inPoint,nPoints): model(nPoints):linInterp1DOut(outPoint,nPoints)<:_,_;
