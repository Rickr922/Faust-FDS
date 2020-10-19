import("stdfaust.lib");

//--------------------------------Model Settings-----------------------------//
nPoints = 30;

//modelType: 1->free-fixed; 2->fixed-free; else->fixed-fixed
modelType = 0;

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
L = 1;                  // String length [m]
Emod = 174e4;              // Young modulus [Pa]
K = sqrt(Emod*I/rho/Area);    // Stiffness parameter
c = sqrt(T/rho/Area);      // Wave speed
sigma1 = 0.01;         // Frequency dependent damping
sigma0 = 0.0005;

//----------------------------------Equations--------------------------------//
nPointInputs = inputs(midPointEq);
spaceDep = (nPointInputs-1)/2;    //n° of spatial side points needed by the update eq

den = 1+sigma0*k;
A = (2*h^4-2*c^2*k^2*h^2-4*sigma1*k*h^2+6*K^2*k^2)/den/h^4;
B = (sigma0*k*h^2-h^2+4*sigma1*k)/den/h^2;
C = (c^2*k^2*h^2+2*sigma1*k*h^2-4*K^2*k^2)/den/h^4;
D = -2*sigma1*k/den/h^2;
E = K^2*k^2/den/h^4;

midPointEq(i, uSide_ll, uSide_l, uSide_r, uSide_rr) = equation
with
{
    fIn = linInterpolForce(i,inPoint,forceModel);
    equation = u
    letrec //u_(l-1)=uSide_l,u_(l+1)=uSide_r
    {
        'u = A*u + B*u' + C*(uSide_r+uSide_l) + D*(uSide_r'+uSide_l') + E*(uSide_rr+uSide_ll) + fIn;
    };
};


//NOT WORKING
freePointEqL(fIn, uSide_ll, uSide_l, uSide_r, uSide_rr) = u
  letrec
  {
    'u = A*u + B*u' + C*(uSide_r+uSide_l) + D*(uSide_r'+uSide_l') + E*(uSide_rr+uSide_ll) + fIn;
  };

freePointEqR(fIn, uSide_ll, uSide_l, uSide_r, uSide_rr) = u
  letrec
  {
    'u = A*u + B*u' + C*(uSide_r+uSide_l) + D*(uSide_r'+uSide_l') + E*(uSide_rr+uSide_ll) + fIn;
  };

fixedPoint(i, uSide_ll, uSide_l, uSide_r, uSide_rr) = 0;

buildStringDamp(nPoints,modelType) = par(i, nPoints, midPointEq(i));

//----------------------------------Controls---------------------------------//
hit = button("hit"):ba.impulsify;
stop = button("Stop");
inPoint=hslider("input point", floor(nPoints/2),0,nPoints-1,0.01);
outPoint=hslider("output point",floor(nPoints/2),0,nPoints-1,0.01):si.smoo;

//----------------------------------Force---------------------------------//
forceModel = hit;
linInterpolForce(i,inPoint,force) = force*select2(
    i==int(inPoint), select2(i==int(inPoint+1),0,fraction),(1-fraction))
with
{
    fraction = ma.frac(inPoint);
};

//NOT WORKING
raisedCosine(0) = 0;
raisedCosine(n) = 0.5*(1+cos(ma.PI*(n-inPoint)/0.1));

//----------------------------------Out Points-------------------------------//
linInterpolOut(outPoint) = par(i,nPoints,_*select2(
        i==int(outPoint), select2(i==int(outPoint+1),0,fraction),(1-fraction))):>_
with
{
    fraction = ma.frac(outPoint);
};

//----------------------------------Build Model-------------------------------//
model = (route(nPoints,nPoints*spaceDep*2,
            par(i,nPoints,
                par(j, spaceDep,
                    (i+1, 2*spaceDep*i + (2*spaceDep-1)*j - 2*(spaceDep^2-spaceDep)),
                    (i+1, 2*spaceDep*i + (2*spaceDep-1)*j - 2*(spaceDep^2-spaceDep) +
                    2*spaceDep^2+spaceDep)))):
                        buildStringDamp(nPoints,modelType))~
                            (si.bus(nPoints)):
                                si.bus(nPoints):linInterpolOut(outPoint);

process = model<:_,_;
//process = buildStringDamp(nPoints,modelType);
