import("stdfaust.lib");

//--------------------------------Model Settings-----------------------------//
nPoints = 80;

//modelType: 1->free-fixed; 2->fixed-free; else->fixed-fixed
modelType = 0;

k = 1/ma.SR;
//Stability condition
coeff = c^2*k^2 + 4*sigma1*k;
h =sqrt((coeff + sqrt((coeff)^2 + 16*k^2*K^2))/2);

//T = 150;                 // Tension [N]
T = hslider("Tension",150,10,1000,0.1);
radius = 3.5560e-04;    // Radius (0.016 gauge) [m]
rho = 8.05*10^3;        // Density [kg/m^3];
Area = ma.PI*radius^2;        // Area of string section
I = (ma.PI*radius^4)/ 4;   // Moment of Inertia
L = 1;                  // String length [m]
Emod = 174e6;              // Young modulus [Pa]
K = sqrt(Emod*I/rho/Area);    // Stiffness parameter
c = sqrt(T/rho/Area);      // Wave speed
sigma1 = 0.01;         // Frequency dependent damping
sigma0 = 0.0005;

//----------------------------------Equations--------------------------------//
nPointInputs = inputs(midPointEq);
nSidePoints = nPointInputs-1; //twice the number of side points (left + right)
spaceDep = nSidePoints/2;    //n° of spatial side points needed by the update eq

den = 1+sigma0*k;
A = (2*h^4-2*c^2*k^2*h^2-4*sigma1*k*h^2+6*K^2*k^2)/den/h^4;
B = (sigma0*k*h^2-h^2+4*sigma1*k)/den/h^2;
C = (c^2*k^2*h^2+2*sigma1*k*h^2-4*K^2*k^2)/den/h^4;
D = -2*sigma1*k/den/h^2;
E = K^2*k^2/den/h^4;

midPointEq(fIn, uSide_ll, uSide_l, uSide_r, uSide_rr) = u //u_(l-1)=uSide_l,u_(l+1)=uSide_r
  letrec
  {
    'u = A*u + B*u' + C*(uSide_r+uSide_l) + D*(uSide_r'+uSide_l') + E*(uSide_rr+uSide_ll) + fIn;
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

fixedPoint(fIn, uSide_ll, uSide_l, uSide_r, uSide_rr) = 0;

buildStringDamp(nPoints,modelType) = par(i, nPoints, midPointEq);

/*
  //modelType: 1->free-fixed; 2->fixed-free; else->fixed-fixed
  par(i,nSidePoints,
        (si.bus(nPointInputs)<:(midPointEq, freePointEqL:select2(modelType==1)))),
  par(i, nPoints-nSidePoints, midPointEq),
  par(i,nSidePoints,
        (si.bus(nPointInputs)<:(midPointEq, freePointEqR:select2(modelType==2))));
          //par(i, nPoints, midPointEq);*/
//----------------------------------Controls---------------------------------//
hit = button("hit"):ba.impulsify;
stop = button("Stop");
inPoint=hslider("input point", floor(nPoints/2),0,nPoints-1,1);
outPoint=hslider("output point",floor(nPoints/2),0,nPoints-1,1);

//----------------------------------Build Model-------------------------------//
model =
    (route(2*nPoints,nPoints*nSidePoints+nPoints,
        par(i,nPoints,
            par(j, spaceDep,
                (i+1, nPointInputs*i + nSidePoints*j - nPointInputs*(spaceDep-1)),
                (i+1, nPointInputs*i + nSidePoints*j - nPointInputs*(spaceDep-1) +
                2*spaceDep^2+nPointInputs))),
        par(i,nPoints,(nPoints+1+i,1+nPointInputs*i))):
            buildStringDamp(nPoints,modelType))~
              si.bus(nPoints):
                si.bus(nPoints):ba.selectn(nPoints,outPoint);

//(si.bus(nPoints):par(i,nPoints,*select2(stop,1,0))):

process = hit:ba.selectoutn(nPoints,inPoint):model<:_,_;
//process = buildStringDamp(nPoints,modelType);
