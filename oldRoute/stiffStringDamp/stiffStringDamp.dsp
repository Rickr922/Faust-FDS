import("stdfaust.lib");

//--------------------------------Model Settings-----------------------------//
nPoints = 6;

//modelType: 1->free-fixed; 2->fixed-free; else->fixed-fixed
modelType = 0;

k = 1/ma.SR;
coeff = c^2*k^2 + 4*sigma1*k;
h =sqrt((coeff + sqrt((coeff)^2 + 16*k^2*K^2))/2);

T = 75;                 // Tension [N]
radius = 3.5560e-04;    // Radius (0.016 gauge) [m]
rho = 8.05*10^3;        // Density [kg/m^3];
Area = ma.PI*radius^2;        // Area of string section
I = (ma.PI*radius^4)/ 4;   // Moment of Inertia
L = 1;                  // String length [m]
Emod = 174e6;              // Young modulus [Pa]
K = sqrt(Emod*I/rho/Area);    // Stiffness parameter
c = sqrt(T/rho/Area);      // Wave speed
sigma1 = 0.005;         // Frequency dependent damping
sigma0 = 0.5;

//----------------------------------Equations--------------------------------//
nPointInputs = inputs(midPointEq);
spaceDep = nPointInputs-1; //n° of spatial side points needed by the update eq
sidePoints = spaceDep/2;

den = 1+sigma0*k;
A = (2*h^4-2*c^2*k^2*h^2-4*sigma1*k*h^2+6*K^2*k^2)/den/h^4;
B = (sigma0*k*h^2-h^2+4*sigma1*k)/den/h^2;
C = (c^2*k^2*h^2+2*sigma1*k*h^2-4*K^2*k^2)/den/h^4;
D = -2*sigma1*k/den/h^2;
E = K^2*k^2/den/h^4;

lambda = c*k/h;

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
  //Nota che non mi servono i fixed points perchè nelle connessioni vuote entra 0!! (e i free points ancora non funzionano)

  //modelType: 1->free-fixed; 2->fixed-free; else->fixed-fixed
  /*par(i,sidePoints,
        (si.bus(nPointInputs)<:(fixedPoint, freePointEqL:select2(modelType==1)))),
  par(i, nPoints-spaceDep, midPointEq),
  par(i,sidePoints,
        (si.bus(nPointInputs)<:(fixedPoint, freePointEqR:select2(modelType==1))));*/

//----------------------------------Controls---------------------------------//
hit = button("hit"):ba.impulsify;
inPoint=hslider("input point", 1,0,nPoints-1,1);
outPoint=hslider("output point",1,0,nPoints-1,1);
                                    //because every module has 2 outs, except
                                    //the boundaries-> -2

//----------------------------------Build Model-------------------------------//
model =
    (route(2*nPoints,nPoints*spaceDep+nPoints,
        par(i,nPoints,(i+1,nPointInputs*i - 5),(i+1,nPointInputs*i - 1),
                      (i+1,nPointInputs*i + 8),(i+1,nPointInputs*i + 12)),
        par(i,nPoints,(nPoints+1+i,1+nPointInputs*i))):
            buildStringDamp(nPoints,modelType))~si.bus(nPoints):
              si.bus(nPoints):ba.selectn(nPoints,outPoint);

process = hit:ba.selectoutn(nPoints,inPoint):model<:_,_;
//process = buildStringDamp(nPoints,modelType);
