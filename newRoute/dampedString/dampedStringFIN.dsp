import("stdfaust.lib");

//--------------------------------Model Settings-----------------------------//
nPoints = 60;

k = 1/ma.SR;
c = 344;
h = c*k;
s0 = 500;

C1 = 1+(2*s0*k*k/h);
C2 = (2*s0*k*k/h)-1;

//----------------------------------Equations--------------------------------//
lambda = c*k/h;
nInputs = inputs(schemeMidPoint);

//west=left, east=right
schemeMidPoint(fIn,u_w,u,u_e) =
        u*2*(1-lambda*lambda)/C1+u'*C2/C1+ lambda*lambda*(u_w+u_e)/C1 + fIn;

schemeFreePointWest(fIn,u_w,u,u_e) =
        u*2*(1-lambda*lambda)/C1+u'*C2/C1+ lambda*lambda*2*u_e/C1 + fIn;

schemeFreePointEast(fIn,u_w,u,u_e) =
        u*2*(1-lambda*lambda)/C1+u'*C2/C1+ lambda*lambda*2*u_w/C1 + fIn;

schemeFixedPoint(fIn,u_w,u,u_e) = 0;

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

linInterp1D(outPoint,nPoints) = par(i,nPoints,_*select2(
        i==int(outPoint), select2(i==int(outPoint+1),0,fraction),(1-fraction)))
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
with {
    connections(x) = Current(x)+nPoints, F(x),
                     Current(x), E(x-1),
                     Current(x), C(x),
                     Current(x), W(x+1);
    Current(x) = x+1;
    F(x) = (1 + 0 + (x*nInputs)) * (x>=0) * (x<nPoints);
    W(x) = (1 + 1 + (x*nInputs)) * (x>=0) * (x<nPoints);
    C(x) = (1 + 2 + (x*nInputs)) * (x>=0) * (x<nPoints);
    E(x) = (1 + 3 + (x*nInputs)) * (x>=0) * (x<nPoints);
};

model(nPoints) = (routing(nPoints,nInputs) : buildScheme(nPoints)) ~ si.bus(nPoints);

process = forceModel<:linInterp1D(inPoint,nPoints): model(nPoints):linInterp1DOut(outPoint,nPoints)<:_,_;
