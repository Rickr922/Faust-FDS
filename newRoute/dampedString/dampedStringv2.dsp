import("stdfaust.lib");

//--------------------------------Model Settings-----------------------------//
nPoints = 5;

k = 1/ma.SR;
c = 344;
h = c*k;
s0 = 500;

C1 = 1+(2*s0*k*k/h);
C2 = (2*s0*k*k/h)-1;

//----------------------------------Equations--------------------------------//
lambda = c*k/h;
nInputs = inputs(schemeMidPoint);

schemeMidPoint(u_w,u,u_e) = u*2*(1-lambda*lambda)/C1+u'*C2/C1+ lambda*lambda*(u_w+u_e)/C1<:_,_,_;

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

//----------------------------------Output-------------------------------//
linInterp(outPoint) = par(i,nPoints,_*select2(
        i==int(outPoint), select2(i==int(outPoint+1),0,fraction),(1-fraction)))
with
{
    fraction = ma.frac(outPoint);
};
linInterpOut(outPoint) = par(i,nPoints,_*select2(
        i==int(outPoint), select2(i==int(outPoint+1),0,fraction),(1-fraction))):>_
with
{
    fraction = ma.frac(outPoint);
};

//----------------------------------Build Model-------------------------------//
buildScheme(nPoints, node) = par (i, nPoints, node);

routing(nPoints,nInputs) =   route(nPoints*nInputs+nPoints, nPoints*nInputs,
                        par(x, nPoints, connections(x)),
                            par(i, nPoints, nPoints*nInputs + 1 + i, C(i)))
with {
    connections(x) =  W(x), E(x-1),
                        C(x), C(x),
                        E(x), W(x+1);

    W(x) = (1 + 0 + (x*nInputs)) * (x>=0) * (x<nPoints);
    C(x) = (1 + 1 + (x*nInputs)) * (x>=0) * (x<nPoints);
    E(x) = (1 + 2 + (x*nInputs)) * (x>=0) * (x<nPoints);
};

model(nPoints) = (routing(nPoints,nInputs) : buildScheme(nPoints, schemeMidPoint)) ~ si.bus(nPoints*nInputs) : par(i,nPoints,!,_,!);

process = forceModel<:linInterp(inPoint): model(nPoints):linInterpOut(outPoint)<:_,_;
