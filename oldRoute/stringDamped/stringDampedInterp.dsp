import("stdfaust.lib");

//--------------------------------Model Settings-----------------------------//
nPoints = 10;

//modelType: 1->free-fixed; 2->fixed-free; else->fixed-fixed
modelType = 0;

k = 1/ma.SR;
c = 344;
h = c*k;
s0 = 500;

//----------------------------------Equations--------------------------------//
nPointInputs = inputs(stringDampMidPoint);
spaceDep = (nPointInputs-1)/2;    //n° of spatial side points needed by the update eq

C1 = 1+(2*s0*k*k/h);
C2 = (2*s0*k*k/h)-1;

lambda = c*k/h;

stringDampMidPoint(i, uSide_l, uSide_r) = equation
with
{
    fIn = linInterpolForce(i,inPoint,forceModel);
    equation = u //u_(l-1)=uSide_l,u_(l+1)=uSide_r
    letrec
    {
        'u = u*2*(1-lambda*lambda)/C1+u'*C2/C1+ lambda*lambda*(uSide_l+uSide_r)/C1 + fIn;
    };
};

stringDampFreePointL(i, uSide_l, uSide_r) = equation
with
{
    fIn = linInterpolForce(i,inPoint,forceModel);
    equation = u
    letrec //u_(l-1)=uSide_l,u_(l+1)=uSide_r
    {
        'u = u*2*(1-lambda*lambda)/C1+u'*C2/C1+ lambda*lambda*2*uSide_r/C1 + fIn;
    };
};

stringDampFreePointR(i, uSide_l, uSide_r) = equation
with
{
    fIn = linInterpolForce(i,inPoint,forceModel);
    equation = u
    letrec     //u_(l-1)=uSide_l,u_(l+1)=uSide_r
    {
        'u = u*2*(1-lambda*lambda)/C1+u'*C2/C1+ lambda*lambda*2*uSide_l/C1 + fIn;
    };
};

fixedPoint(i, uSide_l, uSide_r) = 0;

buildStringDamp(nPoints,modelType) = par(i, nPoints, stringDampMidPoint(i));

/*
  //modelType: 1->free-fixed; 2->fixed-free; else->fixed-fixed
  (si.bus(2*spaceDep)<:(stringDampMidPoint(0), stringDampFreePointL(0):select2(modelType==1))),
  par(i, nPoints-2, stringDampMidPoint(i+1)),
  (si.bus(2*spaceDep)<:(stringDampMidPoint(nPoints-1), stringDampFreePointR(nPoints-1):select2(modelType==2)));*/
  //Nota che non mi servono i fixed points perchè nelle connessioni
  // vuote entra 0!! Praticamente è come se avessi dei virtual boundaries

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
