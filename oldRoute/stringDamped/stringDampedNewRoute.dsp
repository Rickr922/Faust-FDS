import("stdfaust.lib");

//--------------------------------Model Settings-----------------------------//
nPoints = 150;

//modelType: 1->free-fixed; 2->fixed-free; else->fixed-fixed
modelType = 0;

k = 1/ma.SR;
c = 344;
//Stability condition
h = c*k;
s0 = 1000;

//----------------------------------Equations--------------------------------//
nPointInputs = inputs(stringDampMidPoint);
nSidePoints = nPointInputs-1; //twice the number of side points (left + right)
spaceDep = nSidePoints/2;    //n° of spatial side points needed by the update eq

C1 = 1+(2*s0*k*k/h);
C2 = (2*s0*k*k/h)-1;

lambda = c*k/h;

stringDampMidPoint(fIn, uSide_l, uSide_r) = u //u_(l-1)=uSide_l,u_(l+1)=uSide_r
  letrec
  {
    'u = u*2*(1-lambda*lambda)/C1+u'*C2/C1+ lambda*lambda*(uSide_l+uSide_r)/C1 + fIn;
  };

stringDampFreePointL(fIn, uSide_l, uSide_r) = u
  letrec
  {
    'u = u*2*(1-lambda*lambda)/C1+u'*C2/C1+ lambda*lambda*2*uSide_r/C1 + fIn;
  };

stringDampFreePointR(fIn, uSide_l, uSide_r) = u
  letrec
  {
    'u = u*2*(1-lambda*lambda)/C1+u'*C2/C1+ lambda*lambda*2*uSide_l/C1 + fIn;
  };

fixedPoint(fIn, uSide_l, uSide_r) = 0;

buildStringDamp(nPoints,modelType) =
  //modelType: 1->free-fixed; 2->fixed-free; else->fixed-fixed
  (si.bus(nPointInputs)<:(stringDampMidPoint, stringDampFreePointL:select2(modelType==1))),
  par(i, nPoints-2, stringDampMidPoint),
  (si.bus(nPointInputs)<:(stringDampMidPoint, stringDampFreePointR:select2(modelType==2)));
  //Nota che non mi servono i fixed points perchè nelle connessioni
  // vuote entra 0!! Praticamente è come se avessi dei virtual boundaries

//----------------------------------Controls---------------------------------//
hit = button("hit"):ba.impulsify;
inPoint=hslider("input point", 1,0,nPoints-1,1);
outPoint=hslider("output point",1,0,nPoints-1,1);

//----------------------------------Build Model-------------------------------//
model =
    (route(2*nPoints,nPoints*nSidePoints+nPoints,
        par(i,nPoints,
            par(j, spaceDep,
                (i+1, nPointInputs*i + nSidePoints*j - nPointInputs*(spaceDep-1)),
                (i+1, nPointInputs*i + nSidePoints*j - nPointInputs*(spaceDep-1) +
                2*spaceDep^2+nPointInputs))),
        par(i,nPoints,(nPoints+1+i,1+nPointInputs*i))):
            buildStringDamp(nPoints,modelType))~si.bus(nPoints):
              si.bus(nPoints):ba.selectn(nPoints,outPoint);

process = hit:ba.selectoutn(nPoints,inPoint):model<:_,_;
//process = buildStringDamp(nPoints,modelType);
