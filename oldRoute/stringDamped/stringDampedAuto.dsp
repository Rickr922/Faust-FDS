import("stdfaust.lib");

//--------------------------------Model Settings-----------------------------//
nPoints = 4;

//modelType: 1->free-fixed; 2->fixed-free; else->fixed-fixed
modelType = 0;

k = 1/ma.SR;
c = 344;
h = c*k;
s0 = 1000;

//----------------------------------Equations--------------------------------//
spaceDep = 2; //n° of spatial side points needed by the update eq

C1 = 1+(2*s0*k*k/h);
C2 = (2*s0*k*k/h)-1;

lambda = c*k/h;

stringDampMidPoint(fIn, uSide_l, uSide_r) = u //u_(l-1)=uSide_l,u_(l+1)=uSide_r
  letrec
  {
    'u = u*2*(1-lambda*lambda)/C1+u'*C2/C1+ lambda*lambda*(uSide_l+uSide_r)/C1 + fIn;
  }<:par(i,spaceDep,_);

stringDampFreePoint(fIn,uSide) = u
  letrec
  {
    'u = u*2*(1-lambda*lambda)/C1+u'*C2/C1+ lambda*lambda*2*uSide/C1 + fIn;
  }<:par(i,spaceDep-1,_);

fixedPoint(fIn,uSide) = par(i,spaceDep-1,0);

buildStringDamp(nPoints,modelType) =
  //modelType: 1->free-fixed; 2->fixed-free; else->fixed-fixed
  (_,_<:(fixedPoint, stringDampFreePoint:select2(modelType==1))),
  par(i, nPoints-2, stringDampMidPoint),
  (_,_<:(fixedPoint, stringDampFreePoint:select2(modelType==2)));

//----------------------------------Controls---------------------------------//
hit = button("hit"):ba.impulsify;
inPoint=hslider("input point", 1,0,nPoints-1,1);
outPoint=hslider("output point",1,0,(nPoints*spaceDep-2)/spaceDep,1);
                                    //because every module has 2 outs, except
                                    //the boundaries-> -2

//----------------------------------Build Model-------------------------------//
nConnections = nPoints*spaceDep-2+nPoints;
K = (nPoints*spaceDep-4)/2;

model =
    (route(nConnections,nConnections,
        par(i,K+1,(1+2*i,4+3*i),(2+2*i,2+3*i)),
        (2+2*K+1,1),
        par(i,nPoints-1,(2+2*K+(i+2),3*(i+1)))):
            buildStringDamp(nPoints,modelType))~par(i, (nPoints*spaceDep-2), _):
                _,par(i,nPoints-2,!,_),_:ba.selectn(nPoints,outPoint);

process = hit:ba.selectoutn(nPoints,inPoint):model<:_,_;
