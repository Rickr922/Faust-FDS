import("stdfaust.lib");

/*Number of points: 19*/

k = 1/ma.SR;
c = 344;
h = c*k;
s0 = 1000;

nPoints = 19;

C1 = 1+(2*s0*k*k/h);
C2 = (2*s0*k*k/h)-1;

lambda = c*k/h;

//------------------------------------------Equations--------------------------
nSpatialDependency = 2; //n° of spatial side points needed by the update eq

stringDampMidPoint(fIn, uSide_l, uSide_r) = u //u_(l-1)=uSide_l, u_(l+1)=uSide_
  letrec
  {
    'u = u*2*(1-lambda*lambda)/C1+u'*C2/C1+ lambda*lambda*(uSide_l+uSide_r)/C1 + fIn;
  }<:par(i,nSpatialDependency,_);

stringDampFreePoint(fIn,uSide) = u
  letrec
  {
    'u = u*2*(1-lambda*lambda)/C1+u'*C2/C1+ lambda*lambda*2*uSide/C1 + fIn;
  }<:par(i,nSpatialDependency-1,_);

fixedPoint(fIn,uSide) = par(i,nSpatialDependency-1,0);

//----------------------------------------Controls-----------------------------
hit = button("hit"):ba.impulsify;
inPoint=hslider("input point", 9,0,nPoints-1,1);
outPoint=hslider("output point",5,0,(nPoints*nSpatialDependency-2)/nSpatialDependency,1);
                                    //because every module has 2 outs, except
                                    //the boundaries-> -2

//-----------------------------------------Model------------------------------

model =
    (route((nPoints*nSpatialDependency-2+nPoints),(nPoints*nSpatialDependency-2+nPoints),
        //Feedback connections
              (1,4),
        (2,2),(3,7),
        (4,5),(5,10),
        (6,8),(7,13),
        (8,11),(9,16),
        (10,14),(11,19),
        (12,17),(13,22),
        (14,20),(15,25),
        (16,23),(17,28),
        (18,26),(19,31),
        (20,29),(21,34),
        (22,32),(23,37),
        (24,35),(25,40),
        (26,38),(27,43),
        (28,41),(29,46),
        (30,44),(31,49),
        (32,47),(33,52),
        (34,50),(35,55),
        (36,53),
        //Force connections
        (37,1),(38,3),(39,6),(40,9),(41,12),(42,15),(43,18),(44,21),(45,24),(46,27),(47,30),(48,33),(49,36),(50,39),(51,42),(52,45),(53,48),(54,51),(55,54)):
                stringDampFreePoint,
                stringDampMidPoint,
                stringDampMidPoint,
                stringDampMidPoint,
                stringDampMidPoint,
                stringDampMidPoint,
                stringDampMidPoint,
                stringDampMidPoint,
                stringDampMidPoint,
                stringDampMidPoint,
                stringDampMidPoint,
                stringDampMidPoint,
                stringDampMidPoint,
                stringDampMidPoint,
                stringDampMidPoint,
                stringDampMidPoint,
                stringDampMidPoint,
                stringDampMidPoint,
                fixedPoint)~par(i, (nPoints*nSpatialDependency-2), _):
                _,par(i,nPoints-2,!,_),_:ba.selectn(nPoints,outPoint);

process = hit:ba.selectoutn(nPoints,inPoint):model<:_,_;
