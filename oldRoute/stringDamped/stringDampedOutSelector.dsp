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
  }<:_,_;

stringDampFreePoint(uSide) = u
  letrec
  {
    'u = u*2*(1-lambda*lambda)/C1+u'*C2/C1+ lambda*lambda*2*uSide/C1;
  };

fixedPoint(uSide) = 0;

//----------------------------------------Controls-----------------------------

hit = button("hit"):ba.impulsify;
outPoint=hslider("output point",5,0,(nPoints*nSpatialDependency-2)/nSpatialDependency,1);
                                    //because every module has 2 outs, except
                                    //the boundaries-> -2

//-----------------------------------------Model------------------------------

model =
    (route((nPoints*nSpatialDependency-2+1),(nPoints*nSpatialDependency-2+1),
        (1,2),
        (2,1),(3,4),
        (4,3),(5,6),
        (6,5),(7,8),
        (8,7),(9,10),
        (10,9),(11,12),
        (12,11),(13,14),
        (14,13),(15,16),
        (16,15),(17,19),//out 18 left free for the force in->everything else jumps one place
        (18,17),(19,21),
        (20,20),(21,23),
        (22,22),(23,25),
        (24,24),(25,27),
        (26,26),(27,29),
        (28,28),(29,31),
        (30,30),(31,33),
        (32,32),(33,35),
        (34,34),(35,37),
        (36,36),(37,18))://force in
            stringDampFreePoint,
            stringDampMidPoint(0),
            stringDampMidPoint(0),
            stringDampMidPoint(0),
            stringDampMidPoint(0),
            stringDampMidPoint(0),
            stringDampMidPoint(0),
            stringDampMidPoint(0),
            stringDampMidPoint(0),
            stringDampMidPoint   ,
            stringDampMidPoint(0),
            stringDampMidPoint(0),
            stringDampMidPoint(0),
            stringDampMidPoint(0),
            stringDampMidPoint(0),
            stringDampMidPoint(0),
            stringDampMidPoint(0),
            stringDampMidPoint(0),
            fixedPoint)~par(i, (nPoints*nSpatialDependency-2), _):
              _,par(i,nPoints-2,!,_),_:ba.selectn(19,outPoint);

process = hit:model<:_,_;
