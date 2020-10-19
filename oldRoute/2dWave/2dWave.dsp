import("stdfaust.lib");

//--------------------------------Model Settings-----------------------------//
nPointsX = 30;
nPointsY = 10;

k = 1/ma.SR;
c = 344;
h = c * k*sqrt(2);
s0 = 500;

//----------------------------------Equations--------------------------------//
nPointInputs = inputs(waveMidPoint);
nDepX = 2;
nDepY = 2;
nDep = nPointInputs-2;
spaceDep = (nPointInputs-2)/2;    //n° of spatial side points needed by the update eq

lambda = c*k/h;

waveMidPoint(i,j, uSide_w, uSide_e, uSide_s, uSide_n) = equation
with
{
    //fIn = forceModel*(i==inPointX)*(j==inPointY); stairs approx
    fIn = linInterpolation2dIn(i,j,inPointX,inPointY,forceModel);
    equation = u ///u_(l-1)=uSide_w,u_(l+1)=uSide_e,u_(m-1)=uSide_s,u_(m+1)=uSide_n,
    letrec
    {
        'u = u*2*(1-2*lambda*lambda) - u' + lambda*lambda*(uSide_e+uSide_w+uSide_n+uSide_s) + fIn;
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

buildScheme(pointsX,pointsY) =
  par(i, pointsX,
    par(j,pointsY, waveMidPoint(i,j)));
  //Nota che non mi servono i fixed points perchè nelle connessioni
  // vuote entra 0!! Praticamente è come se avessi dei virtual boundaries

//----------------------------------Controls---------------------------------//
hit = button("hit"):ba.impulsify;
stop = button("Stop");
inPointX=hslider("input point x", floor(nPointsX/2),0,nPointsX-1,0.01);
inPointY=hslider("input point y", floor(nPointsY/2),0,nPointsY-1,0.01);
outPointX=hslider("output point x",floor(nPointsX/2),0,nPointsX-1,0.01);
outPointY=hslider("output point y",floor(nPointsY/2),0,nPointsY-1,0.01);

//----------------------------------Force---------------------------------//
forceModel = hit;
linInterpolation2dIn(i,j,pointX,pointY,force) =
force*select2((i==intX) & (j==intY),
        select2((i==(intX+1)) & (j==intY),
            select2((i==intX) & (j==(intY+1)),
                select2((i==(intX+1)) & (j==(intY+1)),
                    0,
                    fractionX*fractionY),
                (1-fractionX)*fractionY),
            fractionX*(1-fractionY)),
        (1-fractionX)*(1-fractionY))
with
{
    fractionX = ma.frac(pointX);
    fractionY = ma.frac(pointY);
    intX = int(pointX);
    intY = int(pointY);
};

//----------------------------------Out Points-------------------------------//
linInterpolation2dOut(pointX,pointY) =
par(i,nPointsX,
    par(j,nPointsY,_*
        select2((i==intX) & (j==intY),
            select2((i==(intX+1)) & (j==intY),
                select2((i==intX) & (j==(intY+1)),
                    select2((i==(intX+1)) & (j==(intY+1)),
                        0,
                        fractionX*fractionY),
                    (1-fractionX)*fractionY),
                fractionX*(1-fractionY)),
            (1-fractionX)*(1-fractionY)))):>_
with
{
    fractionX = ma.frac(pointX);
    fractionY = ma.frac(pointY);
    intX = int(pointX);
    intY = int(pointY);
};

//----------------------------------Build Model-------------------------------//
model =
(route(nPointsX*nPointsY,nPointsX*nPointsY*nDep,
    par(i,nPointsX,
        par(j, nPointsY,
            (j+1+nPointsY*i,
              nDep*(j+nPointsY*i)*(fmod(j+nPointsY*i,nPointsY)!=0)),
            (j+1+nPointsY*i,
              (nDep*(j+nPointsY*i)+nDepX*2+nDepY+1)*(fmod((j+1)+nPointsY*i,nPointsY)!=0)))),
    par(i,nPointsX,
        par(j, nPointsY,
            (j+1+nPointsY*i,4*(j+nPointsY*i)+2-nDep*nPointsY),
            (j+1+nPointsY*i,4*(j+nPointsY*i)+nDep*nPointsY+1)))):
                buildScheme(nPointsX,nPointsY))~
                    si.bus(nPointsX*nPointsY):
                        si.bus(nPointsX*nPointsY):
                            linInterpolation2dOut(outPointX,outPointY);

process = model<:_,_;
//process = fmod(0,3)!=0;
//process=_<:par(i,3,par(j,3,_*(i*j)));
