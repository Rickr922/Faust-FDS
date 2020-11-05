import("stdfaust.lib");

//--------------------------------Model Settings-----------------------------//
k = 1/ma.SR;
c = 344;
h = c * k*sqrt(2);
nPointsX = 10;
nPointsY = 10;

lambda = c*k/h;

alpha = lambda*lambda;
beta = 2*(1-2*lambda*lambda);

midCoeff = 0,alpha,0,
           alpha,beta,alpha,
           0,alpha,0;

r=1;
d=2;

scheme(nPointsX,nPointsY) = par (nPointsX, X,
                par(nPointsY,Y, midCoeff));


//----------------------------------Equations--------------------------------//
schemeMidPoint(fIn,u_nw,u_n,u_ne,u_w,u,u_e,u_sw,u_s,u_se) = 0*u_nw+alpha*u_n+0*u_ne+alpha*u_w+beta*u+alpha*u_e+0*u_sw+alpha*u_s+0*u_se+(-1)*u' + fIn;

buildScheme(X,Y) = par (x, X,
                par(y,Y, schemeMidPoint));

//----------------------------------Controls---------------------------------//
inPointX=hslider("input point x", floor(nPointsX/2),0,nPointsX-1,0.01);
inPointY=hslider("input point y", floor(nPointsY/2),0,nPointsY-1,0.01);
outPointX=hslider("output point x",floor(nPointsX/2),0,nPointsX-1,0.01);
outPointY=hslider("output point y",floor(nPointsY/2),0,nPointsY-1,0.01);
forceModel = button("play") : ba.impulsify;
stop = button("Stop");

//----------------------------------Interpolation---------------------------------//
linInterpolation2D(pointX,pointY) =
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
            (1-fractionX)*(1-fractionY))))
with
{
    fractionX = ma.frac(pointX);
    fractionY = ma.frac(pointY);
    intX = int(pointX);
    intY = int(pointY);
};

//----------------------------------Force---------------------------------//
stairsForce(X,Y,pointX,pointY) = ba.selectoutn(X*Y,pointY+pointX*Y);

//----------------------------------Output-------------------------------//
stairsOutput(X,Y,pointX,pointY) = ba.selectn(X*Y,pointY+pointX*Y);
linInterpolation2DOut(pointX,pointY) = linInterpolation2D(pointX,pointY):>_;
//----------------------------------Build Model-------------------------------//
nInputs = inputs(schemeMidPoint);
routing(X, Y, nInputs,r) = route(X*Y*2, X*Y*nInputs,
                                par(x, X, par(y, Y, connections(x,y))))
with
{
    connections(x,y) = P(x,y) + X*Y, C(x,y,0),
                       par(j,spaceDep,
                         par(i,spaceDep,
                           P(x,y),C(x+i-r,y+j-r,spaceDep^2-(i*spaceDep+j))));

    P(x,y) = x*Y+y+1;
    C(x,y,count) = (1 + count + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
    spaceDep = 2*r+1;
};

model(X,Y,nInputs,r) = (routing(X,Y,nInputs,r) : buildScheme(X,Y)) ~ par(i,X*Y,_*(stop==0));


process = forceModel <: linInterpolation2D(inPointX,inPointY) : model(nPointsX,nPointsY,nInputs,r) : linInterpolation2DOut(outPointX,outPointY)<:_,_;
