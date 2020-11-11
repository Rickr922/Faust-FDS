import("stdfaust.lib");

//--------------------------------Model Settings-----------------------------//
k = 1/ma.SR;
c = 344;
h = c * k*sqrt(2);
nPointsX = 3;
nPointsY = 3;

lambda = c*k/h;

r=1;

alpha = lambda*lambda;
beta = 2*(1-2*lambda*lambda);

scheme(pointsX,pointsY) = par (i, pointsX,
                                par(j,pointsY, coefficients));

//----------------------------------Controls---------------------------------//
inPointX=hslider("input point x", floor(nPointsX/2),0,nPointsX-1,0.01);
inPointY=hslider("input point y", floor(nPointsY/2),0,nPointsY-1,0.01);
outPointX=hslider("output point x",floor(nPointsX/2),0,nPointsX-1,0.01);
outPointY=hslider("output point y",floor(nPointsY/2),0,nPointsY-1,0.01);
forceModel = button("play") : ba.impulsify;
stop = button("Stop");

//----------------------------------Library---------------------------------//
schemeMidPoint(fIn,u_nw,u_n,u_ne,u_w,u,u_e,u_sw,u_s,u_se) =
    0*(u_nw@0)+alpha*(u_n@0)+0*(u_ne@0)+alpha*(u_w@0)+beta*(u@0)+alpha*(u_e@0)+0*(u_sw@0)+alpha*(u_s@0)+0*(u_se@0) +
    0*(u_nw@1)+0*(u_n@1)+0*(u_ne@1)+0*(u_w@1)+(-1)*(u@1)+0*(u_e@1)+0*(u_sw@1)+0*(u_s@1)+0*(u_se@1) +
    fIn;

buildScheme(X,Y) = par (x, X,
                par(y,Y, schemeMidPoint));
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
route2D(X, Y, r) = route(X*Y*2, X*Y*nInputs,
                                par(x, X, par(y, Y, connections(x,y))))
with
{
    connections(x,y) = P(x,y) + X*Y, C(x,y,0),
                       par(j,nNeighborsXY,
                         par(i,nNeighborsXY,
                           P(x,y),C(x+i-r,y+j-r,nNeighbors-(i*nNeighborsXY+j))));

    P(x,y) = x*Y+y+1;
    C(x,y,count) = (1 + count + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
    nNeighborsXY = 2*r+1;
    nNeighbors = nNeighborsXY^2;
    nInputs = nNeighbors+1;
};

model(X,Y,r) = (route2D(X,Y,r) : buildScheme(X,Y)) ~ par(i,X*Y,_/**(stop==0)*/);

 process = forceModel <: linInterpolation2D(inPointX,inPointY) : model(nPointsX,nPointsY,r) : linInterpolation2DOut(outPointX,outPointY)<:_,_;
