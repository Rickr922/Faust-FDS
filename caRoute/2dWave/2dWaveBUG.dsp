import("stdfaust.lib");

//--------------------------------Model Settings-----------------------------//
k = 1/ma.SR;
c = 344;
h = c * k*sqrt(2);
nPointsX = 3;
nPointsY = 3;

lambda = c*k/h;

alpha = lambda*lambda;
beta = 2*(1-2*lambda*lambda);

midCoeff = 0,alpha,0,
           alpha,beta,alpha,
           0,alpha,0;

midCoeffDelay1 = 0,0,0,
                0,-1,0,
                0,0,0;

r=1;
t=1;

coefficients = midCoeff,midCoeffDelay1;

scheme(pointsX,pointsY) = par (i, pointsX,
                                par(j,pointsY, coefficients));

//----------------------------------Controls---------------------------------//
inPointX=hslider("input point x", floor(nPointsX/2),0,nPointsX-1,1);
inPointY=hslider("input point y", floor(nPointsY/2),0,nPointsY-1,1);
outPointX=hslider("output point x",floor(nPointsX/2),0,nPointsX-1,1);
outPointY=hslider("output point y",floor(nPointsY/2),0,nPointsY-1,1);
forceModel = button("play") : ba.impulsify;
stop = button("Stop");

//----------------------------------Library---------------------------------//
schemePoint2D(R,T,fIn) = coeffs,neighbors<:
    sum(t,T+1,
        sum(i,nNeighbors,
            ba.selector(int(i+t*nNeighbors),nNeighbors*(T+1),coeffs)*
            ba.selector(i,nNeighbors,neighbors)@(t)))
                + fIn
with
{
    nNeighbors = (2*R+1)^2;
    neighbors = si.bus(nNeighbors);
    coeffs=si.bus(nNeighbors*(T+1));
};

buildScheme2D(R,T,pointsX,pointsY) =
    par (x, pointsX,
        par(y,pointsY, schemePoint2D(R,T)))
with
{/*
    nPoints = pointsX*pointsY;
    nNeighbors = (2*R+1)^2;
    coeffsLength = int(nNeighbors*(T+1));
    coeffs(x,y,i) = ba.selector((x*pointsY+y)*coeffsLength+i,coeffsLength*nPoints,coefficients);*/
};

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
//nInputs = inputs(schemeMidPoint);
route2D(X, Y, R, T) = route(nPoints*2+nPoints*nCoeffs, nPoints*nInputs,
                                par(x, X, par(y, Y, connections(x,y))))
with
{
    connections(x,y) =  par(k,nCoeffs,nPoints*2+(x*Y+y)*nCoeffs+k,C(x,y,k)),
                        P(x,y) + nPoints + nCoeffs*nPoints, C(x,y,nCoeffs),
                        par(j,nNeighborsXY,
                         par(i,nNeighborsXY,
                           P(x,y),C(x+i-R,y+j-R,nInputs-(i*nNeighborsXY+j))));

    P(x,y) = x*Y+y+1;
    C(x,y,count) = (1 + count + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);

    nNeighborsXY = 2*R+1;
    nNeighbors = nNeighborsXY^2;
    nCoeffs = nNeighbors*(T+1);
    nInputs = nNeighbors+1+nCoeffs;
    nPoints = X*Y;
};

model(X,Y,r,t,coeffs) =
    (route2D(X,Y,r,t) : buildScheme2D(r,t,X,Y)) ~ par(i,X*Y,_/**(stop==0)*/);


process = forceModel <: model(nPointsX,nPointsY,r,t,scheme(nPointsX,nPointsY));
//process = coefficients,10,par(i,9,i):schemePoint2D(1,1);
