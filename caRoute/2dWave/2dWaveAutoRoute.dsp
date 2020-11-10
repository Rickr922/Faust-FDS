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
inPointX=hslider("input point x", floor(nPointsX/2),0,nPointsX-1,0.01);
inPointY=hslider("input point y", floor(nPointsY/2),0,nPointsY-1,0.01);
outPointX=hslider("output point x",floor(nPointsX/2),0,nPointsX-1,0.01);
outPointY=hslider("output point y",floor(nPointsY/2),0,nPointsY-1,0.01);
forceModel = button("play") : ba.impulsify;
stop = button("Stop");

//----------------------------------Library---------------------------------//
schemePoint2D(R,T,coeffs,fIn) = neighbors<:
    sum(t,T+1,
        sum(i,nNeighbors,
            ba.selector(int(i+t*nNeighbors),nNeighbors*(T+1),coeffs)*
            ba.selector(i,nNeighbors,neighbors)@(t)))
                + fIn
with
{
    nNeighbors = (2*R+1)^2;
    neighbors = si.bus(nNeighbors);
};

buildScheme2D(R,T,pointsX,pointsY,coefficients) =
    par (x, pointsX,
        par(y,pointsY, schemePoint2D(R,T,par(i,coeffsLength,coeffs(x,y,i)))))
with
{
    nPoints = pointsX*pointsY;
    nNeighbors = (2*R+1)^2;
    coeffsLength = int(nNeighbors*(T+1));
    coeffs(x,y,i) = ba.selector((x*pointsY+y)*coeffsLength+i,coeffsLength*nPoints,coefficients);
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
route2D(X, Y, r) = route(X*Y*2, X*Y*nInputs,
                                par(x, X, par(y, Y, connections(x,y))))
with
{
    connections(x,y) = P(x,y) + X*Y, C(x,y,0),
                       par(j,nNeighbors,
                         par(i,nNeighbors,
                           P(x,y),C(x+i-r,y+j-r,nNeighbors^2-(i*nNeighbors+j))));

    P(x,y) = x*Y+y+1;
    C(x,y,count) = (1 + count + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
    nNeighbors = 2*r+1;
    nInputs = nNeighbors^2+1;
};

model(X,Y,r,t,coeffs) = (route2D(X,Y,r) : buildScheme2D(r,t,X,Y,coeffs)) ~ par(i,X*Y,_*(stop==0));


process = forceModel <: linInterpolation2D(inPointX,inPointY) : model(nPointsX,nPointsY,r,t,coefficients) : linInterpolation2DOut(outPointX,outPointY)<:_,_;
