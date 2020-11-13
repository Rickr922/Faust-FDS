import("stdfaust.lib");

//--------------------------------Model Settings-----------------------------//
k = 1/ma.SR;
K = 20;
s0 = 2;
s1 = 0.05;
c=344;

coeff = c^2*k^2+4*s1*k;
h = sqrt((coeff+sqrt(coeff*coeff+16*K*K*k*k)));

nPointsX = 5;
nPointsY = 5;

lambda = c*k/h;

//----------------------------------Equations--------------------------------//
mu=K*K*k*k/(h^4);
den = 1+s0*k;
A = 2*(1-10*mu-2*lambda*lambda-4*s1*k*k)/den;
B = (s0*k+4*k*k-1)/den;
C = (8*mu + lambda*lambda + 2*s1*k*k)/den;
D = -2*mu/den;
E = -mu/den;
F = -2*s1*k*k/den;

midCoeff = 0,0,E,0,0,
           0,D,C,D,0,
           E,C,A,C,E,
           0,D,C,D,0,
           0,0,E,0,0;

midCoeffDelay1 = 0,0,0,0,0,
                 0,0,F,0,0,
                 0,F,B,F,0,
                 0,0,F,0,0,
                 0,0,0,0,0;

r=2;
t=1;

coefficients = midCoeff,midCoeffDelay1;

scheme(pointsX,pointsY) = par (i, pointsX,
                                par(j,pointsY, coefficients));

//----------------------------------Controls---------------------------------//
inPointX=hslider("input point x", floor(nPointsX/2),0,nPointsX-1,0.01);
inPointY=hslider("input point y", floor(nPointsY/2),0,nPointsY-1,0.01);
outPointX=hslider("output point x",floor(nPointsX/2),0,nPointsX-1,0.01);
outPointY=hslider("output point y",floor(nPointsY/2),0,nPointsY-1,0.01);
hit = button("play");

//----------------------------------Library---------------------------------//
schemePoint2D(R,T) = routing:operations:>_
with
{
    nNeighbors = (2*R+1)^2;
    routing =
        route(nNeighbors*(T+1)+nNeighbors+1,2*nNeighbors*(T+1)+1,
            (1,1),
            par(t,T+1,
                par(i,nNeighbors,i+t*nNeighbors+2,2*(i+t*nNeighbors)+3,
                                i+nNeighbors*(T+1)+2,2*(i+t*nNeighbors)+2)));
    operations = _,par(t,T+1,
                    par(i,nNeighbors,(_@t),_:*));
};

buildScheme2D(pointsX,pointsY,R,T) =
    par (x, pointsX,
        par(y,pointsY, schemePoint2D(R,T)));

//----------------------------------Interpolation---------------------------------//
linInterpolation2D(X,Y,pointX,pointY) =
par(i,X,
    par(j,Y,_*
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
forceModel = hit:ba.impulsify;
stairsForce(X,Y,pointX,pointY) = ba.selectoutn(X*Y,pointY+pointX*Y);

//----------------------------------Output-------------------------------//
stairsOutput(X,Y,pointX,pointY) = ba.selectn(X*Y,pointY+pointX*Y);
linInterpolation2DOut(X,Y,pointX,pointY) = linInterpolation2D(X,Y,pointX,pointY):>_;
//----------------------------------Build Model-------------------------------//
route2D(X, Y, R, T) = route(nPoints*2+nPoints*nCoeffs, nPoints*nInputs,
                                par(x, X, par(y, Y, connections(x,y))))
with
{
    connections(x,y) =  P(x,y) + nPoints, C(x,y,0),
                        par(k,nCoeffs,(x*Y+y)*nCoeffs+k+1,C(x,y,k+1)),
                        par(j,nNeighborsXY,
                            par(i,nNeighborsXY,
                                P(x,y),C(x+i-R,y+j-R,nInputs-1-(i*nNeighborsXY+j))));

    P(x,y) = x*Y+y+1 + nCoeffs*nPoints;
    C(x,y,count) = (1 + count + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);

    nNeighborsXY = 2*R+1;
    nNeighbors = nNeighborsXY^2;
    nCoeffs = nNeighbors*(T+1);
    nInputs = nNeighbors+1+nCoeffs;
    nPoints = X*Y;
};

model2D(X,Y,r,t,scheme) =
    (route2D(X,Y,r,t,scheme) : buildScheme2D(X,Y,r,t)) ~ si.bus(X*Y);


process =
    forceModel<:linInterpolation2D(nPointsX,nPointsY,inPointX,inPointY):
        model2D(nPointsX,nPointsY,r,t,scheme(nPointsX,nPointsY)):
            linInterpolation2DOut(nPointsX,nPointsY,outPointX,outPointY);
