import("stdfaust.lib");

//--------------------------------Model Settings-----------------------------//
k = 1/ma.SR;
c = 344;
h = c * k*sqrt(2);
nPointsX = 4;
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

//WITH selector
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

//process = 10,par(i,(2*r+1)^2,i):schemePoint2D(r,t,coefficients);
process = par(i,120,i):buildScheme2D(r,t,nPointsX,nPointsY,scheme(nPointsX,nPointsY));
