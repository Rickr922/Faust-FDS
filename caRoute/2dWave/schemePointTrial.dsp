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

r=1;
d=2;

scheme(nPointsX,nPointsY) = par (nPointsX, X,
                par(nPointsY,Y, midCoeff));

schemePoint2D(R,(coeffs),fIn,(connections)) = //ba.take(5,coeffs)*ba.take(5,connections)
sum(i,nNeighbors,ba.take(i+1,coeffs)*ba.take(i+1,connections)) + fIn
    with
    {
        nNeighbors = (2*R+1)^2;
        connections = si.bus(nNeighbors);
    };

/*
buildScheme2D(pointsX,pointsY,(coeffs)) = par (x, pointsX,
                            par(x,pointsY, schemePoint))
    with
    {
        schemePoint()
    }*/

    process = _<:schemePoint2D(r,midCoeff,10);
