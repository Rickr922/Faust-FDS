import("stdfaust.lib");

//--------------------------------Model Settings-----------------------------//
k = 1/ma.SR;
c = 344;
h = c * k*sqrt(2);
nPointsX = 4;
nPointsY = 3;

lambda = c*k/h;

alpha = 1;// lambda*lambda;
beta = 2; //2*(1-2*lambda*lambda);

midCoeff = 0,alpha,0,
           alpha,beta,alpha,
           0,alpha,0,
           0,0,0,
            0,-1,0,
            0,0,0;

r=1;
t=1;

scheme(pointsX,pointsY) = par (i, pointsX,
                                par(j,pointsY, midCoeff));

schemePoint2D(R,T,coeff,fIn) = si.bus(nNeighbors)<:neighbors:
    sum(t,T+1,
        sum(i,nNeighbors,
            ba.take(i+1+int(t*nNeighbors),coeff)*ba.take(i+1,neighbors)@(t)))
                + fIn
    with
    {
        nNeighbors = (2*R+1)^2;
        neighbors = si.bus(nNeighbors*(T+1));
    };

buildScheme2D(R,T,pointsX,pointsY,coefficients) =
    par (x, pointsX,
        par(y,pointsY, schemePoint2D(R,T,par(i,coeffsLength,coeffs(x,y,i)))))
    with
    {
        nPoints = pointsX*pointsY;
        nNeighbors = (2*R+1)^2;
        //coeff(x,y) = ba.subseq(coeffs,int((x*pointsY+y)*coeffLength),coeffLength);
        coeffsLength = int(nNeighbors*(T+1));
        coeffs(x,y,i) = ba.selector((x*pointsY+y)*coeffsLength+i,coeffsLength*nPoints,coefficients);
    };

//process = 10,par(i,(2*r+1)^2,i):schemePoint2D(r,t,midCoeff);
process = par(i,120,i):buildScheme2D(r,t,nPointsX,nPointsY,scheme(nPointsX,nPointsY));
//process = ba.take(20,scheme(nPointsX,nPointsY));
//process=ba.count(scheme(nPointsX,nPointsY));
//process=takeFromCoeff(1,scheme(nPointsX,nPointsY));
//process = scheme(nPointsX,nPointsY);
