import("stdfaust.lib");

//--------------------------------Model Settings-----------------------------//
k = 1/ma.SR;
c = 344;
h = c * k*sqrt(2);
nPointsX = 2;
nPointsY = 2;

lambda = c*k/h;

alpha = 1;// lambda*lambda;
beta = 2; //2*(1-2*lambda*lambda);

midCoeff = 0,alpha,0,
           alpha,beta,alpha,
           0,alpha,0,
           0,0,0,
            0,-1,0,
            0,0,0;

/*midCoeff = 0,alpha,0,
           alpha,beta,alpha,
           0,alpha,0;

midCoeffDelay1 = 0,0,0,
                0,-1,0,
                0,0,0;*/

r=1;
t=1;

coefficients = midCoeff,midCoeffDelay1;

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

//RICORDATI CHE PUOI USARE subseq PER ESTRARRE ELEMENTI DA UNA LISTA

buildScheme2D(R,T,pointsX,pointsY,coeffs) = /*connections:*/
    par (x, pointsX,
        par(y,pointsY, schemePoint2D(R,T,par(i,coeffLength,coeff(x,y,i)))))
    with
    {
        nPoints = pointsX*pointsY;
        nNeighbors = (2*R+1)^2;
        connections = si.bus(nNeighbors*nPoints+nPoints);
        //coeff(x,y) = ba.subseq(coeffs,int((x*pointsY+y)*coeffLength),coeffLength);
        coeff(x,y,i) = ba.selector((x*pointsY+y)*coeffLength+i,coeffLength*nPoints,coeffs);
        coeffLength = int(nNeighbors*(T+1));
    };

//process = 10,par(i,(2*r+1)^2,i):schemePoint2D(r,t,midCoeff);
process = par(i,40,i):buildScheme2D(r,t,2,2,scheme(2,2));
//process = ba.take(20,scheme(nPointsX,nPointsY));
//process=ba.count(scheme(nPointsX,nPointsY));
//process=takeFromCoeff(1,scheme(nPointsX,nPointsY));
//process = scheme(nPointsX,nPointsY);
