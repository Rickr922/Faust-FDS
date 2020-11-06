import("stdfaust.lib");

//--------------------------------Model Settings-----------------------------//
k = 1/ma.SR;
c = 344;
h = c * k*sqrt(2);
//nPointsX = 3;
//nPointsY = 3;

lambda = c*k/h;

alpha = 2;// lambda*lambda;
beta = 3; //2*(1-2*lambda*lambda);

midCoeff = 0,alpha,0,
           alpha,beta,alpha,
           0,alpha,0;

midCoeffDelay1 = 0,0,0,
                0,-1,0,
                0,0,0;

r=1;
d=2;
t=1;

scheme(nPointsX,nPointsY) = par (nPointsX, X,
                                par(nPointsY,Y, midCoeff,midCoeffDelay1));

schemePoint2D(R,coeffs,fIn) =
    sum(i,nNeighbors,ba.take(i+1,coeffs)*ba.take(i+1,connections)) + fIn
    with
    {
        nNeighbors = (2*R+1)^2;
        connections = si.bus(nNeighbors);
    };

//RICORDATI CHE PUOI USARE subseq PER ESTRARRE ELEMENTI DA UNA LISTA

buildScheme2D(pointsX,pointsY,coeffs) = par (x, pointsX,
                            par(x,pointsY, schemePoint))
    with
    {

    };

    process = _<:schemePoint2D(r,midCoeff,10);
