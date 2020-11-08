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

midCoeffDelay1 = 0,0,0,
                0,-1,0,
                0,0,0;

r=1;
t=1;

//coefficients = midCoeff,midCoeffDelay1;

scheme(pointsX,pointsY) = par (nPointsX, pointsX,
                                par(nPointsY,pointsY, midCoeff));

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

buildScheme2D(R,T,pointsX,pointsY,coeffs,signals) =
    par (x, pointsX,
        par(y,pointsY, schemePoint2D(R,T,coeff(x,y),fIn(x,y),neighbors(x,y))))
    with
    {
        nNeighbors = (2*R+1)^2;
        coeff(x,y) = ba.subseq(coeffs,int((x*pointsY+y)*coeffLength),coeffLength);
        coeffLength = int(nNeighbors*(T+1));
        signalsLength = coeffLength+1;
        fIn(x,y) = ba.take(int((x*pointsY+y)*signalsLength),signals);
        neighbors(x,y) = ba.subseq(signals,int((x*pointsY+y)*signalsLength),signalsLength);
    };

process = 10,par(i,(2*r+1)^2,i):schemePoint2D(r,t,midCoeff);
//process = buildScheme2D(r,t,nPointsX,nPointsY,scheme(nPointsX,nPointsY),par(i,10*4,i));
//process = scheme(nPointsX,nPointsY);


/*//METHOD 2
schemePoint2D(c0,fIn,u0) = c0*u0 + fIn;
schemePoint2D((c0,cn),fIn,(u0,un)) = c0*u0 + schemePoint2D(cn,fIn,un);
process = _<:schemePoint2D(midCoeff,10);*/
//questo invece dà proprio errore, sembra che non arrivi alla seconda definizione
