import("stdfaust.lib");

//--------------------------------Model Settings-----------------------------//
k = 1/ma.SR;
c = 344;
h = c * k*sqrt(2);
nPointsX = 30;
nPointsY = 10;

lambda = c*k/h;

//----------------------------------Equations--------------------------------//
schemeMidPoint(fIn,u_n,u_s,u,u_w,u_e) = u*2*(1-2*lambda*lambda) - u' + lambda*lambda*(u_e+u_w+u_n+u_s) + fIn;

buildScheme(X,Y) = par (x, X, 
                par(y,Y, schemeMidPoint));

//----------------------------------Controls---------------------------------//
inPointX=hslider("input point x", floor(nPointsX/2),0,nPointsX-1,0.01);
inPointY=hslider("input point y", floor(nPointsY/2),0,nPointsY-1,0.01);
outPointX=hslider("output point x",floor(nPointsX/2),0,nPointsX-1,0.01);
outPointY=hslider("output point y",floor(nPointsY/2),0,nPointsY-1,0.01);
forceModel = button("play") : ba.impulsify;
stop = button("Stop");

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
nInputs = inputs(schemeMidPoint);
routing(X, Y, x0, y0) = route(X*Y*2, X*Y*nInputs, 
                                par(x, X, par(y, Y, connections(x,y))))   
with 
{
    connections(x,y) =  P(x,y) + X*Y, F(x,y),
                        P(x,y), S(x,y-1),
                        P(x,y), N(x,y+1),
                        P(x,y), C(x,y),
                        P(x,y), E(x-1,y),
                        P(x,y), W(x+1,y);
    P(x,y) = x*Y+y+1;

    F(x,y) = (1 + 0 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
    N(x,y) = (1 + 1 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
    S(x,y) = (1 + 2 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
    C(x,y) = (1 + 3 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
    W(x,y) = (1 + 4 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
    E(x,y) = (1 + 5 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
};

model(X,Y) = (routing(X,Y, 2,2) : buildScheme(X,Y)) ~ par(i,X*Y,_*(stop==0));
    

process = forceModel <: linInterpolation2D(inPointX,inPointY) : model(nPointsX,nPointsY) : linInterpolation2DOut(outPointX,outPointY)<:_,_;
